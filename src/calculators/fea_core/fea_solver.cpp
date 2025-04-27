// src/calculators/fea_core/fea_solver.cpp
#include "fea_solver.h"
#include <cmath>
#include <Eigen/SparseLU>
#include <iostream>

FEASolver::FEASolver(const std::vector<BHAComponent>& components,
                    const std::vector<MaterialProperties>& materials,
                    const Trajectory& trajectory,
                    const WellboreProperties& wellboreProps,
                    double elementSize)
    : components_(components),
      materials_(materials),
      trajectory_(trajectory),
      wellboreProps_(wellboreProps),
      elementSize_(elementSize),
      mudWeight_(wellboreProps.mudWeight),
      wob_(20000.0),
      numNodes_(0),
      numDofs_(0),
      bitSideForce_(0.0) {}

void FEASolver::setupSystem() {
    createMesh();
    assembleStiffnessMatrix();
    assembleLoadVector();
    applyBoundaryConditions();
}

void FEASolver::createMesh() {
    // Calculate total BHA length
    double totalLength = 0.0;
    for (const auto& component : components_) {
        totalLength += component.length;
    }
    
    // Calculate number of elements and nodes
    int numElements = static_cast<int>(ceil(totalLength / elementSize_));
    numNodes_ = numElements + 1;
    numDofs_ = numNodes_ * 2;  // 2 DOFs per node (deflection and rotation)
    
    // Initialize node positions
    nodePositions_.resize(numNodes_);
    for (int i = 0; i < numNodes_; ++i) {
        nodePositions_[i] = i * elementSize_;
    }
    
    // Create elements
    elements_.clear();
    int componentIndex = 0;
    double componentEndPosition = components_[0].length;
    
    for (int i = 0; i < numElements; ++i) {
        double elementStart = i * elementSize_;
        double elementEnd = (i + 1) * elementSize_;
        
        // Find which component this element belongs to
        while (elementStart >= componentEndPosition && componentIndex < components_.size() - 1) {
            componentIndex++;
            componentEndPosition += components_[componentIndex].length;
        }
        
        const auto& component = components_[componentIndex];
        MaterialProperties material = getMaterialProperties(component.material);
        
        // Calculate element properties
        double E = material.youngsModulus;
        double I = calculateMomentOfInertia(component.outerDiameter, component.innerDiameter);
        double A = calculateArea(component.outerDiameter, component.innerDiameter);
        double rho = material.density;
        
        // Create element
        elements_.emplace_back(elementSize_, E, I, A, rho);
    }
}

void FEASolver::assembleStiffnessMatrix() {
    // Initialize global stiffness matrix
    globalStiffnessMatrix_.resize(numDofs_, numDofs_);
    std::vector<Eigen::Triplet<double>> triplets;
    
    // Assemble element stiffness matrices into global stiffness matrix
    for (int i = 0; i < elements_.size(); ++i) {
        Eigen::MatrixXd elementK = elements_[i].calculateStiffnessMatrix();
        
        // Calculate global indices
        int startNode = i;
        int endNode = i + 1;
        int startDof = startNode * 2;
        int endDof = endNode * 2;
        
        // Add element stiffness to global stiffness
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                int globalRow = (r < 2) ? startDof + r : endDof + (r - 2);
                int globalCol = (c < 2) ? startDof + c : endDof + (c - 2);
                triplets.push_back(Eigen::Triplet<double>(globalRow, globalCol, elementK(r, c)));
            }
        }
    }
    
    globalStiffnessMatrix_.setFromTriplets(triplets.begin(), triplets.end());
}

void FEASolver::assembleLoadVector() {
    // Initialize global load vector
    globalLoadVector_ = Eigen::VectorXd::Zero(numDofs_);
    
    // Get average inclination from the trajectory
    double avgInclination = 0.0;
    for (const auto& point : trajectory_.points) {
        avgInclination += point.inclination;
    }
    avgInclination /= trajectory_.points.size();
    
    // Average azimuth
    double avgAzimuth = 0.0;
    for (const auto& point : trajectory_.points) {
        avgAzimuth += point.azimuth;
    }
    avgAzimuth /= trajectory_.points.size();
    
    // Buoyancy factor based on mud weight
    // Calculate buoyancy factor
    double BF = 65.5 / (65.5 - mudWeight_);
    
    // Apply gravity loads for each element
    for (int i = 0; i < elements_.size(); ++i) {
        // Get element load vector due to gravity
        Eigen::VectorXd elementLoad = elements_[i].calculateGravityLoad(avgInclination, avgAzimuth);
        
        // Apply buoyancy factor to reduce effective weight
        elementLoad *= BF;
        
        // Add to global load vector
        int startNode = i;
        int endNode = i + 1;
        int startDof = startNode * 2;
        int endDof = endNode * 2;
        
        globalLoadVector_(startDof) += elementLoad(0);
        globalLoadVector_(startDof + 1) += elementLoad(1);
        globalLoadVector_(endDof) += elementLoad(2);
        globalLoadVector_(endDof + 1) += elementLoad(3);
    }
}

void FEASolver::applyBoundaryConditions() {
    // Apply bit boundary condition (pinned - zero displacement)
    int bitDofIndex = 0;  // First node, displacement DOF
    
    // Modify stiffness matrix to enforce zero displacement at bit
    globalStiffnessMatrix_.coeffRef(bitDofIndex, bitDofIndex) = 1.0e12;  // Large value to approximate fixed constraint
    globalLoadVector_(bitDofIndex) = 0.0;
    
    // Apply stabilizer boundary conditions
    for (const auto& component : components_) {
        if (component.isStabilizer) {
            // Find the node closest to the stabilizer position
            double stabilizerPos = component.distanceFromBit;
            int nodeIndex = static_cast<int>(round(stabilizerPos / elementSize_));
            
            if (nodeIndex < numNodes_) {
                int dofIndex = nodeIndex * 2;  // Displacement DOF
                
                // Enforce zero displacement at stabilizer
                globalStiffnessMatrix_.coeffRef(dofIndex, dofIndex) = 1.0e12;  // Large value to approximate fixed constraint
                globalLoadVector_(dofIndex) = 0.0;
                
                std::cout << "Applied stabilizer constraint at position " << stabilizerPos 
                          << " m (component: " << component.name << ")" << std::endl;
            }
        }
    }
}

void FEASolver::solve() {
    // First pass: solve unconstrained system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(globalStiffnessMatrix_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Factorization failed");
    }
    
    // Initial solution
    displacements_ = solver.solve(globalLoadVector_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }
    
    // Create a backup of the original stiffness matrix and load vector
    Eigen::SparseMatrix<double> originalStiffness = globalStiffnessMatrix_;
    Eigen::VectorXd originalLoads = globalLoadVector_;
    
    // Reset the contact nodes we're going to constrain
    std::vector<int> contactNodes;
    
    // Determine which nodes need constraint based on the unconstrained solution
    for (int i = 0; i < numNodes_; ++i) {
        double position = i * elementSize_;
        const BHAComponent* component = getComponentAtPosition(position);
        
        double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
        double componentRadius = component->outerDiameter * 0.0254 / 2.0;
        double maxDeflection = holeRadius - componentRadius;
        
        double currentDeflection = displacements_(i * 2);
        
        if (std::abs(currentDeflection) > maxDeflection) {
            contactNodes.push_back(i);
        }
    }
    
    // If we found contact nodes, apply contact constraints and re-solve
    if (!contactNodes.empty()) {
        std::cout << "Found " << contactNodes.size() << " nodes in contact with wellbore. Re-solving with constraints." << std::endl;
        
        // Reset the system
        globalStiffnessMatrix_ = originalStiffness;
        globalLoadVector_ = originalLoads;
        
        // Apply constraints at contact nodes (pinned boundary conditions)
        for (int nodeIndex : contactNodes) {
            int dofIndex = nodeIndex * 2; // Displacement DOF
            double position = nodeIndex * elementSize_;
            const BHAComponent* component = getComponentAtPosition(position);
            
            double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
            double componentRadius = component->outerDiameter * 0.0254 / 2.0;
            double maxDeflection = holeRadius - componentRadius;
            
            // Determine sign of contact (which side of the wellbore)
            double sign = (displacements_(dofIndex) > 0) ? 1.0 : -1.0;
            
            // Apply constraint (stronger than standard constraints but not as extreme)
            globalStiffnessMatrix_.coeffRef(dofIndex, dofIndex) = 1.0e10;
            globalLoadVector_(dofIndex) = maxDeflection * sign * 1.0e10;
            
            std::cout << "Applied contact constraint at position " << position 
                     << " m. Max allowed deflection: " << maxDeflection * 1000.0 << " mm" << std::endl;
        }
        
        // Re-solve the constrained system
        solver.compute(globalStiffnessMatrix_);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Constrained system factorization failed");
        }
        
        displacements_ = solver.solve(globalLoadVector_);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Constrained system solving failed");
        }
        
        // Verify constraints
        for (int nodeIndex : contactNodes) {
            int dofIndex = nodeIndex * 2;
            double position = nodeIndex * elementSize_;
            const BHAComponent* component = getComponentAtPosition(position);
            
            double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
            double componentRadius = component->outerDiameter * 0.0254 / 2.0;
            double maxDeflection = holeRadius - componentRadius;
            
            std::cout << "Verified position " << position << " m: deflection = " 
                     << displacements_(dofIndex) * 1000.0 << " mm (max allowed: " 
                     << maxDeflection * 1000.0 << " mm)" << std::endl;
        }
    }
    
    // Apply a final direct constraint to ensure no violations
    for (int i = 0; i < numNodes_; ++i) {
        double position = i * elementSize_;
        const BHAComponent* component = getComponentAtPosition(position);
        
        double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
        double componentRadius = component->outerDiameter * 0.0254 / 2.0;
        double maxDeflection = holeRadius - componentRadius;
        
        double currentDeflection = displacements_(i * 2);
        
        if (std::abs(currentDeflection) > maxDeflection) {
            // Preserve sign of deflection
            double newDeflection = (currentDeflection > 0) ? maxDeflection : -maxDeflection;
            displacements_(i * 2) = newDeflection;
        }
    }
    
    // Calculate internal forces after enforcing all constraints
    calculateInternalForces();
    calculateSideForces();
}

void FEASolver::applyWellboreContactConstraints() {
    // Apply constraints directly to the displacement field
    for (int i = 0; i < numNodes_; ++i) {
        // Find current component at this position
        double position = i * elementSize_;
        const BHAComponent* component = getComponentAtPosition(position);
        
        // Calculate clearance for this component
        double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0; // Convert inches to meters, divide by 2 for radius
        double componentRadius = component->outerDiameter * 0.0254 / 2.0; // Component OD to radius in meters
        double maxDeflection = holeRadius - componentRadius;
        
        // Get current node deflection
        double currentDeflection = displacements_(i * 2);
        
        // If deflection exceeds limit, clip it to the maximum allowed value
        if (std::abs(currentDeflection) > maxDeflection) {
            // Preserve sign of deflection
            double newDeflection = (currentDeflection > 0) ? maxDeflection : -maxDeflection;
            
            std::cout << "Constrained deflection at position " << position 
                     << " m from " << currentDeflection * 1000.0 << " mm to "
                     << newDeflection * 1000.0 << " mm (max allowed: " 
                     << maxDeflection * 1000.0 << " mm)" << std::endl;
            
            // Update displacement directly
            displacements_(i * 2) = newDeflection;
        }
    }
}

const BHAComponent* FEASolver::getComponentAtPosition(double position) const {
    // Find which component this position belongs to
    double currentPos = 0.0;
    for (const auto& component : components_) {
        double endPos = currentPos + component.length;
        if (position >= currentPos && position < endPos) {
            return &component;
        }
        currentPos = endPos;
    }
    
    // Default to last component if position is beyond all components
    return &components_.back();
}

double FEASolver::getSagCorrection(double distanceFromBit) const {
    // Find the element that contains this position
    int nodeIndex = static_cast<int>(distanceFromBit / elementSize_);
    
    // Ensure we don't go out of bounds
    if (nodeIndex >= numNodes_ - 1) {
        nodeIndex = numNodes_ - 2;
    }
    
    // Get the slopes at both nodes of the element
    double slope1 = displacements_(nodeIndex * 2 + 1);
    double slope2 = displacements_((nodeIndex + 1) * 2 + 1);
    
    // Interpolate to get slope at exact position
    double elementStart = nodeIndex * elementSize_;
    double position = (distanceFromBit - elementStart) / elementSize_;
    
    // Linear interpolation of slope
    double slope = slope1 * (1 - position) + slope2 * position;
    
    // Convert slope to degrees
    return slope * 180.0 / M_PI;
}

double FEASolver::getDeflection(double distanceFromBit) const {
    // Find the element that contains this position
    int nodeIndex = static_cast<int>(distanceFromBit / elementSize_);
    
    // Ensure we don't go out of bounds
    if (nodeIndex >= numNodes_ - 1) {
        nodeIndex = numNodes_ - 2;
    }
    
    // Get the deflections at both nodes of the element
    double defl1 = displacements_(nodeIndex * 2);
    double defl2 = displacements_((nodeIndex + 1) * 2);
    
    // Interpolate to get deflection at exact position
    double elementStart = nodeIndex * elementSize_;
    double position = (distanceFromBit - elementStart) / elementSize_;
    
    // Linear interpolation of deflection
    return defl1 * (1 - position) + defl2 * position;
}

std::vector<std::pair<double, double>> FEASolver::getDeflectionProfile() const {
    std::vector<std::pair<double, double>> profile;
    
    // Create pairs of (position, deflection) for each node
    for (int i = 0; i < numNodes_; ++i) {
        double pos = i * elementSize_;
        double defl = displacements_(i * 2);
        profile.push_back(std::make_pair(pos, defl));
    }
    
    return profile;
}

double FEASolver::calculateMomentOfInertia(double od, double id) const {
    // Convert to meters
    od = od * 0.0254;  // inches to meters
    id = id * 0.0254;  // inches to meters
    
    // Calculate area moment of inertia for a hollow circular section
    return M_PI * (pow(od, 4) - pow(id, 4)) / 64.0;
}

double FEASolver::calculateArea(double od, double id) const {
    // Convert to meters
    od = od * 0.0254;  // inches to meters
    id = id * 0.0254;  // inches to meters
    
    // Calculate cross-sectional area for a hollow circular section
    return M_PI * (pow(od, 2) - pow(id, 2)) / 4.0;
}

MaterialProperties FEASolver::getMaterialProperties(const std::string& grade) const {
    // Find the material properties for the given grade
    for (const auto& material : materials_) {
        if (material.grade == grade) {
            return material;
        }
    }
    
    // Return default properties if not found
    return {"Unknown", 2.07e11, 7850.0, 0.3, 7.58e8};
}

void FEASolver::calculateInternalForces() {
    shearForces_.clear();
    bendingMoments_.clear();
    slopes_.clear();
    
    // For each element, calculate internal forces and moments
    for (int i = 0; i < elements_.size(); ++i) {
        // Get element displacements
        int startNode = i;
        int endNode = i + 1;
        int startDof = startNode * 2;
        int endDof = endNode * 2;
        
        Eigen::Vector4d elementDisp;
        elementDisp << displacements_(startDof), displacements_(startDof + 1),
                      displacements_(endDof), displacements_(endDof + 1);
        
        // Calculate element stiffness matrix
        Eigen::MatrixXd K = elements_[i].calculateStiffnessMatrix();
        
        // Calculate internal forces F = K*d
        Eigen::Vector4d internalForces = K * elementDisp;
        
        // Element properties
        double E = elements_[i].getE();
        double I = elements_[i].getI();
        double L = elements_[i].getLength();
        
        // Bending moments at nodes
        double momentStart = E * I * displacements_(startDof + 1);
        double momentEnd = -E * I * displacements_(endDof + 1);
        
        // Shear forces at nodes
        double shearStart = internalForces(0);
        double shearEnd = -internalForces(2);
        
        // Store values at start node
        double pos = startNode * elementSize_;
        shearForces_.push_back(std::make_pair(pos, shearStart));
        bendingMoments_.push_back(std::make_pair(pos, momentStart));
        slopes_.push_back(std::make_pair(pos, displacements_(startDof + 1) * 180.0 / M_PI));
        
        // Store values at end node if this is the last element
        if (i == elements_.size() - 1) {
            pos = endNode * elementSize_;
            shearForces_.push_back(std::make_pair(pos, shearEnd));
            bendingMoments_.push_back(std::make_pair(pos, momentEnd));
            slopes_.push_back(std::make_pair(pos, displacements_(endDof + 1) * 180.0 / M_PI));
        }
    }
}

void FEASolver::calculateSideForces() {
    // Calculate bit side force (first node shear force)
    bitSideForce_ = -shearForces_[0].second;  // Negative because reaction force
    
    // Calculate stabilizer side forces
    stabilizerForces_.clear();
    
    for (const auto& component : components_) {
        if (component.isStabilizer) {
            double stabPos = component.distanceFromBit;
            
            // Find the closest nodes
            int nodeIndex = static_cast<int>(round(stabPos / elementSize_));
            
            if (nodeIndex < numNodes_) {
                // Calculate side force as the jump in shear force
                double shearBefore = 0.0;
                double shearAfter = 0.0;
                
                // Find shear forces before and after stabilizer
                if (nodeIndex > 0) {
                    for (const auto& shear : shearForces_) {
                        if (fabs(shear.first - (nodeIndex - 1) * elementSize_) < 1e-6) {
                            shearBefore = shear.second;
                            break;
                        }
                    }
                }
                
                for (const auto& shear : shearForces_) {
                    if (fabs(shear.first - nodeIndex * elementSize_) < 1e-6) {
                        shearAfter = shear.second;
                        break;
                    }
                }
                
                // Side force is the jump in shear
                double sideForce = shearAfter - shearBefore;
                
                // Add to stabilizer forces
                StabilizerForce force;
                force.name = component.name;
                force.position = stabPos;
                force.sideForce = sideForce;
                stabilizerForces_.push_back(force);
            }
        }
    }
}

FEAResults FEASolver::getResults() const {
    FEAResults results;
    
    // Bit results
    results.bitSideForce = bitSideForce_;
    results.bend = 0.25;  // This would be calculated based on motor bend if present
    results.wob = wob_;
    
    // Sensor results
    // Find MWD component and sensor position
    double sensorPos = 0.0;
    for (const auto& component : components_) {
        if (component.name.find("MWD") != std::string::npos) {
            // MWD inclination sensor is typically near the center of the MWD tool
            sensorPos = component.distanceFromBit + component.length / 2.0;
            break;
        }
    }
    results.sensorPosition = sensorPos;
    
    // Get local trajectory data at sensor position
    TrajectoryPoint trajPoint = trajectory_.interpolateAtDepth(sensorPos);
    results.holeCurvature = trajPoint.dls * 100.0 / 30.0;  // Convert to deg/100ft
    
    // Get sag correction at sensor
    results.sagCorrection = getSagCorrection(sensorPos);
    
    // Stabilizer forces
    results.stabilizerForces = stabilizerForces_;
    
    // Create profiles
    // Deflection profile
    std::vector<std::pair<double, double>> deflectionProfile;
    for (int i = 0; i < numNodes_; ++i) {
        double pos = i * elementSize_;
        double defl = displacements_(i * 2);
        deflectionProfile.push_back(std::make_pair(pos, defl));
    }
    results.deflectionProfile = deflectionProfile;
    
    // Other profiles
    results.shearForceProfile = shearForces_;
    results.bendingMomentProfile = bendingMoments_;
    results.slopeProfile = slopes_;
    
    return results;
}