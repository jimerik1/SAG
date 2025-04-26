// src/calculators/fea_core/fea_solver.cpp
#include "fea_solver.h"
#include <cmath>
#include <Eigen/SparseLU>

FEASolver::FEASolver(const std::vector<BHAComponent>& components,
                    const std::vector<MaterialProperties>& materials,
                    const Trajectory& trajectory,
                    double elementSize,
                    double mudWeight,
                    double wob)
    : components_(components),
      materials_(materials),
      trajectory_(trajectory),
      elementSize_(elementSize),
      mudWeight_(mudWeight),
      wob_(wob),
      numNodes_(0),
      numDofs_(0) {}

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

void FEASolver::solve() {
    // Solve the system Kd = F
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(globalStiffnessMatrix_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Factorization failed");
    }
    
    displacements_ = solver.solve(globalLoadVector_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }
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

void FEASolver::solve() {
    // Solve the system Kd = F
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(globalStiffnessMatrix_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Factorization failed");
    }
    
    displacements_ = solver.solve(globalLoadVector_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }
    
    // Calculate internal forces and side forces after solving
    calculateInternalForces();
    calculateSideForces();
}
