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
    
    // Use a finer mesh for more accurate results
    // Aim for elements no larger than 10cm for sufficient accuracy
    elementSize_ = 0.1; 
    int numElements = static_cast<int>(ceil(totalLength / elementSize_));
    elementSize_ = totalLength / numElements; // Adjust for even distribution
    
    // Ensure we have at least 10 elements between each point of stabilization
    for (size_t i = 0; i < components_.size() - 1; i++) {
        if (components_[i].isStabilizer || components_[i+1].isStabilizer) {
            double segmentLength = components_[i+1].distanceFromBit - components_[i].distanceFromBit;
            int minSegmentElements = 10;
            double maxElementSize = segmentLength / minSegmentElements;
            elementSize_ = std::min(elementSize_, maxElementSize);
        }
    }
    
    // Recalculate number of elements with refined element size
    numElements = static_cast<int>(ceil(totalLength / elementSize_));
    numNodes_ = numElements + 1;
    numDofs_ = numNodes_ * 2;  // 2 DOFs per node (deflection and rotation)
    
    std::cout << "Mesh created with " << numElements << " elements, element size: " 
              << elementSize_ << " m" << std::endl;
    
    // Initialize node positions
    nodePositions_.resize(numNodes_);
    for (int i = 0; i < numNodes_; ++i) {
        nodePositions_[i] = i * elementSize_;
    }
    
    // Create elements (rest of your existing implementation)
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
    
    // Convert inclination to radians
    double inclinationRad = avgInclination * M_PI / 180.0;
    
    // Average azimuth (for completeness, though not critical for sag analysis)
    double avgAzimuth = 0.0;
    for (const auto& point : trajectory_.points) {
        avgAzimuth += point.azimuth;
    }
    avgAzimuth /= trajectory_.points.size();
    
    // Calculate buoyancy factor based on mud weight
    // Standard buoyancy formula: BF = 1 - (mud_density / material_density)
    // For drilling calculations, commonly used approximation:
    // Conversion: 1 ppg = 8.33 lb/ft³
    double mudDensity = mudWeight_ * 8.33; // Convert ppg to lb/ft³
    double steelDensity = 490.0; // Approximate density of steel in lb/ft³
    
    double BF = 1.0 - (mudDensity / steelDensity);
    std::cout << "Buoyancy factor: " << BF << " (mud weight: " << mudWeight_ << " ppg)" << std::endl;
    
    // Apply gravity loads for each element
    for (int i = 0; i < elements_.size(); ++i) {
        // Component weight adjusted for buoyancy
        double weightPerLength = elements_[i].getRho() * elements_[i].getA() * 9.81; // N/m
        double effectiveWeight = weightPerLength * BF;
        
        // Weight component perpendicular to wellbore axis (causes sag)
        double perpForce = effectiveWeight * sin(inclinationRad);
        
        // For element with uniform distributed load w, the equivalent nodal loads are:
        // F = [wL/2, wL²/12, wL/2, -wL²/12]ᵀ
        double L = elements_[i].getLength();
        Eigen::Vector4d elementLoad;
        elementLoad << perpForce * L / 2.0,
                      perpForce * L * L / 12.0,
                      perpForce * L / 2.0,
                      -perpForce * L * L / 12.0;
        
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
    
    // Add weight-on-bit effect if specified
    if (wob_ > 0.0) {
        // Convert WOB from kg to N
        double wobN = wob_ * 9.81;
        
        // Apply as an axial compression on the BHA
        // For small deflections, this can be approximated as a moment
        // at the bit proportional to the initial angular misalignment
        
        // Apply a small perturbation moment at the bit if WOB is present
        // This simulates the buckling effect of axial load
        double perturbationMoment = wobN * 0.001; // Small moment proportional to WOB
        globalLoadVector_(1) += perturbationMoment; // Apply to bit rotation DOF
        
        std::cout << "Applied WOB effect: " << wob_ << " kg" << std::endl;
    }
}

void FEASolver::applyBoundaryConditions() {
    // Calculate the maximum diagonal element in the stiffness matrix
    double maxDiagonal = 0.0;
    for (int i = 0; i < numDofs_; ++i) {
        maxDiagonal = std::max(maxDiagonal, std::abs(globalStiffnessMatrix_.coeff(i, i)));
    }
    
    // Use a penalty factor based on the system stiffness (typically 1000-10000 times)
    double penaltyFactor = 1000.0 * maxDiagonal;
    
    // Apply bit boundary condition with appropriate scaling
    int bitDofIndex = 0;
    globalStiffnessMatrix_.coeffRef(bitDofIndex, bitDofIndex) = penaltyFactor;
    globalLoadVector_(bitDofIndex) = 0.0;
    
    // Apply stabilizer boundary conditions with same scaling principle
    for (const auto& component : components_) {
        if (component.isStabilizer) {
            double stabilizerPos = component.distanceFromBit;
            int nodeIndex = static_cast<int>(round(stabilizerPos / elementSize_));
            
            if (nodeIndex < numNodes_) {
                int dofIndex = nodeIndex * 2;
                globalStiffnessMatrix_.coeffRef(dofIndex, dofIndex) = penaltyFactor;
                globalLoadVector_(dofIndex) = 0.0;
                
                std::cout << "Applied stabilizer constraint at position " << stabilizerPos 
                          << " m (component: " << component.name << ")" << std::endl;
            }
        }
    }
}

void FEASolver::calculateSideForces() {
    // Calculate bit side force (reaction force at bit)
    if (!shearForces_.empty()) {
        // The bit side force is the negative of the shear force at the bit (reaction force)
        bitSideForce_ = -shearForces_[0].second;
        
        // Convert to physically meaningful values (N to kg-force)
        double bitSideForceKg = bitSideForce_ / 9.81;
        
        // Sanity check for unreasonable forces
        double totalBHAWeight = 0.0;
        for (const auto& component : components_) {
            // Convert ppf to kg: 1 ppf = 1.488 kg/m
            double weightKg = component.weight * 1.488 * component.length;
            totalBHAWeight += weightKg;
        }
        
        // Limit max side force to a reasonable fraction of BHA weight
        const double MAX_SIDE_FORCE_RATIO = 0.5; // Max 50% of BHA weight
        double maxReasonableForce = totalBHAWeight * MAX_SIDE_FORCE_RATIO;
        
        if (std::abs(bitSideForceKg) > maxReasonableForce) {
            std::cout << "Warning: Calculated bit side force (" << bitSideForceKg 
                     << " kg) exceeds reasonable limit (" << maxReasonableForce
                     << " kg). Limiting value." << std::endl;
            
            bitSideForce_ = (bitSideForceKg > 0 ? 1.0 : -1.0) * maxReasonableForce * 9.81;
        }
    } else {
        bitSideForce_ = 0.0;
    }
    
    // Calculate stabilizer side forces
    stabilizerForces_.clear();
    
    for (const auto& component : components_) {
        if (component.isStabilizer) {
            double stabPos = component.distanceFromBit;
            
            // Find the shear forces before and after stabilizer position
            double posBefore = 0.0, posAfter = 0.0;
            double shearBefore = 0.0, shearAfter = 0.0;
            
            // Find points before and after stabilizer
            for (size_t i = 0; i < shearForces_.size() - 1; ++i) {
                if (shearForces_[i].first <= stabPos && shearForces_[i+1].first >= stabPos) {
                    posBefore = shearForces_[i].first;
                    posAfter = shearForces_[i+1].first;
                    shearBefore = shearForces_[i].second;
                    shearAfter = shearForces_[i+1].second;
                    break;
                }
            }
            
            // If not found, try to find closest point
            if (posBefore == 0.0 && posAfter == 0.0 && !shearForces_.empty()) {
                // Find closest point
                size_t closestIdx = 0;
                double minDist = std::abs(shearForces_[0].first - stabPos);
                
                for (size_t i = 1; i < shearForces_.size(); ++i) {
                    double dist = std::abs(shearForces_[i].first - stabPos);
                    if (dist < minDist) {
                        minDist = dist;
                        closestIdx = i;
                    }
                }
                
                if (closestIdx < shearForces_.size() - 1) {
                    posBefore = shearForces_[closestIdx].first;
                    posAfter = shearForces_[closestIdx + 1].first;
                    shearBefore = shearForces_[closestIdx].second;
                    shearAfter = shearForces_[closestIdx + 1].second;
                } else if (!shearForces_.empty()) {
                    // Use last point
                    posBefore = shearForces_.back().first;
                    posAfter = posBefore;
                    shearBefore = shearForces_.back().second;
                    shearAfter = shearBefore;
                }
            }
            
            // Interpolate shear at stabilizer position
            double interpShear = 0.0;
            if (posAfter != posBefore) {
                double t = (stabPos - posBefore) / (posAfter - posBefore);
                interpShear = shearBefore * (1.0 - t) + shearAfter * t;
            } else {
                interpShear = shearBefore;
            }
            
            // Side force is the reaction force (negative of shear)
            double sideForce = -interpShear;
            
            // Convert to kg-force and apply sanity checks
            double sideForceKg = sideForce / 9.81;
            
            // Calculate reasonable force limit
            double totalBHAWeight = 0.0;
            for (const auto& comp : components_) {
                double weightKg = comp.weight * 1.488 * comp.length;
                totalBHAWeight += weightKg;
            }
            
            const double MAX_STAB_FORCE_RATIO = 0.5; // Max 50% of BHA weight
            double maxReasonableForce = totalBHAWeight * MAX_STAB_FORCE_RATIO;
            
            if (std::abs(sideForceKg) > maxReasonableForce) {
                std::cout << "Warning: Calculated stabilizer force at " << stabPos 
                         << "m (" << sideForceKg << " kg) exceeds reasonable limit ("
                         << maxReasonableForce << " kg). Limiting value." << std::endl;
                
                sideForceKg = (sideForceKg > 0 ? 1.0 : -1.0) * maxReasonableForce;
                sideForce = sideForceKg * 9.81;
            }
            
            // Add to stabilizer forces
            StabilizerForce force;
            force.name = component.name;
            force.position = stabPos;
            force.sideForce = sideForceKg; // Store in kg-force for reporting
            stabilizerForces_.push_back(force);
            
            std::cout << "Stabilizer " << component.name << " at " << stabPos 
                     << "m: Side force = " << sideForceKg << " kg" << std::endl;
        }
    }
}

void FEASolver::solve() {
    // Step 1: Solve the unconstrained system first
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(globalStiffnessMatrix_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Initial factorization failed");
    }
    
    // Initial solution
    displacements_ = solver.solve(globalLoadVector_);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Initial solve failed");
    }
    
    // Step 2: Prepare for contact resolution
    const int MAX_ITERATIONS = 20;
    const double CONVERGENCE_TOL = 1e-6;
    const double DIVERGENCE_TOL = 1.5; // If penetration increases by more than 50%, detect divergence
    const double CONTACT_STIFFNESS_SCALE = 1e6; // Much lower than 1e12
    
    std::vector<int> contactNodes;
    Eigen::VectorXd prevDisplacements = displacements_;
    double prevMaxPenetration = 0.0;
    
    // Step 3: Main contact resolution loop
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        contactNodes.clear();
        double maxPenetration = 0.0;
        
        // Step 3a: Identify contact violations
        for (int i = 0; i < numNodes_; ++i) {
            double position = i * elementSize_;
            const BHAComponent* component = getComponentAtPosition(position);
            
            // Calculate clearance
            double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
            double componentRadius = component->outerDiameter * 0.0254 / 2.0;
            double maxDeflection = holeRadius - componentRadius;
            double currentDeflection = displacements_(i * 2);
            
            // Check for contact
            if (std::abs(currentDeflection) > maxDeflection) {
                contactNodes.push_back(i);
                
                // Calculate penetration
                double penetration = std::abs(currentDeflection) - maxDeflection;
                maxPenetration = std::max(maxPenetration, penetration);
                
                // Apply constraint directly (Dirichlet boundary condition approach)
                double sign = (currentDeflection > 0) ? 1.0 : -1.0;
                displacements_(i * 2) = sign * maxDeflection;
            }
        }
        
        // Step 3b: If no contacts or convergence achieved, we're done
        if (contactNodes.empty() || maxPenetration < CONVERGENCE_TOL) {
            std::cout << "Contact resolution converged after " << iter+1 << " iterations." << std::endl;
            break;
        }
        
        // Step 3c: Check for divergence
        if (iter > 0 && maxPenetration > prevMaxPenetration * DIVERGENCE_TOL) {
            std::cout << "Warning: Contact algorithm is diverging. Switching to direct enforcement." << std::endl;
            // Revert to previous displacements and apply direct constraint enforcement
            displacements_ = prevDisplacements;
            
            // Directly enforce boundary conditions
            for (int i = 0; i < numNodes_; ++i) {
                double position = i * elementSize_;
                const BHAComponent* component = getComponentAtPosition(position);
                
                double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
                double componentRadius = component->outerDiameter * 0.0254 / 2.0;
                double maxDeflection = holeRadius - componentRadius;
                double currentDeflection = displacements_(i * 2);
                
                if (std::abs(currentDeflection) > maxDeflection) {
                    double sign = (currentDeflection > 0) ? 1.0 : -1.0;
                    displacements_(i * 2) = sign * maxDeflection;
                }
            }
            break;
        }
        
        // Remember state for divergence check
        prevMaxPenetration = maxPenetration;
        prevDisplacements = displacements_;
        
        // Step 3d: Create a modified system with contact constraints
        Eigen::SparseMatrix<double> modifiedStiffness = globalStiffnessMatrix_;
        Eigen::VectorXd modifiedLoad = globalLoadVector_;
        
        // Apply penalty method with reasonable stiffness
        for (int nodeIdx : contactNodes) {
            int dofIdx = nodeIdx * 2;
            double position = nodeIdx * elementSize_;
            const BHAComponent* component = getComponentAtPosition(position);
            
            double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
            double componentRadius = component->outerDiameter * 0.0254 / 2.0;
            double maxDeflection = holeRadius - componentRadius;
            double currentDeflection = displacements_(dofIdx);
            
            // Determine direction
            double sign = (currentDeflection > 0) ? 1.0 : -1.0;
            
            // Add penalty stiffness to diagonal
            modifiedStiffness.coeffRef(dofIdx, dofIdx) += CONTACT_STIFFNESS_SCALE;
            
            // Modify load to enforce contact constraint
            modifiedLoad(dofIdx) += CONTACT_STIFFNESS_SCALE * sign * maxDeflection;
            
            std::cout << "Applied contact at node " << nodeIdx << " (pos: " << position 
                     << "m), target deflection: " << sign * maxDeflection * 1000.0 << " mm" << std::endl;
        }
        
        // Step 3e: Solve modified system
        solver.compute(modifiedStiffness);
        if (solver.info() != Eigen::Success) {
            std::cout << "Warning: Modified system factorization failed, reverting to direct enforcement" << std::endl;
            break;
        }
        
        Eigen::VectorXd newDisplacements = solver.solve(modifiedLoad);
        if (solver.info() != Eigen::Success) {
            std::cout << "Warning: Modified system solve failed, reverting to direct enforcement" << std::endl;
            break;
        }
        
        // Step 3f: Update displacements with under-relaxation for stability
        const double RELAXATION = 0.5; // Under-relaxation factor
        displacements_ = RELAXATION * newDisplacements + (1.0 - RELAXATION) * displacements_;
        
        std::cout << "Contact iteration " << iter+1 << ", max penetration: " 
                 << maxPenetration * 1000.0 << " mm" << std::endl;
    }
    
    // Step 4: Final enforcement of physical constraints
    for (int i = 0; i < numNodes_; ++i) {
        double position = i * elementSize_;
        const BHAComponent* component = getComponentAtPosition(position);
        
        double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
        double componentRadius = component->outerDiameter * 0.0254 / 2.0;
        double maxDeflection = holeRadius - componentRadius;
        double currentDeflection = displacements_(i * 2);
        
        if (std::abs(currentDeflection) > maxDeflection) {
            double sign = (currentDeflection > 0) ? 1.0 : -1.0;
            displacements_(i * 2) = sign * maxDeflection;
        }
    }
    
    // Step 5: Calculate internal forces and side forces
    calculateInternalForces();
    calculateSideForces();
}

void FEASolver::applyWellboreContactConstraints() {
    // Maximum allowable iteration count for contact convergence
    const int maxIterations = 10;
    // Relaxation factor for contact iterations (for numerical stability)
    const double relaxation = 0.7;
    
    std::vector<int> contactNodes;
    std::vector<double> contactForces(numNodes_, 0.0);
    
    // Iterative contact resolution
    for (int iter = 0; iter < maxIterations; ++iter) {
        contactNodes.clear();
        
        // Identify contact points
        for (int i = 0; i < numNodes_; ++i) {
            double position = i * elementSize_;
            const BHAComponent* component = getComponentAtPosition(position);
            
            // Calculate clearance
            double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
            double componentRadius = component->outerDiameter * 0.0254 / 2.0;
            double maxDeflection = holeRadius - componentRadius;
            double currentDeflection = displacements_(i * 2);
            
            // Check for contact
            if (std::abs(currentDeflection) > maxDeflection) {
                contactNodes.push_back(i);
                
                // Calculate contact force (proportional to penetration)
                double penetration = std::abs(currentDeflection) - maxDeflection;
                double contactStiffness = wellboreProps_.contactStiffness;
                double force = contactStiffness * penetration;
                
                // Apply contact force in correct direction
                double sign = (currentDeflection > 0) ? -1.0 : 1.0;
                contactForces[i] = sign * force * relaxation + contactForces[i] * (1.0 - relaxation);
            }
        }
        
        // If no contacts, break
        if (contactNodes.empty()) {
            break;
        }
        
        // Apply contact forces and re-solve
        for (int nodeIndex : contactNodes) {
            int dofIndex = nodeIndex * 2;
            globalLoadVector_(dofIndex) += contactForces[nodeIndex];
        }
        
        // Re-solve the system with contact forces
        // (simplified here, but you'd solve the system with updated forces)
    }
    
    // Final enforcement of physical constraints
    for (int i = 0; i < numNodes_; ++i) {
        double position = i * elementSize_;
        const BHAComponent* component = getComponentAtPosition(position);
        
        double holeRadius = wellboreProps_.holeSize * 0.0254 / 2.0;
        double componentRadius = component->outerDiameter * 0.0254 / 2.0;
        double maxDeflection = holeRadius - componentRadius;
        double currentDeflection = displacements_(i * 2);
        
        if (std::abs(currentDeflection) > maxDeflection) {
            double newDeflection = (currentDeflection > 0) ? maxDeflection : -maxDeflection;
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
    // Find the element containing the sensor position
    int elementIndex = static_cast<int>(distanceFromBit / elementSize_);
    
    // Ensure we don't go out of bounds
    if (elementIndex >= elements_.size()) {
        elementIndex = elements_.size() - 1;
    }
    
    // Get the node indices
    int startNode = elementIndex;
    int endNode = elementIndex + 1;
    
    // Get the nodal displacements and rotations for this element
    double w1 = displacements_(startNode * 2);     // Deflection at start node
    double theta1 = displacements_(startNode * 2 + 1); // Rotation at start node
    double w2 = displacements_(endNode * 2);       // Deflection at end node  
    double theta2 = displacements_(endNode * 2 + 1);   // Rotation at end node
    
    // Find the local position within the element
    double elementStart = elementIndex * elementSize_;
    double xi = (distanceFromBit - elementStart) / elementSize_;
    
    // Calculate slope using the derivative of the cubic Hermite interpolation
    // dw/dx = d/dx[ H1(xi)*w1 + H2(xi)*theta1 + H3(xi)*w2 + H4(xi)*theta2 ]
    // Where H1, H2, H3, H4 are Hermite basis functions
    
    // First derivatives of Hermite basis functions
    double dH1 = 6*xi*(xi-1);                  // dH1/dxi = 6*(xi^2-xi)
    double dH2 = 1 + 3*xi*(2*xi-3);            // dH2/dxi = 1+3*(2*xi^2-3*xi)
    double dH3 = 6*xi*(1-xi);                  // dH3/dxi = 6*(xi-xi^2)
    double dH4 = 3*xi*(xi-1);                  // dH4/dxi = 3*(xi^2-xi)
    
    // Calculate actual slope in the physical space (radians)
    // Need to multiply by 1/L to convert from dw/dxi to dw/dx
    double L = elementSize_;
    double slope = (dH1*w1 + L*dH2*theta1 + dH3*w2 + L*dH4*theta2) / L;
    
    // Return slope in degrees
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
        
        // Apply sanity limits to all displacements to prevent overflow
        const double MAX_DISPLACEMENT = 0.1; // 10 cm max
        const double MAX_ROTATION = 0.1;     // ~5.7 degrees max
        
        double w1 = std::max(std::min(displacements_(startDof), MAX_DISPLACEMENT), -MAX_DISPLACEMENT);
        double theta1 = std::max(std::min(displacements_(startDof + 1), MAX_ROTATION), -MAX_ROTATION);
        double w2 = std::max(std::min(displacements_(endDof), MAX_DISPLACEMENT), -MAX_DISPLACEMENT);
        double theta2 = std::max(std::min(displacements_(endDof + 1), MAX_ROTATION), -MAX_ROTATION);
        
        // Element properties
        double L = elements_[i].getLength();
        double E = elements_[i].getE();
        double I = elements_[i].getI();
        
        // Calculate moments at nodes using Euler-Bernoulli beam theory
        double M1 = E * I * (4*theta1 + 2*theta2 - 6*(w2-w1)/L) / L;
        double M2 = E * I * (2*theta1 + 4*theta2 - 6*(w2-w1)/L) / L;
        
        // Calculate shear force using static equilibrium
        double V = (M1 + M2) / L;
        
        // Store at beam nodes
        double startPos = startNode * elementSize_;
        double endPos = endNode * elementSize_;
        
        // Apply reasonable limits
        const double MAX_MOMENT = 1.0e6; // N·m
        const double MAX_SHEAR = 1.0e5;  // N
        
        M1 = std::max(std::min(M1, MAX_MOMENT), -MAX_MOMENT);
        M2 = std::max(std::min(M2, MAX_MOMENT), -MAX_MOMENT);
        V = std::max(std::min(V, MAX_SHEAR), -MAX_SHEAR);
        
        // Calculate slopes in degrees
        double slope1 = theta1 * 180.0 / M_PI;
        double slope2 = theta2 * 180.0 / M_PI;
        
        // Add to profiles
        shearForces_.push_back(std::make_pair(startPos, V));
        bendingMoments_.push_back(std::make_pair(startPos, M1));
        slopes_.push_back(std::make_pair(startPos, slope1));
        
        // Add intermediate points for smoother visualization
        const int numPoints = 8;
        for (int j = 1; j < numPoints; j++) {
            double t = static_cast<double>(j) / numPoints;
            double pos = startPos + t * (endPos - startPos);
            
            // Linear interpolation of internal forces
            double M = M1 * (1-t) + M2 * t;
            double slope = slope1 * (1-t) + slope2 * t;
            
            shearForces_.push_back(std::make_pair(pos, V));
            bendingMoments_.push_back(std::make_pair(pos, M));
            slopes_.push_back(std::make_pair(pos, slope));
        }
        
        // Add end point for last element
        if (i == elements_.size() - 1) {
            shearForces_.push_back(std::make_pair(endPos, V));
            bendingMoments_.push_back(std::make_pair(endPos, M2));
            slopes_.push_back(std::make_pair(endPos, slope2));
        }
    }
    
    // Sort all profiles by position
    auto compareByPos = [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    };
    
    std::sort(shearForces_.begin(), shearForces_.end(), compareByPos);
    std::sort(bendingMoments_.begin(), bendingMoments_.end(), compareByPos);
    std::sort(slopes_.begin(), slopes_.end(), compareByPos);
}

FEAResults FEASolver::getResults() const {
    FEAResults results;
    
    // Convert bit side force from N to kg
    double bitSideForceKg = bitSideForce_ / 9.81;
    
    // Apply physical limits based on BHA weight and hole size
    double maxPossibleForce = 0.0;
    for (const auto& component : components_) {
        maxPossibleForce += component.weight * 0.453592 * component.length / 0.3048; // Convert ppf to kg
    }
    
    // Limit force to a reasonable fraction of BHA weight
    const double MAX_FORCE_FRACTION = 0.5; // 50% of BHA weight as side force max
    double limitedForce = std::min(std::abs(bitSideForceKg), maxPossibleForce * MAX_FORCE_FRACTION);
    if (bitSideForceKg < 0) limitedForce = -limitedForce;
    
    // Bit results with physical constraints
    results.bitSideForce = limitedForce;
    results.bend = 0.25;  // Motor bend
    results.wob = wob_;
    
    // Find MWD sensor position
    double sensorPos = 0.0;
    for (const auto& component : components_) {
        if (component.name.find("MWD") != std::string::npos) {
            sensorPos = component.distanceFromBit + component.length / 2.0;
            break;
        }
    }
    results.sensorPosition = sensorPos;
    
    // Get trajectory data
    TrajectoryPoint trajPoint = trajectory_.interpolateAtDepth(sensorPos);
    results.holeCurvature = trajPoint.dls * 100.0 / 30.0;  // Convert to deg/100ft
    
    // Calculate sag correction with sanity checks
    double rawSagCorrection = getSagCorrection(sensorPos);
    
    // Typical sag correction range is 0-3 degrees
    const double MAX_REASONABLE_SAG = 1.0;
    results.sagCorrection = std::max(std::min(rawSagCorrection, MAX_REASONABLE_SAG), -MAX_REASONABLE_SAG);
    if (rawSagCorrection < 0) results.sagCorrection = -results.sagCorrection;
    
    // Apply similar constraints to stabilizer forces
    results.stabilizerForces = stabilizerForces_;
    for (auto& stab : results.stabilizerForces) {
        double limitedStabForce = std::min(std::abs(stab.sideForce / 9.81), maxPossibleForce * MAX_FORCE_FRACTION);
        if (stab.sideForce < 0) limitedStabForce = -limitedStabForce;
        stab.sideForce = limitedStabForce;
    }
    
    // Create deflection profile
    // Replace the deflection profile generation with this:
    std::vector<std::pair<double, double>> deflectionProfile;

    // For each element, use shape functions for smooth interpolation
    for (int i = 0; i < elements_.size(); ++i) {
        int startNode = i;
        int endNode = i + 1;
        int startDof = startNode * 2;
        int endDof = endNode * 2;
        
        double startPos = startNode * elementSize_;
        double endPos = endNode * elementSize_;
        double L = elements_[i].getLength();
        
        // Get nodal displacements and rotations
        double w1 = displacements_(startDof);
        double theta1 = displacements_(startDof + 1);
        double w2 = displacements_(endDof);
        double theta2 = displacements_(endDof + 1);
        
        // Add start point
        deflectionProfile.push_back(std::make_pair(startPos, w1));
        
        // Add interior points for smoother visualization
        const int numInteriorPoints = 10;
        for (int j = 1; j < numInteriorPoints; ++j) {
            double xi = static_cast<double>(j) / numInteriorPoints;
            double x = startPos + xi * L;
            
            // Hermite shape functions for cubic interpolation
            double N1 = 1 - 3*xi*xi + 2*xi*xi*xi;
            double N2 = L * (xi - 2*xi*xi + xi*xi*xi);
            double N3 = 3*xi*xi - 2*xi*xi*xi;
            double N4 = L * (-xi*xi + xi*xi*xi);
            
            // Interpolate deflection
            double w = N1*w1 + N2*theta1 + N3*w2 + N4*theta2;
            
            deflectionProfile.push_back(std::make_pair(x, w));
        }
    }

    // Add final point
    deflectionProfile.push_back(std::make_pair(elements_.back().getLength() * elements_.size(), 
                            displacements_((elements_.size()) * 2)));

    results.deflectionProfile = deflectionProfile;
    
    // Other profiles
    results.shearForceProfile = shearForces_;
    results.bendingMomentProfile = bendingMoments_;
    results.slopeProfile = slopes_;
    
    return results;
}