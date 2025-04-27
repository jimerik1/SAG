// src/calculators/fea_core/fea_solver.h
#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "bha_component.h"
#include "trajectory.h"
#include "fea_results.h"
#include "../beam_elements/beam_element.h"
#include "../../utils/file_readers.h" // Add this include to access WellboreProperties

class FEASolver {
public:
    FEASolver(const std::vector<BHAComponent>& components,
              const std::vector<MaterialProperties>& materials,
              const Trajectory& trajectory,
              const WellboreProperties& wellboreProps,
              double elementSize = 0.1);
    
    // Set up the FEA system
    void setupSystem();
    
    // Solve the FEA system
    void solve();
    
    // Get sag correction at a specific distance from bit
    double getSagCorrection(double distanceFromBit) const;
    
    // Get deflection at a specific distance from bit
    double getDeflection(double distanceFromBit) const;
    
    // Get complete deflection profile
    std::vector<std::pair<double, double>> getDeflectionProfile() const;
    
    // Get complete results
    FEAResults getResults() const;
    
private:
    // Input data
    std::vector<BHAComponent> components_;
    std::vector<MaterialProperties> materials_;
    Trajectory trajectory_;
    double elementSize_;
    double mudWeight_;
    double wob_;
    WellboreProperties wellboreProps_;

    
    // FEA system
    std::vector<BeamElement> elements_;
    Eigen::SparseMatrix<double> globalStiffnessMatrix_;
    Eigen::VectorXd globalLoadVector_;
    Eigen::VectorXd displacements_;
    
    void applyWellboreContactConstraints();
    const BHAComponent* getComponentAtPosition(double position) const;

    
    // Mesh information
    int numNodes_;
    int numDofs_;
    std::vector<double> nodePositions_;
    
    // Results
    std::vector<std::pair<double, double>> shearForces_;
    std::vector<std::pair<double, double>> bendingMoments_;
    std::vector<std::pair<double, double>> slopes_;
    double bitSideForce_;
    std::vector<StabilizerForce> stabilizerForces_;
    
    // Create mesh
    void createMesh();
    
    // Assemble global stiffness matrix
    void assembleStiffnessMatrix();
    
    // Assemble global load vector
    void assembleLoadVector();
    
    // Apply boundary conditions
    void applyBoundaryConditions();
    
    // Calculate internal forces from displacements
    void calculateInternalForces();
    
    // Calculate side forces at stabilizers and bit
    void calculateSideForces();
    
    // Calculate material properties for each element
    MaterialProperties getMaterialProperties(const std::string& grade) const;
    
    // Calculate area moment of inertia
    double calculateMomentOfInertia(double od, double id) const;
    
    // Calculate cross-sectional area
    double calculateArea(double od, double id) const;
};