// src/calculators/beam_elements/beam_element.h
#pragma once
#include <Eigen/Dense>
#include <vector>

class BeamElement {
public:
    BeamElement(double length, double E, double I, double A, double rho);
    
    // Calculate element stiffness matrix
    Eigen::MatrixXd calculateStiffnessMatrix() const;
    
    // Calculate element mass matrix
    Eigen::MatrixXd calculateMassMatrix() const;
    
    // Calculate element load vector due to gravity
    Eigen::VectorXd calculateGravityLoad(double inclination, double azimuth) const;
    
    // Calculate element load vector due to buoyancy
    Eigen::VectorXd calculateBuoyancyLoad(double mudWeight, double inclination, double azimuth) const;
    
private:
    double length_;     // Element length
    double E_;          // Young's modulus
    double I_;          // Area moment of inertia
    double A_;          // Cross-sectional area
    double rho_;        // Density
};