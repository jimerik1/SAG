// src/calculators/beam_elements/beam_element.cpp
#include "beam_element.h"
#include <cmath>

BeamElement::BeamElement(double length, double E, double I, double A, double rho)
    : length_(length), E_(E), I_(I), A_(A), rho_(rho) {}

Eigen::MatrixXd BeamElement::calculateStiffnessMatrix() const {
    // For a 2D Euler-Bernoulli beam element with 2 nodes (4 DOFs per element)
    // DOFs: [y1, θ1, y2, θ2] where y is deflection and θ is rotation
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(4, 4);
    
    double L = length_;
    double EI = E_ * I_;
    
    // Assemble stiffness matrix for beam bending
    // This is the standard stiffness matrix for Euler-Bernoulli beam
    K(0, 0) = 12 * EI / pow(L, 3);
    K(0, 1) = 6 * EI / pow(L, 2);
    K(0, 2) = -12 * EI / pow(L, 3);
    K(0, 3) = 6 * EI / pow(L, 2);
    
    K(1, 0) = 6 * EI / pow(L, 2);
    K(1, 1) = 4 * EI / L;
    K(1, 2) = -6 * EI / pow(L, 2);
    K(1, 3) = 2 * EI / L;
    
    K(2, 0) = -12 * EI / pow(L, 3);
    K(2, 1) = -6 * EI / pow(L, 2);
    K(2, 2) = 12 * EI / pow(L, 3);
    K(2, 3) = -6 * EI / pow(L, 2);
    
    K(3, 0) = 6 * EI / pow(L, 2);
    K(3, 1) = 2 * EI / L;
    K(3, 2) = -6 * EI / pow(L, 2);
    K(3, 3) = 4 * EI / L;
    
    return K;
}

Eigen::VectorXd BeamElement::calculateGravityLoad(double inclination, double azimuth) const {
    // Convert inclination and azimuth to radians
    double inc_rad = inclination * M_PI / 180.0;
    
    // Calculate component of gravity perpendicular to wellbore axis
    // This is the force that causes sag
    double gravity = 9.81;  // m/s²
    double weight_per_length = rho_ * A_ * gravity;  // N/m
    
    // Weight component perpendicular to wellbore axis
    double perpendicular_force = weight_per_length * sin(inc_rad);
    
    // For a 2D beam element, distributed load creates a load vector
    Eigen::VectorXd F = Eigen::VectorXd::Zero(4);
    
    // For uniformly distributed load w, the equivalent nodal loads are:
    // F = [wL/2, wL²/12, wL/2, -wL²/12]ᵀ
    F(0) = perpendicular_force * length_ / 2.0;
    F(1) = perpendicular_force * pow(length_, 2) / 12.0;
    F(2) = perpendicular_force * length_ / 2.0;
    F(3) = -perpendicular_force * pow(length_, 2) / 12.0;
    
    return F;
}