// src/calculators/fea_core/bha_component.h
#pragma once
#include <string>

struct BHAComponent {
    std::string name;        // Component name
    double length;           // Length in meters
    double outerDiameter;    // Outer diameter in inches
    double innerDiameter;    // Inner diameter in inches
    double weight;           // Weight in ppf (pounds per foot)
    std::string material;    // Material grade (e.g., "SAE 4145")
    double distanceFromBit;  // Distance from bit in meters
    bool isStabilizer;       // Flag for stabilizer components
};

struct MaterialProperties {
    std::string grade;       // Material grade name
    double youngsModulus;    // Young's modulus in Pa
    double density;          // Density in kg/m³
    double poissonRatio;     // Poisson's ratio (dimensionless)
    double yieldStrength;    // Yield strength in Pa
};
