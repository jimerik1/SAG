// src/calculators/fea_core/fea_results.h
#pragma once
#include <vector>
#include <string>
#include <utility>

struct StabilizerForce {
    std::string name;
    double position;  // distance from bit
    double sideForce; // in kg or N
};

class FEAResults {
public:
    // Bit results
    double bitSideForce;  // kg or N
    double bend;          // degrees
    double wob;           // Weight on bit (kg or N)
    
    // Sensor results
    double sensorPosition;     // m
    double holeCurvature;      // deg/100ft
    double sagCorrection;      // degrees
    
    // Stabilizer forces
    std::vector<StabilizerForce> stabilizerForces;
    
    // Profile data
    std::vector<std::pair<double, double>> deflectionProfile;    // (position, deflection)
    std::vector<std::pair<double, double>> shearForceProfile;    // (position, shear)
    std::vector<std::pair<double, double>> bendingMomentProfile; // (position, moment)
    std::vector<std::pair<double, double>> slopeProfile;         // (position, slope in degrees)
    
    // Output summary to console
    void printSummary() const;
    
    // Write detailed results to files
    void writeProfilesToFiles(const std::string& directory) const;
};