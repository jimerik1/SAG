// src/utils/file_readers.h
#pragma once
#include <string>
#include <vector>
#include "../calculators/fea_core/bha_component.h"
#include "../calculators/fea_core/trajectory.h"

// Define the WellboreProperties struct
struct WellboreProperties {
    double holeSize;             // in inches
    double mudWeight;            // in pounds per gallon (ppg)
    double contactStiffness;     // in N/m
    int maxIterations;           // number of iterations for contact solution
    double convergenceTolerance; // in meters
};

class FileReaders {
public:
    // Load trajectory data from file
    static Trajectory loadTrajectoryFromFile(const std::string& filename);
    
    // Load BHA components from file
    static std::vector<BHAComponent> loadBHAComponentsFromFile(const std::string& filename);
    
    // Load material properties from file
    static std::vector<MaterialProperties> loadMaterialPropertiesFromFile(const std::string& filename);
    
    // Load wellbore properties from file
    static WellboreProperties loadWellborePropertiesFromFile(const std::string& filename);
};