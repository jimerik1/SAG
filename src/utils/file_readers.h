// src/calculators/utils/file_readers.h
#pragma once
#include <string>
#include <vector>
#include "../calculators/fea_core/bha_component.h"
#include "../calculators/fea_core/trajectory.h"

class FileReaders {
public:
    // Load trajectory data from file
    static Trajectory loadTrajectoryFromFile(const std::string& filename);
    
    // Load BHA components from file
    static std::vector<BHAComponent> loadBHAComponentsFromFile(const std::string& filename);
    
    // Load material properties from file
    static std::vector<MaterialProperties> loadMaterialPropertiesFromFile(const std::string& filename);
};