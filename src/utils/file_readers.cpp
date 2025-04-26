// src/calculators/utils/file_readers.cpp
#include "file_readers.h"
#include <fstream>
#include <sstream>
#include <iostream>

Trajectory FileReaders::loadTrajectoryFromFile(const std::string& filename) {
    Trajectory trajectory;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open trajectory file: " + filename);
    }
    
    // Skip header lines (first 2 lines)
    std::string line;
    std::getline(file, line); // Column names
    std::getline(file, line); // Units
    
    // Read data
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        TrajectoryPoint point;
        
        // Parse MD
        std::getline(ss, token, '\t');
        // Replace comma with dot for decimal point
        std::replace(token.begin(), token.end(), ',', '.');
        point.md = std::stod(token);
        
        // Parse inclination
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        point.inclination = std::stod(token);
        
        // Parse azimuth
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        point.azimuth = std::stod(token);
        
        // Parse TVD
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        point.tvd = std::stod(token);
        
        // Skip EW, NS, V.SEC
        std::getline(ss, token, '\t');
        std::getline(ss, token, '\t');
        std::getline(ss, token, '\t');
        
        // Parse DLS
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        point.dls = std::stod(token);
        
        trajectory.points.push_back(point);
    }
    
    return trajectory;
}

std::vector<BHAComponent> FileReaders::loadBHAComponentsFromFile(const std::string& filename) {
    std::vector<BHAComponent> components;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open BHA components file: " + filename);
    }
    
    // Skip header lines (first 2 lines)
    std::string line;
    std::getline(file, line); // Column names
    std::getline(file, line); // Units
    
    // Read data
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        BHAComponent component;
        
        // Parse component name
        std::getline(ss, token, '\t');
        component.name = token;
        
        // Parse length
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        component.length = std::stod(token);
        
        // Parse outer diameter
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        component.outerDiameter = std::stod(token);
        
        // Parse inner diameter
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        component.innerDiameter = std::stod(token);
        
        // Parse weight
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        component.weight = std::stod(token);
        
        // Parse material grade
        std::getline(ss, token, '\t');
        component.material = token;
        
        // Parse is_stabilizer flag
        std::getline(ss, token, '\t');
        component.isStabilizer = (token == "True" || token == "true" || token == "TRUE" || token == "1");
        
        // Parse distance from bit
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        component.distanceFromBit = std::stod(token);
        
        components.push_back(component);
    }
    
    return components;
}

std::vector<MaterialProperties> FileReaders::loadMaterialPropertiesFromFile(const std::string& filename) {
    std::vector<MaterialProperties> materials;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open material properties file: " + filename);
    }
    
    // Skip header lines (first 2 lines)
    std::string line;
    std::getline(file, line); // Column names
    std::getline(file, line); // Units
    
    // Read data
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        MaterialProperties material;
        
        // Parse grade
        std::getline(ss, token, '\t');
        material.grade = token;
        
        // Parse Young's modulus
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        material.youngsModulus = std::stod(token);
        
        // Parse density
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        material.density = std::stod(token);
        
        // Parse Poisson's ratio
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        material.poissonRatio = std::stod(token);
        
        // Parse yield strength
        std::getline(ss, token, '\t');
        std::replace(token.begin(), token.end(), ',', '.');
        material.yieldStrength = std::stod(token);
        
        materials.push_back(material);
    }
    
    return materials;
}