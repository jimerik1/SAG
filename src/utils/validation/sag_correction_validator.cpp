// src/utils/validation/sag_correction_validator.cpp
#include "sag_correction_validator.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

SagCorrectionValidator::SagCorrectionValidator(const FEASolver& solver)
    : solver_(solver) {
}

ValidationResult SagCorrectionValidator::validateAgainstLabData(const std::string& labDataFile) {
    // Read laboratory data file
    std::ifstream file(labDataFile);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open lab data file: " + labDataFile);
    }
    
    // Parse lab data
    std::vector<std::pair<double, double>> labData; // (position, inclination)
    std::string line;
    
    // Skip header
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        // Parse position
        std::getline(ss, token, ',');
        std::replace(token.begin(), token.end(), ',', '.');
        double position = std::stod(token);
        
        // Parse inclination
        std::getline(ss, token, ',');
        std::replace(token.begin(), token.end(), ',', '.');
        double inclination = std::stod(token);
        
        labData.push_back(std::make_pair(position, inclination));
    }
    
    // Get model predictions
    std::vector<std::pair<double, double>> predictions;
    for (const auto& point : labData) {
        double position = point.first;
        double predictedInc = solver_.getSagCorrection(position);
        predictions.push_back(std::make_pair(position, predictedInc));
    }
    
    // Calculate error metrics
    return calculateErrorMetrics(predictions, labData);
}

ValidationResult SagCorrectionValidator::validateAgainstDualInclination(
    const std::vector<TrajectoryPoint>& mwdSurveys,
    const std::vector<TrajectoryPoint>& rssSurveys) {
    
    // Extract MWD inclinations
    std::vector<std::pair<double, double>> mwdData;
    for (const auto& survey : mwdSurveys) {
        mwdData.push_back(std::make_pair(survey.md, survey.inclination));
    }
    
    // Extract RSS inclinations
    std::vector<std::pair<double, double>> rssData;
    for (const auto& survey : rssSurveys) {
        rssData.push_back(std::make_pair(survey.md, survey.inclination));
    }
    
    // Find matching depth points
    std::vector<std::pair<double, double>> matchedMwd;
    std::vector<std::pair<double, double>> matchedRss;
    
    for (const auto& mwd : mwdData) {
        // Find closest RSS survey
        auto closest = std::min_element(rssData.begin(), rssData.end(),
            [&mwd](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                return std::abs(a.first - mwd.first) < std::abs(b.first - mwd.first);
            });
        
        if (closest != rssData.end() && std::abs(closest->first - mwd.first) < 0.5) {
            matchedMwd.push_back(mwd);
            matchedRss.push_back(*closest);
        }
    }
    
    // Get model predictions for MWD points
    std::vector<std::pair<double, double>> predictions;
    for (const auto& point : matchedMwd) {
        double position = point.first;
        double predictedInc = solver_.getSagCorrection(position);
        predictions.push_back(std::make_pair(position, predictedInc + point.second));
    }
    
    // Calculate error metrics against RSS readings
    return calculateErrorMetrics(predictions, matchedRss);
}

ValidationResult SagCorrectionValidator::validateAgainstTheoreticalSolution() {
    // Implement a theoretical closed-form solution for a simple case
    // For example, a beam with one end fixed and the other end free
    
    // Theoretical solution for a cantilever beam with uniform load
    // deflection = (w * x^2 * (6*L^2 - 4*L*x + x^2)) / (24 * E * I)
    // where w is the weight per unit length, L is the length, E is Young's modulus,
    // I is moment of inertia, and x is the position along the beam
    
    // Generate theoretical deflection profile
    std::vector<std::pair<double, double>> theoretical;
    std::vector<std::pair<double, double>> predicted;
    
    double L = 10.0; // 10m length
    double E = 2.1e11; // Young's modulus for steel
    double OD = 0.12; // 120mm outer diameter
    double ID = 0.04; // 40mm inner diameter
    double I = M_PI * (pow(OD, 4) - pow(ID, 4)) / 64.0; // Moment of inertia
    double rho = 7850.0; // Density of steel
    double A = M_PI * (pow(OD, 2) - pow(ID, 2)) / 4.0; // Cross-sectional area
    double w = rho * A * 9.81; // Weight per unit length
    
    // Calculate theoretical deflections
    for (double x = 0; x <= L; x += 0.1) {
        double deflection = (w * pow(x, 2) * (6 * pow(L, 2) - 4 * L * x + pow(x, 2))) / (24 * E * I);
        double angle = atan(deflection / x) * 180.0 / M_PI; // Convert to degrees
        theoretical.push_back(std::make_pair(x, angle));
    }
    
    // Get model predictions
    for (const auto& point : theoretical) {
        double position = point.first;
        double predictedInc = solver_.getSagCorrection(position);
        predicted.push_back(std::make_pair(position, predictedInc));
    }
    
    // Calculate error metrics
    return calculateErrorMetrics(predicted, theoretical);
}

ValidationResult SagCorrectionValidator::calculateErrorMetrics(
    const std::vector<std::pair<double, double>>& predicted,
    const std::vector<std::pair<double, double>>& actual) {
    
    ValidationResult result;
    result.errorProfile.clear();
    
    double sumError = 0.0;
    double sumSquaredError = 0.0;
    result.maxError = 0.0;
    
    for (size_t i = 0; i < predicted.size() && i < actual.size(); i++) {
        double error = std::abs(predicted[i].second - actual[i].second);
        sumError += error;
        sumSquaredError += error * error;
        
        if (error > result.maxError) {
            result.maxError = error;
        }
        
        result.errorProfile.push_back(std::make_pair(predicted[i].first, error));
    }
    
    result.averageError = sumError / predicted.size();
    result.rmseError = std::sqrt(sumSquaredError / predicted.size());
    
    return result;
}

void SagCorrectionValidator::writeResultsToFile(const ValidationResult& result, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open output file: " + filename);
    }
    
    file << "Validation Results" << std::endl;
    file << "================" << std::endl;
    file << "Average Error: " << result.averageError << " degrees" << std::endl;
    file << "Maximum Error: " << result.maxError << " degrees" << std::endl;
    file << "RMSE: " << result.rmseError << " degrees" << std::endl;
    file << std::endl;
    
    file << "Error Profile" << std::endl;
    file << "=============" << std::endl;
    file << "Position(m),Error(deg)" << std::endl;
    
    for (const auto& point : result.errorProfile) {
        file << point.first << "," << point.second << std::endl;
    }
    
    file.close();
}