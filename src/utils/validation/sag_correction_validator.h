// src/utils/validation/sag_correction_validator.h
#pragma once
#include "../../calculators/fea_core/fea_solver.h"
#include "../../calculators/fea_core/trajectory.h"
#include <string>
#include <vector>
#include <utility>

struct ValidationResult {
    double averageError;
    double maxError;
    double rmseError;
    std::vector<std::pair<double, double>> errorProfile; // (position, error)
};

class SagCorrectionValidator {
public:
    SagCorrectionValidator(const FEASolver& solver);
    
    // Validate against laboratory data
    ValidationResult validateAgainstLabData(const std::string& labDataFile);
    
    // Validate against dual inclination data
    ValidationResult validateAgainstDualInclination(
        const std::vector<TrajectoryPoint>& mwdSurveys,
        const std::vector<TrajectoryPoint>& rssSurveys);
    
    // Validate against theoretical closed-form solution
    ValidationResult validateAgainstTheoreticalSolution();
    
    // Write validation results to file
    void writeResultsToFile(const ValidationResult& result, const std::string& filename);
    
private:
    const FEASolver& solver_;
    
    // Helper method to calculate error metrics
    ValidationResult calculateErrorMetrics(
        const std::vector<std::pair<double, double>>& predicted,
        const std::vector<std::pair<double, double>>& actual);
};