// src/calculators/fea_core/tvd_assurance_workflow.h
#pragma once
#include <vector>
#include "trajectory.h"
#include "bha_component.h"
#include "fea_solver.h"
#include "slide_sheet.h"

class TvdAssuranceWorkflow {
public:
    TvdAssuranceWorkflow(
        const std::vector<TrajectoryPoint>& staticSurveys,
        const std::vector<TrajectoryPoint>& continuousInclination,
        const std::vector<BHAComponent>& bhaComponents,
        const std::vector<MaterialProperties>& materials,
        const SlideSheet& slideSheet,
        double mudWeight = 12.0,
        double wob = 20000.0);
    
    std::vector<TrajectoryPoint> process();
    
    // Get the final corrected trajectory
    std::vector<TrajectoryPoint> getCorrectedTrajectory() const;
    
    // Get delta inclination for each static survey
    std::vector<std::pair<double, double>> getDeltaInclinations() const;
    
private:
    // Input data
    std::vector<TrajectoryPoint> staticSurveys_;
    std::vector<TrajectoryPoint> contIncl_;
    std::vector<BHAComponent> bhaComponents_;
    std::vector<MaterialProperties> materials_;
    SlideSheet slideSheet_;
    double mudWeight_;
    double wob_;
    
    // Output data
    std::vector<TrajectoryPoint> correctedTrajectory_;
    std::vector<std::pair<double, double>> deltaInclinations_; // (md, delta_inc)
    
    // Convergence parameters
    double convergenceTolerance_;
    int maxIterations_;
    
    // Methods for sub-algorithms
    std::vector<TrajectoryPoint> novelMsa(const std::vector<TrajectoryPoint>& rawSurveys);
    std::vector<TrajectoryPoint> hdTrajectory(
        const std::vector<TrajectoryPoint>& correctedSurveys,
        const std::vector<TrajectoryPoint>& contIncl);
    std::vector<TrajectoryPoint> sagCorrection(
        const std::vector<TrajectoryPoint>& hdTraj,
        const std::vector<TrajectoryPoint>& approxTraj);
    std::vector<TrajectoryPoint> approxTrajectory(
        const std::vector<TrajectoryPoint>& noSagTraj);
    
    // Helper methods
    double calculateTrajDifference(
        const std::vector<TrajectoryPoint>& traj1,
        const std::vector<TrajectoryPoint>& traj2);
};