// src/calculators/fea_core/tvd_assurance_workflow.cpp
#include "tvd_assurance_workflow.h"
#include <cmath>
#include <iostream>
#include <algorithm>

TvdAssuranceWorkflow::TvdAssuranceWorkflow(
    const std::vector<TrajectoryPoint>& staticSurveys,
    const std::vector<TrajectoryPoint>& continuousInclination,
    const std::vector<BHAComponent>& bhaComponents,
    const std::vector<MaterialProperties>& materials,
    const SlideSheet& slideSheet,
    double mudWeight,
    double wob)
    : staticSurveys_(staticSurveys),
      contIncl_(continuousInclination),
      bhaComponents_(bhaComponents),
      materials_(materials),
      slideSheet_(slideSheet),
      mudWeight_(mudWeight),
      wob_(wob),
      convergenceTolerance_(0.01),
      maxIterations_(10) {
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::process() {
    std::cout << "Starting TVD assurance workflow..." << std::endl;
    
    // Step 1: Correct static surveys for accelerometer errors
    std::cout << "Step 1: Correcting accelerometer errors using Novel MSA..." << std::endl;
    std::vector<TrajectoryPoint> correctedSurveys = novelMsa(staticSurveys_);
    
    // Step 2: Calculate initial HD trajectory
    std::cout << "Step 2: Creating high-definition trajectory..." << std::endl;
    std::vector<TrajectoryPoint> hdTraj = hdTrajectory(correctedSurveys, contIncl_);
    
    // Step 3: Initial sag correction
    std::cout << "Step 3: Initial sag correction..." << std::endl;
    std::vector<TrajectoryPoint> noSagTraj = sagCorrection(hdTraj, hdTraj);
    
    // Step 4: Initial approximated trajectory
    std::cout << "Step 4: Approximating trajectory with slide sheet data..." << std::endl;
    std::vector<TrajectoryPoint> approxTraj = approxTrajectory(noSagTraj);
    
    // Steps 5-7: Iterative refinement
    std::cout << "Step 5: Starting iterative refinement..." << std::endl;
    
    for (int i = 0; i < maxIterations_; i++) {
        // Save previous iteration
        std::vector<TrajectoryPoint> prevNoSagTraj = noSagTraj;
        
        // Update with refined approximation
        std::cout << "  Iteration " << i+1 << "..." << std::endl;
        noSagTraj = sagCorrection(hdTraj, approxTraj);
        approxTraj = approxTrajectory(noSagTraj);
        
        // Check convergence
        double change = calculateTrajDifference(noSagTraj, prevNoSagTraj);
        std::cout << "  Change: " << change << std::endl;
        
        if (change < convergenceTolerance_) {
            std::cout << "Converged after " << i+1 << " iterations." << std::endl;
            break;
        }
    }
    
    // Store final corrected trajectory
    correctedTrajectory_ = approxTraj;
    
    // Calculate delta inclinations for each static survey
    deltaInclinations_.clear();
    for (const auto& survey : staticSurveys_) {
        double md = survey.md;
        double originalInc = survey.inclination;
        
        // Find corresponding corrected inclination
        auto it = std::find_if(approxTraj.begin(), approxTraj.end(),
            [md](const TrajectoryPoint& pt) { return std::abs(pt.md - md) < 0.01; });
        
        if (it != approxTraj.end()) {
            double correctedInc = it->inclination;
            double deltaInc = correctedInc - originalInc;
            deltaInclinations_.push_back(std::make_pair(md, deltaInc));
        }
    }
    
    std::cout << "TVD assurance workflow completed." << std::endl;
    
    return correctedTrajectory_;
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::getCorrectedTrajectory() const {
    return correctedTrajectory_;
}

std::vector<std::pair<double, double>> TvdAssuranceWorkflow::getDeltaInclinations() const {
    return deltaInclinations_;
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::novelMsa(const std::vector<TrajectoryPoint>& rawSurveys) {
    // This would be a full implementation of the Novel MSA algorithm
    // For now, we'll use a placeholder that returns the input with small corrections
    std::vector<TrajectoryPoint> corrected = rawSurveys;
    
    // Apply a simple correction to illustrate the concept
    // In a real implementation, this would be replaced with the full Novel MSA algorithm
    for (auto& point : corrected) {
        if (point.inclination > 85.0) {
            // Apply more correction in the horizontal section
            point.inclination += 0.04 * sin(point.md * 0.01); // Simulated correction
        }
    }
    
    return corrected;
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::hdTrajectory(
    const std::vector<TrajectoryPoint>& correctedSurveys,
    const std::vector<TrajectoryPoint>& contIncl) {
    
    // This would be a full implementation of the HD trajectory algorithm
    // For now, create a high-density trajectory by merging and interpolating
    std::vector<TrajectoryPoint> hdTraj = correctedSurveys;
    
    // Add continuous inclination data to create higher density
    for (const auto& point : contIncl) {
        // Check if this md already exists in hdTraj
        auto it = std::find_if(hdTraj.begin(), hdTraj.end(),
            [&point](const TrajectoryPoint& p) { return std::abs(p.md - point.md) < 0.01; });
        
        if (it == hdTraj.end()) {
            // Find neighboring static surveys
            auto next = std::find_if(correctedSurveys.begin(), correctedSurveys.end(),
                [&point](const TrajectoryPoint& p) { return p.md > point.md; });
            
            if (next != correctedSurveys.end() && next != correctedSurveys.begin()) {
                auto prev = next - 1;
                
                // Blend continuous inclination with interpolated static surveys
                double weight = (point.md - prev->md) / (next->md - prev->md);
                double staticInc = prev->inclination + weight * (next->inclination - prev->inclination);
                
                // Create a new point, weighting static higher than continuous
                TrajectoryPoint newPoint = point;
                newPoint.inclination = 0.8 * staticInc + 0.2 * point.inclination;
                
                hdTraj.push_back(newPoint);
            }
        }
    }
    
    // Sort by measured depth
    std::sort(hdTraj.begin(), hdTraj.end(),
        [](const TrajectoryPoint& a, const TrajectoryPoint& b) { return a.md < b.md; });
    
    return hdTraj;
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::sagCorrection(
    const std::vector<TrajectoryPoint>& hdTraj,
    const std::vector<TrajectoryPoint>& approxTraj) {
    
    // Create trajectory from the HD trajectory points
    Trajectory trajectory;
    trajectory.points = hdTraj;
    
    // Create FEA solver with the current approximated trajectory for local curvature
    double elementSize = 0.1;  // 10cm elements
    
    FEASolver solver(bhaComponents_, materials_, trajectory, elementSize, mudWeight_, wob_);
    
    // Set up and solve FEA system
    solver.setupSystem();
    solver.solve();
    
    // Get BHA sag at each trajectory point
    std::vector<TrajectoryPoint> correctedTraj = hdTraj;
    
    for (auto& point : correctedTraj) {
        double md = point.md;
        double sagCorrection = solver.getSagCorrection(md);
        
        // Apply the correction (sag is added to the raw inclination)
        point.inclination += sagCorrection;
    }
    
    return correctedTraj;
}

std::vector<TrajectoryPoint> TvdAssuranceWorkflow::approxTrajectory(
    const std::vector<TrajectoryPoint>& noSagTraj) {
    
    // This would be the full implementation of the piecewise linear approximation
    // For now, we'll use a simplified version that relies on slide sheet data
    
    // Create a new trajectory with linear segments based on slide sheet data
    std::vector<TrajectoryPoint> approxTraj;
    
    for (const auto& segment : slideSheet_.segments) {
        double startMD = segment.startMD;
        double endMD = segment.endMD;
        bool isSlide = segment.isSlide;
        
        // Find inclination values at segment boundaries
        double startInc = 0.0;
        double endInc = 0.0;
        
        // Get inclination at start point
        auto startIt = std::find_if(noSagTraj.begin(), noSagTraj.end(),
            [startMD](const TrajectoryPoint& p) { return std::abs(p.md - startMD) < 0.1; });
        
        if (startIt != noSagTraj.end()) {
            startInc = startIt->inclination;
        } else {
            // Interpolate
            auto next = std::find_if(noSagTraj.begin(), noSagTraj.end(),
                [startMD](const TrajectoryPoint& p) { return p.md > startMD; });
            
            if (next != noSagTraj.end() && next != noSagTraj.begin()) {
                auto prev = next - 1;
                double weight = (startMD - prev->md) / (next->md - prev->md);
                startInc = prev->inclination + weight * (next->inclination - prev->inclination);
            }
        }
        
        // Get inclination at end point
        auto endIt = std::find_if(noSagTraj.begin(), noSagTraj.end(),
            [endMD](const TrajectoryPoint& p) { return std::abs(p.md - endMD) < 0.1; });
        
        if (endIt != noSagTraj.end()) {
            endInc = endIt->inclination;
        } else {
            // Interpolate
            auto next = std::find_if(noSagTraj.begin(), noSagTraj.end(),
                [endMD](const TrajectoryPoint& p) { return p.md > endMD; });
            
            if (next != noSagTraj.end() && next != noSagTraj.begin()) {
                auto prev = next - 1;
                double weight = (endMD - prev->md) / (next->md - prev->md);
                endInc = prev->inclination + weight * (next->inclination - prev->inclination);
            }
        }
        
        // Create trajectory points along this segment
        double stepSize = 1.0; // 1m steps
        int numSteps = static_cast<int>((endMD - startMD) / stepSize) + 1;
        
        for (int i = 0; i < numSteps; i++) {
            double md = startMD + i * stepSize;
            if (md > endMD) md = endMD;
            
            // Linear interpolation of inclination
            double weight = (md - startMD) / (endMD - startMD);
            double inc = startInc + weight * (endInc - startInc);
            
            TrajectoryPoint point;
            point.md = md;
            point.inclination = inc;
            
            // Copy other values from the closest point in noSagTraj
            auto closest = std::min_element(noSagTraj.begin(), noSagTraj.end(),
                [md](const TrajectoryPoint& a, const TrajectoryPoint& b) {
                    return std::abs(a.md - md) < std::abs(b.md - md);
                });
            
            if (closest != noSagTraj.end()) {
                point.azimuth = closest->azimuth;
                point.tvd = closest->tvd;
                point.dls = closest->dls;
            }
            
            approxTraj.push_back(point);
        }
    }
    
    // Sort by measured depth
    std::sort(approxTraj.begin(), approxTraj.end(),
        [](const TrajectoryPoint& a, const TrajectoryPoint& b) { return a.md < b.md; });
    
    return approxTraj;
}

double TvdAssuranceWorkflow::calculateTrajDifference(
    const std::vector<TrajectoryPoint>& traj1,
    const std::vector<TrajectoryPoint>& traj2) {
    
    // Calculate RMS difference in inclination between two trajectories
    double sumSquaredDiff = 0.0;
    int count = 0;
    
    // For each point in traj1, find the nearest point in traj2
    for (const auto& point1 : traj1) {
        auto closest = std::min_element(traj2.begin(), traj2.end(),
            [&point1](const TrajectoryPoint& a, const TrajectoryPoint& b) {
                return std::abs(a.md - point1.md) < std::abs(b.md - point1.md);
            });
        
        if (closest != traj2.end() && std::abs(closest->md - point1.md) < 0.5) {
            double diff = point1.inclination - closest->inclination;
            sumSquaredDiff += diff * diff;
            count++;
        }
    }
    
    if (count == 0) return 1000.0; // Large value if no matching points
    
    return std::sqrt(sumSquaredDiff / count);
}