// src/calculators/fea_core/trajectory.cpp
#include "trajectory.h"
#include <algorithm>
#include <cmath>

double Trajectory::getCurvatureAtDepth(double md) const {
    // Find closest points
    if (points.empty()) {
        return 0.0;
    }
    
    // If md is beyond the trajectory, return the curvature of the last segment
    if (md > points.back().md) {
        return points.back().dls;
    }
    
    // If md is before the trajectory, return the curvature of the first segment
    if (md < points.front().md) {
        return points.front().dls;
    }
    
    // Find the trajectory points that bracket the requested depth
    auto it = std::lower_bound(points.begin(), points.end(), md, 
                               [](const TrajectoryPoint& point, double depth) {
                                   return point.md < depth;
                               });
    
    if (it == points.begin()) {
        return it->dls;
    }
    
    auto prev = std::prev(it);
    
    // Interpolate curvature (dogleg severity)
    double ratio = (md - prev->md) / (it->md - prev->md);
    return prev->dls * (1 - ratio) + it->dls * ratio;
}

TrajectoryPoint Trajectory::interpolateAtDepth(double md) const {
    // Create a new point with interpolated values
    TrajectoryPoint result;
    result.md = md;
    
    // Handle empty trajectory or out-of-bounds md
    if (points.empty()) {
        result.inclination = 0.0;
        result.azimuth = 0.0;
        result.tvd = md;
        result.dls = 0.0;
        return result;
    }
    
    // If md is beyond the trajectory, extrapolate from the last segment
    if (md > points.back().md) {
        result.inclination = points.back().inclination;
        result.azimuth = points.back().azimuth;
        // Simple extrapolation for TVD
        double lastIncRad = points.back().inclination * M_PI / 180.0;
        result.tvd = points.back().tvd + (md - points.back().md) * cos(lastIncRad);
        result.dls = points.back().dls;
        return result;
    }
    
    // If md is before the trajectory, extrapolate from the first segment
    if (md < points.front().md) {
        result.inclination = points.front().inclination;
        result.azimuth = points.front().azimuth;
        // Simple extrapolation for TVD
        double firstIncRad = points.front().inclination * M_PI / 180.0;
        result.tvd = points.front().tvd - (points.front().md - md) * cos(firstIncRad);
        result.dls = points.front().dls;
        return result;
    }
    
    // Find the trajectory points that bracket the requested depth
    auto it = std::lower_bound(points.begin(), points.end(), md, 
                               [](const TrajectoryPoint& point, double depth) {
                                   return point.md < depth;
                               });
    
    if (it == points.begin()) {
        return *it;
    }
    
    auto prev = std::prev(it);
    
    // Interpolate all values
    double ratio = (md - prev->md) / (it->md - prev->md);
    result.inclination = prev->inclination * (1 - ratio) + it->inclination * ratio;
    result.azimuth = prev->azimuth * (1 - ratio) + it->azimuth * ratio;
    result.tvd = prev->tvd * (1 - ratio) + it->tvd * ratio;
    result.dls = prev->dls * (1 - ratio) + it->dls * ratio;
    
    return result;
}