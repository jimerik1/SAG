// src/calculators/fea_core/trajectory.h
#pragma once
#include <vector>

struct TrajectoryPoint {
    double md;         // Measured depth (m)
    double inclination; // Inclination angle (degrees)
    double azimuth;    // Azimuth angle (degrees)
    double tvd;        // True vertical depth (m)
    double dls;        // Dogleg severity (deg/30m)
};

class Trajectory {
public:
    std::vector<TrajectoryPoint> points;
    
    // Calculate wellbore curvature at a specific measured depth
    double getCurvatureAtDepth(double md) const;
    
    // Interpolate trajectory at arbitrary measured depth
    TrajectoryPoint interpolateAtDepth(double md) const;
};
