// src/calculators/fea_core/fea_results.cpp
#include "fea_results.h"
#include <iostream>
#include <fstream>
#include <iomanip>

void FEAResults::printSummary() const {
    std::cout << "======== BHA SAG CALCULATION RESULTS ========" << std::endl;
    std::cout << "Bit Side Force = " << bitSideForce << " kg" << std::endl;
    std::cout << "Bend = " << bend << " degs" << std::endl;
    std::cout << "WOB = " << wob << " kg" << std::endl;
    std::cout << "Sensor Position = " << sensorPosition << " m" << std::endl;
    std::cout << "Hole Curvature = " << holeCurvature << " deg/100ft" << std::endl;
    std::cout << "Sag Correction = " << sagCorrection << " deg" << std::endl;
    
    std::cout << std::endl << "Stabilizer Forces:" << std::endl;
    for (const auto& stab : stabilizerForces) {
        std::cout << stab.name << " at " << stab.position << "m: Side Force = " 
                  << stab.sideForce << " kg" << std::endl;
    }
    
    std::cout << std::endl << "Number of data points in profiles:" << std::endl;
    std::cout << "Deflection: " << deflectionProfile.size() << std::endl;
    std::cout << "Shear Force: " << shearForceProfile.size() << std::endl;
    std::cout << "Bending Moment: " << bendingMomentProfile.size() << std::endl;
    std::cout << "Slope: " << slopeProfile.size() << std::endl;
}

void FEAResults::writeProfilesToFiles(const std::string& directory) const {
    // Write deflection profile
    std::ofstream deflFile(directory + "/deflection_profile.txt");
    deflFile << "Position(m)\tDeflection(m)" << std::endl;
    for (const auto& point : deflectionProfile) {
        deflFile << point.first << "\t" << point.second << std::endl;
    }
    deflFile.close();
    
    // Write shear force profile
    std::ofstream shearFile(directory + "/shear_profile.txt");
    shearFile << "Position(m)\tShear_Force(N)" << std::endl;
    for (const auto& point : shearForceProfile) {
        shearFile << point.first << "\t" << point.second << std::endl;
    }
    shearFile.close();
    
    // Write bending moment profile
    std::ofstream momentFile(directory + "/moment_profile.txt");
    momentFile << "Position(m)\tBending_Moment(Nm)" << std::endl;
    for (const auto& point : bendingMomentProfile) {
        momentFile << point.first << "\t" << point.second << std::endl;
    }
    momentFile.close();
    
    // Write slope profile
    std::ofstream slopeFile(directory + "/slope_profile.txt");
    slopeFile << "Position(m)\tSlope(deg)" << std::endl;
    for (const auto& point : slopeProfile) {
        slopeFile << point.first << "\t" << point.second << std::endl;
    }
    slopeFile.close();
}