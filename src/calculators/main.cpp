// src/calculators/main.cpp
#include <iostream>
#include <fstream>
#include <filesystem>
#include "fea_core/fea_solver.h"
#include "fea_core/fea_results.h"
#include "../utils/file_readers.h"

// Function to generate simple HTML visualization
void generateVisualization(const FEAResults& results, const std::string& outputFile) {
    std::ofstream htmlFile(outputFile);
    htmlFile << "<!DOCTYPE html>\n";
    htmlFile << "<html lang=\"en\">\n";
    htmlFile << "<head>\n";
    htmlFile << "    <meta charset=\"UTF-8\">\n";
    htmlFile << "    <title>BHA Sag Analysis Results</title>\n";
    htmlFile << "    <script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>\n";
    htmlFile << "    <style>\n";
    htmlFile << "        .chart-container { width: 80%; margin: 20px auto; }\n";
    htmlFile << "        .results-container { width: 80%; margin: 20px auto; font-family: Arial, sans-serif; }\n";
    htmlFile << "        table { width: 100%; border-collapse: collapse; }\n";
    htmlFile << "        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n";
    htmlFile << "        th { background-color: #f2f2f2; }\n";
    htmlFile << "    </style>\n";
    htmlFile << "</head>\n";
    htmlFile << "<body>\n";
    
    // Summary results
    htmlFile << "    <div class=\"results-container\">\n";
    htmlFile << "        <h1>BHA Sag Analysis Results</h1>\n";
    htmlFile << "        <table>\n";
    htmlFile << "            <tr><th>Parameter</th><th>Value</th></tr>\n";
    htmlFile << "            <tr><td>Bit Side Force</td><td>" << results.bitSideForce << " kg</td></tr>\n";
    htmlFile << "            <tr><td>Bend</td><td>" << results.bend << " degs</td></tr>\n";
    htmlFile << "            <tr><td>WOB</td><td>" << results.wob << " kg</td></tr>\n";
    htmlFile << "            <tr><td>Sensor Position</td><td>" << results.sensorPosition << " m</td></tr>\n";
    htmlFile << "            <tr><td>Hole Curvature</td><td>" << results.holeCurvature << " deg/100ft</td></tr>\n";
    htmlFile << "            <tr><td>Sag Correction</td><td>" << results.sagCorrection << " deg</td></tr>\n";
    htmlFile << "        </table>\n";
    
    // Stabilizer forces
    htmlFile << "        <h2>Stabilizer Forces</h2>\n";
    htmlFile << "        <table>\n";
    htmlFile << "            <tr><th>Stabilizer</th><th>Position</th><th>Side Force</th></tr>\n";
    for (const auto& stab : results.stabilizerForces) {
        htmlFile << "            <tr><td>" << stab.name << "</td><td>" << stab.position 
                << " m</td><td>" << stab.sideForce << " kg</td></tr>\n";
    }
    htmlFile << "        </table>\n";
    htmlFile << "    </div>\n";
    
    // Charts for the profiles
    htmlFile << "    <div class=\"chart-container\">\n";
    htmlFile << "        <canvas id=\"deflectionChart\"></canvas>\n";
    htmlFile << "    </div>\n";
    
    htmlFile << "    <div class=\"chart-container\">\n";
    htmlFile << "        <canvas id=\"slopeChart\"></canvas>\n";
    htmlFile << "    </div>\n";
    
    htmlFile << "    <div class=\"chart-container\">\n";
    htmlFile << "        <canvas id=\"shearChart\"></canvas>\n";
    htmlFile << "    </div>\n";
    
    htmlFile << "    <div class=\"chart-container\">\n";
    htmlFile << "        <canvas id=\"momentChart\"></canvas>\n";
    htmlFile << "    </div>\n";
    
    // JavaScript for charts
    htmlFile << "    <script>\n";
    
    // Function to prepare data
    htmlFile << "        function prepareData(data) {\n";
    htmlFile << "            return {\n";
    htmlFile << "                labels: data.map(point => point[0]),\n";
    htmlFile << "                datasets: [{\n";
    htmlFile << "                    label: 'Data',\n";
    htmlFile << "                    data: data.map(point => point[1]),\n";
    htmlFile << "                    borderColor: 'rgb(75, 192, 192)',\n";
    htmlFile << "                    tension: 0.1\n";
    htmlFile << "                }]\n";
    htmlFile << "            };\n";
    htmlFile << "        }\n";
    
    // Deflection data
    htmlFile << "        const deflectionData = [\n";
    for (const auto& point : results.deflectionProfile) {
        htmlFile << "            [" << point.first << ", " << point.second << "],\n";
    }
    htmlFile << "        ];\n";
    
    // Slope data
    htmlFile << "        const slopeData = [\n";
    for (const auto& point : results.slopeProfile) {
        htmlFile << "            [" << point.first << ", " << point.second << "],\n";
    }
    htmlFile << "        ];\n";
    
    // Shear force data
    htmlFile << "        const shearData = [\n";
    for (const auto& point : results.shearForceProfile) {
        htmlFile << "            [" << point.first << ", " << point.second << "],\n";
    }
    htmlFile << "        ];\n";
    
    // Bending moment data
    htmlFile << "        const momentData = [\n";
    for (const auto& point : results.bendingMomentProfile) {
        htmlFile << "            [" << point.first << ", " << point.second << "],\n";
    }
    htmlFile << "        ];\n";
    
    // Create charts
    htmlFile << "        // Deflection Chart\n";
    htmlFile << "        new Chart(document.getElementById('deflectionChart'), {\n";
    htmlFile << "            type: 'line',\n";
    htmlFile << "            data: prepareData(deflectionData),\n";
    htmlFile << "            options: {\n";
    htmlFile << "                scales: {\n";
    htmlFile << "                    y: { reverse: true },\n";
    htmlFile << "                    x: { title: { display: true, text: 'Distance from Bit (m)' } }\n";
    htmlFile << "                },\n";
    htmlFile << "                plugins: { title: { display: true, text: 'BHA Deflected Shape' } }\n";
    htmlFile << "            }\n";
    htmlFile << "        });\n";
    
    htmlFile << "        // Slope Chart\n";
    htmlFile << "        new Chart(document.getElementById('slopeChart'), {\n";
    htmlFile << "            type: 'line',\n";
    htmlFile << "            data: prepareData(slopeData),\n";
    htmlFile << "            options: {\n";
    htmlFile << "                plugins: { title: { display: true, text: 'BHA Slope (degrees)' } },\n";
    htmlFile << "                scales: { x: { title: { display: true, text: 'Distance from Bit (m)' } } }\n";
    htmlFile << "            }\n";
    htmlFile << "        });\n";
    
    htmlFile << "        // Shear Force Chart\n";
    htmlFile << "        new Chart(document.getElementById('shearChart'), {\n";
    htmlFile << "            type: 'line',\n";
    htmlFile << "            data: prepareData(shearData),\n";
    htmlFile << "            options: {\n";
    htmlFile << "                plugins: { title: { display: true, text: 'Shear Force' } },\n";
    htmlFile << "                scales: { x: { title: { display: true, text: 'Distance from Bit (m)' } } }\n";
    htmlFile << "            }\n";
    htmlFile << "        });\n";
    
    htmlFile << "        // Bending Moment Chart\n";
    htmlFile << "        new Chart(document.getElementById('momentChart'), {\n";
    htmlFile << "            type: 'line',\n";
    htmlFile << "            data: prepareData(momentData),\n";
    htmlFile << "            options: {\n";
    htmlFile << "                plugins: { title: { display: true, text: 'Bending Moment' } },\n";
    htmlFile << "                scales: { x: { title: { display: true, text: 'Distance from Bit (m)' } } }\n";
    htmlFile << "            }\n";
    htmlFile << "        });\n";
    
    htmlFile << "    </script>\n";
    htmlFile << "</body>\n";
    htmlFile << "</html>\n";
    
    htmlFile.close();
}

int main() {
    try {
        // Load data from files
        std::string trajectoryFile = "calculators/data/trajectory.txt";
        std::string bhaComponentsFile = "calculators/data/bha_components.txt";
        std::string materialPropertiesFile = "calculators/data/material_properties.txt";
        
        std::cout << "Loading trajectory data..." << std::endl;
        Trajectory trajectory = FileReaders::loadTrajectoryFromFile(trajectoryFile);
        std::cout << "Loaded " << trajectory.points.size() << " trajectory points." << std::endl;
        
        std::cout << "Loading BHA components..." << std::endl;
        std::vector<BHAComponent> components = FileReaders::loadBHAComponentsFromFile(bhaComponentsFile);
        std::cout << "Loaded " << components.size() << " BHA components." << std::endl;
        
        std::cout << "Loading material properties..." << std::endl;
        std::vector<MaterialProperties> materials = FileReaders::loadMaterialPropertiesFromFile(materialPropertiesFile);
        std::cout << "Loaded " << materials.size() << " material definitions." << std::endl;
        
        // Create FEA solver
        double elementSize = 0.1;  // 10cm elements
        double mudWeight = 12.0;   // 12 ppg mud weight
        double wob = 20000.0;      // 20k lbs WOB
        
        std::cout << "Creating FEA solver..." << std::endl;
        FEASolver solver(components, materials, trajectory, elementSize, mudWeight, wob);
        
        // Set up and solve FEA system
        std::cout << "Setting up FEA system..." << std::endl;
        solver.setupSystem();
        
        std::cout << "Solving FEA system..." << std::endl;
        solver.solve();
        
        // Get comprehensive results
        FEAResults results = solver.getResults();
        
        // Print summary to console
        results.printSummary();
        
        // Create output directory if it doesn't exist
        std::filesystem::create_directories("calculators/results");
        
        // Write detailed results to files
        results.writeProfilesToFiles("calculators/results");
        
        // Generate visualization
        generateVisualization(results, "calculators/results/visualization.html");
        
        std::cout << "Calculation complete. Results written to 'calculators/results/' directory." << std::endl;
        std::cout << "Open visualization.html in a web browser to view interactive charts." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}