// src/utils/file_readers_slide_sheet.cpp
#include "file_readers.h"
#include "../calculators/fea_core/slide_sheet.h"
#include <fstream>
#include <sstream>
#include <algorithm>

// Add this method to your FileReaders class in file_readers.h
/*
public:
    // Load slide sheet data from file
    static SlideSheet loadSlideSheetFromFile(const std::string& filename);
*/

SlideSheet FileReaders::loadSlideSheetFromFile(const std::string& filename) {
    SlideSheet slideSheet;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open slide sheet file: " + filename);
    }
    
    // Skip header line
    std::string line;
    std::getline(file, line);
    
    // Read data
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        // Parse start MD
        std::getline(ss, token, ',');
        std::replace(token.begin(), token.end(), ',', '.');
        double startMD = std::stod(token);
        
        // Parse end MD
        std::getline(ss, token, ',');
        std::replace(token.begin(), token.end(), ',', '.');
        double endMD = std::stod(token);
        
        // Parse drilling mode (slide/rotate)
        std::getline(ss, token, ',');
        bool isSlide = (token == "slide" || token == "Slide" || token == "SLIDE");
        
        // Parse toolface (if available)
        double toolface = 0.0;
        if (std::getline(ss, token, ',')) {
            std::replace(token.begin(), token.end(), ',', '.');
            if (!token.empty()) {
                toolface = std::stod(token);
            }
        }
        
        // Parse RPM (if available)
        double rpm = 0.0;
        if (std::getline(ss, token, ',')) {
            std::replace(token.begin(), token.end(), ',', '.');
            if (!token.empty()) {
                rpm = std::stod(token);
            }
        }
        
        // Parse ROP (if available)
        double rop = 0.0;
        if (std::getline(ss, token, ',')) {
            std::replace(token.begin(), token.end(), ',', '.');
            if (!token.empty()) {
                rop = std::stod(token);
            }
        }
        
        // Add segment to slide sheet
        slideSheet.addSegment(startMD, endMD, isSlide, toolface, rpm, rop);
    }
    
    return slideSheet;
}

// Implementation of the slide sheet methods
void SlideSheet::addSegment(double startMD, double endMD, bool isSlide, 
                            double toolface, double rpm, double rop) {
    SlideSegment segment;
    segment.startMD = startMD;
    segment.endMD = endMD;
    segment.isSlide = isSlide;
    segment.toolface = toolface;
    segment.rpm = rpm;
    segment.rop = rop;
    
    segments.push_back(segment);
    
    // Sort segments by startMD
    std::sort(segments.begin(), segments.end(), 
              [](const SlideSegment& a, const SlideSegment& b) {
                  return a.startMD < b.startMD;
              });
}

std::vector<SlideSegment> SlideSheet::getSegmentsInRange(double startMD, double endMD) const {
    std::vector<SlideSegment> result;
    
    for (const auto& segment : segments) {
        // Check if segment overlaps with the specified range
        if (segment.endMD >= startMD && segment.startMD <= endMD) {
            result.push_back(segment);
        }
    }
    
    return result;
}

bool SlideSheet::isInSlide(double md) const {
    for (const auto& segment : segments) {
        if (md >= segment.startMD && md <= segment.endMD && segment.isSlide) {
            return true;
        }
    }
    
    return false;
}

double SlideSheet::getToolfaceAt(double md) const {
    for (const auto& segment : segments) {
        if (md >= segment.startMD && md <= segment.endMD && segment.isSlide) {
            return segment.toolface;
        }
    }
    
    return 0.0;
}

double SlideSheet::getSlidePercentage() const {
    double totalLength = 0.0;
    double slideLength = 0.0;
    
    for (const auto& segment : segments) {
        double length = segment.endMD - segment.startMD;
        totalLength += length;
        
        if (segment.isSlide) {
            slideLength += length;
        }
    }
    
    return (totalLength > 0.0) ? (slideLength / totalLength * 100.0) : 0.0;
}