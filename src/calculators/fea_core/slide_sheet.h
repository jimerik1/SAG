// src/calculators/fea_core/slide_sheet.h
#pragma once
#include <vector>
#include <string>

struct SlideSegment {
    double startMD;      // Starting measured depth
    double endMD;        // Ending measured depth
    bool isSlide;        // true for slide, false for rotate
    double toolface;     // Toolface angle in degrees (for slides)
    double rpm;          // RPM (for rotary mode)
    double rop;          // Rate of penetration
};

class SlideSheet {
public:
    SlideSheet() = default;
    
    // Add a segment to the slide sheet
    void addSegment(double startMD, double endMD, bool isSlide, 
                    double toolface = 0.0, double rpm = 0.0, double rop = 0.0);
    
    // Get segments that overlap with a specific measured depth range
    std::vector<SlideSegment> getSegmentsInRange(double startMD, double endMD) const;
    
    // Check if a given measured depth is within a slide segment
    bool isInSlide(double md) const;
    
    // Get the toolface at a given measured depth (returns 0.0 if not in slide)
    double getToolfaceAt(double md) const;
    
    // Get the total slide percentage
    double getSlidePercentage() const;
    
    // Public member holding all segments
    std::vector<SlideSegment> segments;
};