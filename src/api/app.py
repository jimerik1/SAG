# src/api/app.py
from flask import Flask, request, jsonify
import ctypes
import numpy as np
import os
from .parsers import parse_trajectory, parse_bha_components

app = Flask(__name__)

# Load the C++ library
lib_path = os.path.join(os.path.dirname(__file__), '../../build/libfea_solver.so')
fea_lib = ctypes.CDLL(lib_path)

# Define the API endpoint
@app.route('/calculate-sag', methods=['POST'])
def calculate_sag():
    # Get input data
    data = request.json
    
    # Validate input
    if not all(k in data for k in ['trajectory', 'bha_components', 'mud_weight']):
        return jsonify({'error': 'Missing required input fields'}), 400
    
    # Parse input data
    trajectory = parse_trajectory(data['trajectory'])
    bha_components = parse_bha_components(data['bha_components'])
    mud_weight = float(data['mud_weight'])
    element_size = float(data.get('element_size', 0.1))  # Default 10cm
    wob = float(data.get('weight_on_bit', 0.0))
    
    # Find MWD component
    mwd_distance = None
    for component in bha_components:
        if 'mwd' in component['name'].lower():
            mwd_distance = component['distanceFromBit']
            break
    
    if mwd_distance is None:
        return jsonify({'error': 'MWD component not found in BHA'}), 400
    
    try:
        # Call the C++ library to calculate sag
        # This would be replaced with your actual C++ binding code
        
        # Mock result for demonstration
        sag_correction = 0.25  # degrees
        
        # Prepare response
        response = {
            'sag_correction': sag_correction,
            'mwd_position': mwd_distance,
            'computation_details': {
                'element_size': element_size,
                'number_of_elements': len(trajectory) * 10  # Approximate
            }
        }
        
        return jsonify(response)
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)