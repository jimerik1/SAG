# src/api/parsers.py
import numpy as np
from typing import List, Dict, Any

def parse_trajectory(trajectory_data: List[Dict[str, Any]]) -> np.ndarray:
    """
    Parse trajectory data from API input
    
    Returns a structured numpy array with trajectory points
    """
    trajectory = np.zeros(len(trajectory_data), 
                         dtype=[('md', float), ('inc', float), ('azi', float), 
                                ('tvd', float), ('dls', float)])
    
    for i, point in enumerate(trajectory_data):
        trajectory[i] = (
            point['md'],
            point['inclination'],
            point['azimuth'],
            point.get('tvd', 0.0),  # TVD might be calculated if not provided
            point.get('dls', 0.0)   # DLS might be calculated if not provided
        )
    
    # If DLS is not provided, calculate it
    if np.all(trajectory['dls'] == 0):
        # Calculate DLS between points
        # Implementation depends on your specific DLS calculation method
        pass
        
    return trajectory

def parse_bha_components(bha_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Parse BHA component data from API input
    """
    components = []
    distance_from_bit = 0.0
    
    # Sort components from bit upward
    # Assuming first component is the bit
    for component in bha_data:
        component_dict = {
            'name': component['name'],
            'length': float(component['length']),
            'outerDiameter': float(component['body_od']),
            'innerDiameter': float(component['body_id']),
            'weight': float(component['nominal_weight']),
            'material': component['grade'],
            'distanceFromBit': distance_from_bit,
            'isStabilizer': 'stabilizer' in component['name'].lower()
        }
        
        distance_from_bit += component_dict['length']
        components.append(component_dict)
    
    return components
