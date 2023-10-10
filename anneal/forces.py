import math
import numpy as np

# Calculates lennard-jones force between two objects
def lennard_jones_force(distance_vector, epsilon, sigma, epsilon_min):
    distance = math.sqrt(sum(v * v for v in distance_vector))
    distance = max(distance, epsilon_min)  # Ensure distance is not too close to zero
    force_magnitude = -4*epsilon*(12 * (sigma**12/distance**13) - 6*(sigma**6/distance**7))
    force = [force_magnitude * (v/distance) for v in distance_vector]
    return force

def r_cubed_force(distance_vector, k):
    distance = np.linalg.norm(np.array(distance_vector))
    force_magnitude = k / distance ** 3
    force = [force_magnitude * (v/distance) for v in distance_vector]
    return force

def r_squared_force(distance_vector, k):
    distance = np.linalg.norm(np.array(distance_vector))
    force_magnitude = k / distance ** 2
    force = [force_magnitude * (v/distance) for v in distance_vector]
    return force

def over_r_force(distance_vector, k):
    distance = np.linalg.norm(np.array(distance_vector))
    force_magnitude = k / distance
    force = [force_magnitude * (v/distance) for v in distance_vector]
    return force

