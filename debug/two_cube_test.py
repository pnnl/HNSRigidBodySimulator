import pybullet as p
import pybullet_data
import numpy as np
import math

def system_KE(meshes):
    tKE = 0
    rKE = 0
    for i in range(len(meshes)):
        mesh = meshes[i]
        mass = p.getDynamicsInfo(mesh, -1)[0]
        inertia = p.getDynamicsInfo(mesh, -1)[2][0]
        lin_velocity, ang_velocity = p.getBaseVelocity(mesh)
        tKE += 0.5 * mass * (lin_velocity[0] ** 2 + lin_velocity[1] ** 2 + lin_velocity[2] ** 2)
        rKE += 0.5 * inertia * (ang_velocity[0] ** 2 + ang_velocity[1] ** 2 + ang_velocity[2] ** 2)
    return tKE + rKE

def system_PE(meshes, epsilon, sigma, epsilon_min):
    f = open("pairs.out", "w")
    pe = 0.0
    num_meshes = len(meshes)
    min_dist = 100000
    
    for i in range(num_meshes):
        for j in range(i + 1, num_meshes):
            pos1, _ = p.getBasePositionAndOrientation(meshes[i])
            pos2, _ = p.getBasePositionAndOrientation(meshes[j])
            distance_vector = np.array(pos2) - np.array(pos1)
            distance = np.linalg.norm(distance_vector)
            distance = max(distance, epsilon_min)
            if distance < min_dist:
               min_dist = distance
            temp = 4*epsilon*((sigma/distance)**12 - (sigma/distance)**6)
            pe += temp
            f.write(str(distance) +  "  " + str(temp) + "\n")
    
    return pe, min_dist

# Calculates lennard-jones force between two objects
def lennard_jones_force(distance_vector, epsilon, sigma, epsilon_min):
    distance = math.sqrt(sum(v * v for v in distance_vector))
    distance = max(distance, epsilon_min)  # Ensure distance is not too close to zero
    force_magnitude = -4*epsilon*(12 * (sigma**12/distance**13) - 6*(sigma**6/distance**7))
    #force_magnitude = (-48*epsilon/distance)*((sigma**12/distance**12) - (sigma**6/distance**6))
    #force_magnitude = -4 * epsilon * ((6 * sigma ** 6 / distance ** 7) - (12 * sigma ** 12/distance ** 13))
    force = [force_magnitude * (v/distance) for v in distance_vector]
    return force

def scale_velocity(body_a, body_b, scale_factor):

    a_lin_vel = p.getBaseVelocity(body_a)[0]
    #a_ang_vel = p.getDynamicsInfo(body_a)[1]
    b_lin_vel = p.getBaseVelocity(body_b)[0]
    #b_ang_vel = p.getDynamicsInfo(body_b)[1]

    p.resetBaseVelocity(body_a, linearVelocity = np.array(a_lin_vel) / math.sqrt(scale_factor))
    p.resetBaseVelocity(body_b, linearVelocity = np.array(b_lin_vel) / math.sqrt(scale_factor))
    

# Connect to the physics server
p.connect(p.GUI)  # Use p.DIRECT for headless mode

# Set additional search path for URDF data
p.setAdditionalSearchPath(pybullet_data.getDataPath())

# Load the URDF file for the cube
cube_urdf_path = "general.urdf"

# Create the first cube and set its initial position
cube1_position = [0, 0, 2]  # Place the first cube 5 units above the ground
cube1_orientation = [0, 0, 0, 1]  # No angular offset
cube1_id = p.loadURDF(cube_urdf_path, basePosition=cube1_position, baseOrientation=cube1_orientation)
p.resetBaseVelocity(cube1_id, linearVelocity=[0, 0, 0], angularVelocity=[0, 0, 0])  # Move towards the second cube
p.changeDynamics(cube1_id, -1, mass=1, restitution=1.0, lateralFriction=0.0, spinningFriction=0.0, rollingFriction=0.0, linearDamping=0.0, angularDamping=0.0)

# Create the second cube and set its initial position
cube2_position = [0, 0, -2]  # Place the second cube 5 units below the ground
cube2_orientation = [0, 0, 0, 1]  # No angular offset
cube2_id = p.loadURDF(cube_urdf_path, basePosition=cube2_position, baseOrientation=cube2_orientation)
p.resetBaseVelocity(cube2_id, linearVelocity=[0, 0, 0], angularVelocity=[0, 0, 0])  # Move towards the second cube
p.changeDynamics(cube2_id, -1, mass=1, restitution=1.0, lateralFriction=0.0, spinningFriction=0.0, rollingFriction=0.0, linearDamping=0.0, angularDamping=0.0)


meshes = [cube1_id, cube2_id]


e = open("e.out", "w")
c = open("col.out", "w")
v = open("vels.out", "w")
dists = open("two_dist.out", "w")

epsilon = 1000.0
sigma = 0.8908987181403393
e_min = 0.0

kinetic = system_KE(meshes)
energy_prev = kinetic
energy_current = energy_prev
contacts_per_step = []

p.setTimeStep(1/500)
# Simulation loop
for t in range(1000):
    num_meshes = 2
    # Calculates and applies LJ forces to each mesh
    forces = np.zeros((num_meshes, 3))

    for i in range(num_meshes):
        for j in range(i + 1, num_meshes): 
            pos1, _ = p.getBasePositionAndOrientation(meshes[i])
            pos2, _ = p.getBasePositionAndOrientation(meshes[j])
                
            distance_vector = [pos2[k] - pos1[k] for k in range(3)]  
            force = lennard_jones_force(distance_vector, epsilon, sigma, e_min)
            forces[i] = np.add(forces[i], force)  
            forces[j] = np.subtract(forces[j], force) 

    for i in range(num_meshes):
        com_position, _ = p.getBasePositionAndOrientation(meshes[i])
        p.applyExternalForce(meshes[i], -1, forces[i], com_position, p.WORLD_FRAME)

    kinetic = system_KE(meshes)
    potential, min_dist = system_PE(meshes, epsilon, sigma, e_min)
    pos1, _ = p.getBasePositionAndOrientation(meshes[0])
    pos2, _ = p.getBasePositionAndOrientation(meshes[1])
    distance_vector = np.array(pos2) - np.array(pos1)
    distance = np.linalg.norm(distance_vector)
    e.write(f"{t} {potential} {kinetic} {potential + kinetic} {distance}\n")
    lin1, ang1 = p.getBaseVelocity(meshes[0])
    lin2, ang2 = p.getBaseVelocity(meshes[1])
    lin1 = np.linalg.norm(np.array(lin1))
    ang1 = np.linalg.norm(np.array(ang1))
    lin2 = np.linalg.norm(np.array(lin2))
    ang2 = np.linalg.norm(np.array(ang2))
    v.write(f"{t} {lin1} {ang1} {lin2} {ang2}\n")
    contacts = p.getContactPoints()
    if contacts:
           c.write(f"{t} {potential}\n")
    
    dists.write(f"{distance}\n")
    energy_prev = energy_current
    p.stepSimulation()


# Disconnect from the physics server
p.disconnect()
