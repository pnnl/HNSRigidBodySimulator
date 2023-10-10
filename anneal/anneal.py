import math
import numpy as np
import random as r
import pybullet as p
import pybullet_data
import time
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import struct
import forces as f

def main():
    create_environment()

    inp = open("input.in", "r") # input file
    out = open("output.out", "w") #writes positions and quaternions
    ener = open("energy.out", "w") #writes energies. Plot with script "p1"
    dists = open("distances.out", "w") #writes the minimum distance each time-step

    num_meshes = int(inp.readline())
    meshes = np.zeros(num_meshes, dtype=int) 
    path = inp.readline().strip()
    box_size = float(inp.readline()) # box dimensions are [-box_size, box_size] in all three dimensions
    scale = float(inp.readline().strip())

    # Creates meshes at random positions in box
    for i in range(num_meshes):
        meshes[i] = create_mesh(mesh_path=path, position=[r.uniform(-box_size, box_size), r.uniform(-box_size, box_size), r.uniform(-box_size, box_size)], scale=scale)

    # Sets Lennard-Jones parameters
    params = inp.readline().split()
    epsilon = float(params[0])
    sigma = float(params[1])

    # Initializes structure factor calculations
    num_points = int(inp.readline())  # Number of points in each dimension for position delta function
    structure_factor = np.zeros((num_points, num_points, num_points))
    kx = np.zeros(num_points)
    ky = np.zeros(num_points)
    kz = np.zeros(num_points)
    count = 0
    ray_size = float(inp.readline())

    timesteps = int(inp.readline())
    freq = int(inp.readline()) # frequency at which to collect structure factor information

    kinetic, _, _ = system_KE(meshes)
    energy_prev = kinetic
    energy_current = energy_prev
    contacts_per_step = []

    prev_vels = [] # previous velocities, used to resest velocitys
    for mesh in meshes:
        prev_vels.append(p.getBaseVelocity(mesh)) 
    
    all_ke = []
    all_ke.append(system_KE(meshes)[0])

    p.setTimeStep(1/10000)
    for t in range(timesteps):
        forces = np.zeros((num_meshes, 3))

        # Calculates and applies LJ forces to each mesh
        for i in range(num_meshes):
            for j in range(i + 1, num_meshes): 
                pos1, _ = p.getBasePositionAndOrientation(meshes[i])
                pos2, _ = p.getBasePositionAndOrientation(meshes[j])
                distance_vector = [pos2[k] - pos1[k] for k in range(3)]  
                force = f.lennard_jones_force(distance_vector, epsilon, sigma, 1e-6)
                forces[i] = np.add(forces[i], force)  
                forces[j] = np.subtract(forces[j], force) 
        
        for i in range(num_meshes):
            com_position, _ = p.getBasePositionAndOrientation(meshes[i])
            p.applyExternalForce(meshes[i], -1, forces[i], com_position, p.WORLD_FRAME)

        # Enforces periodic boundary conditions (PBC) on the mesh positions
        box = [[-box_size, box_size], [-box_size, box_size], [-box_size, box_size]]
        for i in range(num_meshes):
            pos, orn = p.getBasePositionAndOrientation(meshes[i])
            pos = pbc(np.array(pos), box)
            velocity = p.getBaseVelocity(meshes[i])
            linear_vel = velocity[0]  
            angular_vel = velocity[1] 
            p.resetBasePositionAndOrientation(meshes[i], pos, orn) # Resets position
            p.resetBaseVelocity(meshes[i], linearVelocity=linear_vel, angularVelocity=angular_vel) # Resets velocity

        # Uncomment to check min-distance
        #minDist = checkDist(num_meshes, meshes)
        #dists.write(f"{minDist}\n")

        # scale velocities to conserve energies
        contacts_per_step.append(p.getContactPoints())
        recent_contacts = contacts_per_step[-1]
        bodies_to_update = []
        kinetic, _, _ = system_KE(meshes)
        energy_current = kinetic
        if energy_current > energy_prev:
            num_collisions = 0
            if recent_contacts != None:
                for contact in recent_contacts:
                    body_pair = [contact[1], contact[2]]
                    if body_pair not in bodies_to_update:
                        bodies_to_update.append(body_pair)
                        num_collisions += 1

            if num_collisions > 0:
                for pair in bodies_to_update:
                    scale_velocity(pair[0], pair[1], prev_vels)

        # Calls structure factor calculation
        if (t % freq == 0):
            structure_factor_temp, kx_temp, ky_temp, kz_temp = calculate_structure_factor(box_size, num_points, ray_size)
            structure_factor += structure_factor_temp
            kx += kx_temp
            ky += ky_temp
            kz += kz_temp
            count += 1

        # store all velocities in list
        prev_vels = []
        for mesh in meshes:
            prev_vels.append(p.getBaseVelocity(mesh)) 
        
        energy_prev = energy_current

        # Write positions and orientations of all objects to a file
        for mesh in meshes:
            pos, orn = p.getBasePositionAndOrientation(mesh)
            out.write(f"{pos[0]}  {pos[1]}  {pos[2]} {orn[0]}  {orn[1]}  {orn[2]}  {orn[3]}\n")

        #Resets velocities to zero if KE begins decreasing
        all_ke.append(system_KE(meshes)[0])
        if all_ke[t-1] > all_ke[t]:
            for i in range (num_meshes):
                p.resetBaseVelocity(meshes[i], linearVelocity=[0,0,0], angularVelocity=[0,0,0])

        ke,_,_= system_KE(meshes)
        pe,_ = system_PE(meshes, epsilon, sigma, epsilon_min=1.0)
        te = ke + pe

        ener.write(f"{t}  {pe}  {ke}  {te}\n")

        p.stepSimulation()
    
    # Calculates average values
    structure_factor /= count
    kx /= count
    ky /= count
    kz /= count

    # Plot the magnitutde of k vs. S(k)
    plot = open("s_of_k.out", "w") # s_of_k output file
    k_mag = []
    s_of_k = []
    for x in range(num_points):
       for y in range(num_points):
          for z in range(num_points):
             k_mag.append(math.sqrt(kx[x] ** 2 + ky[y] ** 2 + kz[z] ** 2))
             s_of_k.append(structure_factor[x][y][z])

             plot.write(str(k_mag[-1]) + "  " + str(s_of_k[-1]) + "\n")

    plt.scatter(k_mag, s_of_k)
    plt.show()
    
    p.disconnect()


def scale_velocity(body_a, body_b, prev_vels):

    a_lin_vel = p.getBaseVelocity(body_a)[0]
    a_ang_vel = p.getBaseVelocity(body_a)[1]
    b_lin_vel = p.getBaseVelocity(body_b)[0]
    b_ang_vel = p.getBaseVelocity(body_b)[1]

    curr_KE, _, _ = system_KE([body_a, body_b])

    tKE = 0
    rKE = 0
    mass = p.getDynamicsInfo(body_a, -1)[0]
    inertia = p.getDynamicsInfo(body_a, -1)[2][0]
    lin_velocity_a, ang_velocity_a = prev_vels[body_a]
    lin_velocity_b, ang_velocity_b = prev_vels[body_b]
    tKE += 0.5 * mass * (lin_velocity_a[0] ** 2 + lin_velocity_a[1] ** 2 + lin_velocity_a[2] ** 2)
    tKE += 0.5 * mass * (lin_velocity_b[0] ** 2 + lin_velocity_b[1] ** 2 + lin_velocity_b[2] ** 2)
    rKE += 0.5 * inertia * (ang_velocity_a[0] ** 2 + ang_velocity_a[1] ** 2 + ang_velocity_a[2] ** 2)
    rKE += 0.5 * inertia * (ang_velocity_b[0] ** 2 + ang_velocity_b[1] ** 2 + ang_velocity_b[2] ** 2)

    prev_KE = tKE + rKE

    scale_factor = curr_KE / prev_KE

    p.resetBaseVelocity(body_a, linearVelocity = np.array(a_lin_vel) / math.sqrt(scale_factor), angularVelocity = np.array(a_ang_vel) / math.sqrt(scale_factor))
    p.resetBaseVelocity(body_b, linearVelocity = np.array(b_lin_vel) / math.sqrt(scale_factor), angularVelocity = np.array(b_ang_vel) / math.sqrt(scale_factor))

# Initializes the PyBullet simulation environment
def create_environment():
    p.connect(p.GUI) # p.GUI for visualization, p.DIRECT for no visualiztion
    p.setAdditionalSearchPath(pybullet_data.getDataPath())

# Creates a mesh given an object file, position, mass, and scale
def create_mesh(mesh_path, position, mass=2, scale=1):
    mesh = p.loadURDF("general.urdf", basePosition= position, baseOrientation=[0, 0, 0, 1], useFixedBase=0, globalScaling=scale)

    p.changeDynamics(mesh, -1, mass=mass, restitution=1.0, lateralFriction=0.0, spinningFriction=0.0, rollingFriction=0.0, linearDamping=0.0, angularDamping=0.0)
    
    ndof = p.getNumJoints(mesh) + 6
    
    # Generate initial velocities using the gaussian function
    initial_temp = 500
    initial_velocities = gaussian(initial_temp, ndof, 1)
    linear_vel = initial_velocities[:3]
    angular_vel = initial_velocities[3:]  
    p.resetBaseVelocity(mesh, linearVelocity=linear_vel, angularVelocity=angular_vel)

    return mesh

def rethermalize(temp, meshes):
    for mesh in meshes:
            ndof = p.getNumJoints(mesh) + 6
    
            # Generate initial velocities using the gaussian function
            initial_velocities = gaussian(temp, ndof, 1)
            linear_vel = initial_velocities[:3]
            angular_vel = initial_velocities[3:]  
            p.resetBaseVelocity(mesh, linearVelocity=linear_vel, angularVelocity=angular_vel) 

def system_KE(meshes):
    tKE = 0
    rKE = 0
    lin_total = 0
    ang_total = 0
    for i in range(len(meshes)):
        mesh = meshes[i]
        mass = p.getDynamicsInfo(mesh, -1)[0]
        inertia = p.getDynamicsInfo(mesh, -1)[2][0]
        lin_velocity, ang_velocity = p.getBaseVelocity(mesh)
        tKE += 0.5 * mass * (lin_velocity[0] ** 2 + lin_velocity[1] ** 2 + lin_velocity[2] ** 2)
        rKE += 0.5 * inertia * (ang_velocity[0] ** 2 + ang_velocity[1] ** 2 + ang_velocity[2] ** 2)
        lin_total += np.linalg.norm(np.array(lin_velocity))
        ang_total += np.linalg.norm(np.array(ang_velocity))
    return tKE + rKE, lin_total, ang_total#[0]

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

# Generates initial particle velocities with temperature T according to a Maxwell-Boltzmann distribution
def gaussian(temp, natoms, ndims):
 sqrtT = temp**0.5
 seed=7654321
 np.random.seed(seed)
 vels=np.zeros((natoms,ndims))
 for ii in range(natoms):
   for jj in range(ndims):
     vels[ii,jj] = sqrtT * r.gauss(0.0, 1.0)
 sumv = np.zeros((ndims))
 for ii in range(natoms):
   for jj in range(ndims):
     sumv[jj] = sumv[jj] + vels[ii,jj]
 sumv = sumv / natoms
 for ii in range(natoms):
   for jj in range(ndims):
     vels[ii,jj] = vels[ii,jj] - sumv[jj]
 return vels

# Sets position based on periodic boundary conditions
def pbc(pos, box):
    pos = pos.reshape((3,1))
    natoms, ndims = pos.shape
    for j in range(ndims):
        boxl = box[j][1] - box[j][0]
        for i in range(natoms):
            pos[i, j] = pos[i, j] - boxl * int(np.round(pos[i, j] / boxl))
    return pos

# Calculates structure factor of system
def calculate_structure_factor(box_size, num_points, ray_size):
    # Creates position delta function f(r)
    lattice_size = box_size / num_points

    fr = np.zeros((num_points, num_points, num_points))
    for x in range(num_points):
       for y in range(num_points):
          for z in range(num_points):
             point = (x * lattice_size - 10, y * lattice_size - 10, z * lattice_size - 10)
             fr[x][y][z] = inside(point, ray_size)

    # Used to test if f(r) function is working correctly
    fofr = open("fr.out", "w")
    for x in range(num_points):
       for y in range(num_points):
          for z in range(num_points):
             if fr[x][y][z] == 1:
                fofr.write(str(x) + "  " + str(y) + "  " + str(z) + "\n")
                                  
    # Performs 3D Fourier transform on f(r), then magnitude squared
    structure_factor = np.abs(fft.fftn(fr))**2

    # Calculate the magnitudes of spatial frequencies
    kx = fft.fftfreq(num_points, d=(2 * box_size) / num_points)
    ky = fft.fftfreq(num_points, d=(2 * box_size) / num_points)
    kz = fft.fftfreq(num_points, d=(2 * box_size) / num_points)

    return structure_factor, kx, ky, kz

# Checks if a point is inside an object
def inside(point, ray_extend):
    # Cast a ray from the point outward
    ray_from = point
    ray_to = [point[0], point[1], point[2] - ray_extend] 

    # Perform the ray test
    results = p.rayTest(rayFromPosition=ray_from, rayToPosition=ray_to)
    if results[0][0] != -1:
        return 1
    return 0

def checkDist(num_meshes, meshes):
    min_dist = 100000
    for i in range(num_meshes):
        for j in range(num_meshes):
            if i != j:
                pos_1, _ = p.getBasePositionAndOrientation(meshes[i])
                pos_2, _ = p.getBasePositionAndOrientation(meshes[j])
                dist_vector = np.array(pos_1) - np.array(pos_2)
                dist = np.linalg.norm(dist_vector)
                if dist < min_dist:
                    min_dist = dist
    return min_dist

if __name__ == "__main__":
    main()
