import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

trajectories = np.loadtxt('output.out')

n = 0 # particle index
n_steps = 4200
n_particles = 40

trajectory_n = np.zeros((n_steps, 7))

for i in range(n_steps):
    trajectory_n[i] = trajectories[i * n_particles + n]

pos_x = np.zeros(4200)
pos_y = np.zeros(4200)
pos_z = np.zeros(4200)

for i in range(4200):
    pos_x[i] = trajectory_n[i,0]
    pos_y[i] = trajectory_n[i,1]
    pos_z[i] = trajectory_n[i,2]

fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")

ax.plot_wireframe(pos_x,pos_y,pos_z)

plt.show()
# plt.plot(pos_x)
# plt.show()

