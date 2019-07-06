import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



scater = np.loadtxt("Gene395_071_400_Ce2temperature.txt", dtype='float')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(len(scater[: , 0])):
    ax.scatter(scater[i , 0], scater[i , 1], scater[i , 2], c=( i / len(scater[: , 0]), 0.0, 1 - i/len(scater[: , 0])) )
ax.set_xlabel('Prt')
ax.set_ylabel('Valor cebeci')
ax.set_zlabel('L2 norm')
ax.view_init(30, 100)
plt.show()