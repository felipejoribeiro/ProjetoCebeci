import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



scater = np.loadtxt("Gene395_071_400.txt", dtype='float')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(scater[: , 0], scater[: , 1], scater[: , 2])
ax.set_xlabel('Prt')
ax.set_ylabel('Valor cebeci')
ax.set_zlabel('L2 norm')
plt.show()