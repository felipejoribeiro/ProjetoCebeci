import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import rc
from matplotlib import cm
import os
from os import path


outpath = os.getcwd()    # get current directory


######################################################################################### EXEMPLO FORMATACAO




tamanho = 6
aspectratio = 1/1

dadosin = np.loadtxt("saida3 (1).txt", dtype='float')
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=20)

plt.figure(figsize=(tamanho , tamanho * aspectratio))
xs = dadosin[: , 0]
ys = dadosin[: , 2]
plt.ylim((0, 0.003))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=20)
plt.xlabel('xlabel', fontsize=15 )
plt.ylabel('ylabel', fontsize=20 )
plt.xlabel('CFL')
plt.ylabel('Norma L2')
plt.plot(xs, ys, 'k', linewidth = '2' , label = "L2")
plt.show()



plt.figure(figsize=(tamanho , tamanho * aspectratio))
xs = dadosin[: , 0]
ys = dadosin[: , 3]
plt.plot(xs, ys, 'k', linewidth = '2' ,label = "L1" )
plt.ylim((0, 0.003))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.rc('xtick', labelsize=5)
plt.rc('ytick', labelsize=5)
plt.xlabel('xlabel', fontsize=5 )
plt.ylabel('ylabel', fontsize=5 )
plt.xlabel('CFL')
plt.ylabel('Norma L1')
plt.show()



plt.figure(figsize=(tamanho , tamanho * aspectratio))
xs = dadosin[: , 0]
ys = dadosin[: , 3]
plt.plot(xs, ys, 'k', linewidth = '2' ,label = "L1" )
plt.ylim((0, 0.003))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.rc('xtick', labelsize=5)
plt.rc('ytick', labelsize=5)
plt.xlabel('xlabel', fontsize=5 )
plt.ylabel('ylabel', fontsize=5 )
plt.xlabel('CFL')
plt.ylabel('Norma L1')
plt.show()




#########################################################################################################################################




############################################# EXEMPLO FUNCIONAL INTELIGIVEL SURF:


X = np.zeros((21,21), dtype=np.float64)
Y = np.zeros((21,21), dtype=np.float64)
Z = np.zeros((21,21), dtype=np.float64)
for i in range(21):
    X[i , :] = [-10 , -9 , -8 , -7 , -6 , -5 , -4 , -3 , -2 , -1 , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10]
    Y[: , i] = [-10 , -9 , -8 , -7 , -6 , -5 , -4 , -3 , -2 , -1 , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10]
for i in range(21):
   for ii in range(21):
       Z[i , ii] = X[i , ii] * X[i , ii] + Y[i , ii] * Y[i , ii]

# Normalizando para [0 : 1]
Z = (Z-Z.min())/(Z.max()-Z.min())
# Definindo cores
colors = cm.viridis(Z)
rcount, ccount, _ = colors.shape
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount, facecolors=colors, shade=False)
surf.set_facecolor((0,0,0,0))
plt.show()


############################################ EXEMPLO FUNCIONAL INTELIGIVEL SCATTER:

# Fixing random state for reproducibility
np.random.seed(19680801)


def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
    xs = randrange(n, 23, 32)
    ys = randrange(n, 0, 100)
    zs = randrange(n, zlow, zhigh)
    ax.scatter(xs, ys, zs, c=c, marker=m)

plt.show()