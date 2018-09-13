import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import os
from os import path


outpath = os.getcwd()    # get current directory


# Dados colhidos nos testes de ajuste do modelo Prandtl turbulento em função do Reynolds turbulento e do número de Prandtl.
##########################################################################################################################################

         #Para N = 100 e Re_t = 395
X = [ -0.04 , -0.035 , -0.03 , -0.025 , -0.020 , -0.015 , -0.01  , -0.008 , -0.007 , -0.006 , -0.005    ,   0    , 0.005 , 0.01 , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04]  # V
Y = [ 1.13  , 0.88  ,  0.65 , 0.45   ,  0.3715,  0.47  , 0.67 ,  0.7763, 0.8263  , 0.8771  , 0.9286   , 1.1944   , 1.46 , 1.75 , 2.03  ,  2.32,  2.61 , 2.91 , 3.20  , 3.51]   # L2
plt.plot( X, Y, 'k', linewidth = '1' , color = "b" , label = "N = 100" )

         #Para N = 400 e Re_t = 395
XX = [ -0.04 , -0.035 , -0.03 , -0.025 ,-0.02 , -0.018 , -0.016 , -0.014 , -0.012 , -0.01  , -0.008 , -0.006 , -0.004 , -0.002 , 0, 0.005 , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04]   # V
YY = [ 1.88  , 1.60   , 1.33  , 1.05    ,0.78  ,  0.67  ,  0.58  ,  0.48  ,  0.41  ,  0.35  ,  0.33   , 0.35  ,   0.40  ,  0.48 , 0.57, 0.83 , 1.41 , 1.71  , 2.01  , 2.32 , 2.63  , 2.94]    # L2
plt.plot( XX, YY, 'k', linewidth = '1', color = "green" , label = "N = 400")

         #Para N = 1000
X1000 = [ -0.04 , -0.035 , -0.03 , -0.025 , -0.020 , -0.015 , -0.01 , -0.005 , 0 , 0.005 , 0.01  , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04] # V
Y1000 = [  1.93 ,  1.656 ,  1.37 ,  1.101 ,  0.829 ,  0.575 ,  0.37 , 0.3609 ,0.54 ,0.79 , 1.07  , 1.370 , 1.66 , 1.97  , 2.27 , 2.587 , 2.90] # L2
plt.figure(1)
plt.plot( X1000, Y1000, 'k', linewidth = '1', color = 'red' , label = 'N = 1000')
plt.title("Erros no juste de V, independeência de malha")
plt.legend()
plt.xlabel('Valores para o número de Cebeci')
plt.ylabel('Norma L2')
plt.show()
##########################################################################################################################################

# Leitura e Plotagem e gráfico tridimensional com análise geral dos resultados do algorítmo para DNS's gerais para o método ajustado (Prandtl e Cebeci).
##########################################################################################################################################

orto = np.loadtxt("results/ResultadosGeraisOrtodoxos.txt", dtype='float')
C905 = np.loadtxt("results/ResultadosPrtFixoIdealVs26.txt", dtype='float')
RcCm = np.loadtxt("results/ResultadosPrtFixoIdealVsMod.txt", dtype='float')
Mode = np.loadtxt("results/ResultadosGeraisModelados.txt", dtype='float')

fig = plt.figure(figsize=(11, 6))                                                                # Determining object image
ax = fig.add_subplot(111, projection='3d')                                        # Adding 3D axis

xs = orto[: , 0]
ys = orto[: , 1]
zs = orto[: , 3]
ax.scatter(xs, ys, zs , color = "b" , label = 'Literature')

xs = Mode[: , 0]
ys = Mode[: , 1]
zs = Mode[: , 3]
ax.scatter(xs, ys, zs , color = 'black', label = 'Prt(Rey, Pr) and A(Rey)')

xs = C905[: , 0]
ys = C905[: , 1]
zs = C905[: , 3]
ax.scatter(xs, ys, zs , color = 'orange', label = 'Prt = 0.905 and A = 26')

xs = RcCm[: , 0]
ys = RcCm[: , 1]
zs = RcCm[: , 3]
ax.scatter(xs, ys, zs , color = 'purple', label = 'Prt = 0.905 and A(Ret)')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.title('Análise geral dos métodos')
ax.set_xlabel('Reynolds turbulento')
ax.set_ylabel('Prandtl')
ax.set_zlabel('norma L2')
ax.legend()
# cria o video
# ii = 0
# for i in range (0, 360 , 2):
#     ax.view_init(azim=i)
#     plt.draw()
#     print (ax.azim)
#     ii = ii + 1
#     m = outpath + "/plot/plot_{0:03d}.png".format(ii)
#     print(m)
#     plt.savefig(m)
# Cortes do gráfico tridimensional acima.
# Para Pr = 0.71
plt.figure(2)
i = [orto[1 , 0] , orto[3 , 0] , orto[10 , 0] ,  orto[11 , 0]]
ii = [orto[1 , 3] , orto[3 , 3] , orto[10 , 3] ,  orto[11 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'red' , label = 'ortodoxo')
i = [Mode[1 , 0] , Mode[3 , 0] , Mode[10 , 0] ,  Mode[11 , 0]]
ii = [Mode[1 , 3] , Mode[3 , 3] , Mode[10 , 3] ,  Mode[11 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'black' , label = 'Prt(Rey, Pr) and A(Rey)')
i = [C905[1 , 0] , C905[3 , 0] , C905[10 , 0] ,  C905[11 , 0]]
ii = [C905[1 , 3] , C905[3 , 3] , C905[10 , 3] ,  C905[11 , 3] ]
plt.plot( i, ii, 'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905 and A = 26')
i = [RcCm[1 , 0] , RcCm[3 , 0] , RcCm[10 , 0] ,  RcCm[11 , 0]]
ii = [RcCm[1 , 3] , RcCm[3 , 3] , RcCm[10 , 3] ,  RcCm[11 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'orange' , label = 'Prt = 0.905 and A(Ret)')
plt.legend()
plt.xlabel('Reynolds turbulento')
plt.ylabel('Norma L2')
# Para Ret = 395
plt.figure(3)
i = [orto[2 , 1] , orto[3 , 1] , orto[4 , 1] ,  orto[5 , 1] , orto[6 , 1] , orto[7 , 1] , orto[8 , 1]]
ii = [orto[2 , 3] , orto[3 , 3] , orto[4 , 3] ,  orto[5 , 3], orto[6 , 3] , orto[7 , 3] , orto[8 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'red' , label = 'ortodoxo')
i = [Mode[2 , 1] , Mode[3 , 1] , Mode[4 , 1] ,  Mode[5 , 1] , Mode[6 , 1] , Mode[7 , 1] , Mode[8 , 1]]
ii = [Mode[2 , 3] , Mode[3 , 3] , Mode[4 , 3] ,  Mode[5 , 3], Mode[6 , 3] , Mode[7 , 3] , Mode[8 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'black' , label = 'Prt(Rey, Pr) and A(Rey)')
i = [C905[2 , 1] , C905[3 , 1] , C905[4 , 1] ,  C905[5 , 1] , C905[6 , 1] , C905[7 , 1] , C905[8 , 1]]
ii = [C905[2 , 3] , C905[3 , 3] , C905[4 , 3] ,  C905[5 , 3], C905[6 , 3] , C905[7 , 3] , C905[8 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905 and A = 26')
i = [RcCm[2 , 1] , RcCm[3 , 1] , RcCm[4 , 1] ,  RcCm[5 , 1] , RcCm[6 , 1] , RcCm[7 , 1] , RcCm[8 , 1]]
ii = [RcCm[2 , 3] , RcCm[3 , 3] , RcCm[4 , 3] ,  RcCm[5 , 3], RcCm[6 , 3] , RcCm[7 , 3] , RcCm[8 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'orange' , label = 'Prt = 0.905 and A(Ret)')
plt.legend()
plt.xlabel('Prandtl')
plt.ylabel('Norma L2')
plt.show()
##########################################################################################################################################


# # Fit dos valores numéricos do número de Cebeci:

# #ajuste Prandtl:
# xPrandtl = np.array([ (150.0) , (395.0) , (640.0)  , (1020.0)])
# yPrandtl = np.array([(0.884375) , (0.90000000) , (0.90937499) , (0.912499)])
# zPrandtl = np.polyfit(xPrandtl , yPrandtl , 3)
# #print(zPrandtl)
# plt.figure(3)
# plt.plot( xPrandtl , yPrandtl ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
# xPrandtlvontinum = np.arange( (150) , (1020) , 0.1)
# yPrandtlmodelado =  4.58738633e-12 * np.multiply(xPrandtlvontinum, np.multiply(xPrandtlvontinum, xPrandtlvontinum))-5.74977771e-08 * np.multiply(xPrandtlvontinum, xPrandtlvontinum) + 9.40210329e-05 * xPrandtlvontinum + 8.71550063e-01
# plt.plot( xPrandtlvontinum , yPrandtlmodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
# plt.plot( xPrandtl , (np.zeros(4) + 0.906) ,  'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905')
# print()
# plt.legend()
# plt.title('Ajuste do valor de Prandtl turbulento')
# plt.xlabel('Ret')
# plt.ylabel('Prt')


# #ajuste cebeci:
# xcebeci = np.array([ math.log(150.0) , math.log(395.0) , math.log(640.0)  , math.log(1020.0)])
# ycebeci = np.array([math.log(28.61621093) , math.log(25.67382812) , math.log(25.00122070) , math.log(25.002197265)])
# zcebeci = np.polyfit(xcebeci , ycebeci , 2)
# #print(zcebeci)
# plt.figure(4)
# plt.plot( (math.exp(1))**(xcebeci) , (math.exp(1)) ** (ycebeci) ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
# xcebecivontinum = np.arange( math.log(150) , math.log(1020) , 0.1)
# ycebecimodelado =   (math.exp(1))**(0.04510771 * np.multiply(xcebecivontinum, xcebecivontinum) - 0.60942948 * xcebecivontinum + 5.27533404)
# plt.plot( (math.exp(1))**xcebecivontinum , ycebecimodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
# plt.plot( (math.exp(1))**(xcebeci) , (np.zeros(4) + 26) ,  'k', linewidth = '0.8', color = 'blue' , label = 'C = 26')
# print()
# plt.legend()
# plt.title('Ajuste do valor de Cebeci')
# plt.xlabel('Ret')
# plt.ylabel('A')
# plt.show()
##########################################################################################################################################




















############################################## EXEMPLO FUNCIONAL INTELIGIVEL SURF:


# X = np.zeros((21,21), dtype=np.float64)
# Y = np.zeros((21,21), dtype=np.float64)
# Z = np.zeros((21,21), dtype=np.float64)
# for i in range(21):
#     X[i , :] = [-10 , -9 , -8 , -7 , -6 , -5 , -4 , -3 , -2 , -1 , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10]
#     Y[: , i] = [-10 , -9 , -8 , -7 , -6 , -5 , -4 , -3 , -2 , -1 , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10]
# for i in range(21):
#    for ii in range(21):
#        Z[i , ii] = X[i , ii] * X[i , ii] + Y[i , ii] * Y[i , ii]

# # Normalizando para [0 : 1]
# Z = (Z-Z.min())/(Z.max()-Z.min())
# # Definindo cores
# colors = cm.viridis(Z)
# rcount, ccount, _ = colors.shape
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount, facecolors=colors, shade=False)
# surf.set_facecolor((0,0,0,0))
# plt.show()


############################################# EXEMPLO FUNCIONAL INTELIGIVEL SCATTER:

# # Fixing random state for reproducibility
# np.random.seed(19680801)


# def randrange(n, vmin, vmax):
#     '''
#     Helper function to make an array of random numbers having shape (n, )
#     with each number distributed Uniform(vmin, vmax).
#     '''
#     return (vmax - vmin)*np.random.rand(n) + vmin

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# n = 100

# # For each set of style and range settings, plot n random points in the box
# # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
# for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#     xs = randrange(n, 23, 32)
#     ys = randrange(n, 0, 100)
#     zs = randrange(n, zlow, zhigh)
#     ax.scatter(xs, ys, zs, c=c, marker=m)

# plt.show()