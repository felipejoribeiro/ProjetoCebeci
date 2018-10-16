import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import rc
from matplotlib import cm
import os
from os import path

outpath = os.getcwd()    # get current directory


# Dados colhidos nos testes de ajuste do modelo Prandtl turbulento em função do Reynolds turbulento e do número de Prandtl.
##########################################################################################################################################

         #Para N = 100 e Re_t = 395
X = [ -0.04 , -0.035 , -0.03 , -0.025 , -0.020 , -0.015 , -0.01  , -0.008 , -0.007 , -0.006 , -0.005    ,   0    , 0.005 , 0.01 , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04]  # V
Y = [ 1.13  , 0.88  ,  0.65 , 0.45   ,  0.3715,  0.47  , 0.67 ,  0.7763, 0.8263  , 0.8771  , 0.9286   , 1.1944   , 1.46 , 1.75 , 2.03  ,  2.32,  2.61 , 2.91 , 3.20  , 3.51]   # L2
plt.plot( X, Y, 'k', linewidth = 1 , label = "N = 100" )

         #Para N = 400 e Re_t = 395
XX = [ -0.04 , -0.035 , -0.03 , -0.025 ,-0.02 , -0.018 , -0.016 , -0.014 , -0.012 , -0.01  , -0.008 , -0.006 , -0.004 , -0.002 , 0, 0.005 , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04]   # V
YY = [ 1.88  , 1.60   , 1.33  , 1.05    ,0.78  ,  0.67  ,  0.58  ,  0.48  ,  0.41  ,  0.35  ,  0.33   , 0.35  ,   0.40  ,  0.48 , 0.57, 0.83 , 1.41 , 1.71  , 2.01  , 2.32 , 2.63  , 2.94]    # L2
plt.plot( XX, YY, 'k', linewidth = 1, linestyle='-.' , label = "N = 400")

         #Para N = 1000
X1000 = [ -0.04 , -0.035 , -0.03 , -0.025 , -0.020 , -0.015 , -0.01 , -0.005 , 0 , 0.005 , 0.01  , 0.015 , 0.02 , 0.025 , 0.03 , 0.035 , 0.04] # V
Y1000 = [  1.93 ,  1.656 ,  1.37 ,  1.101 ,  0.829 ,  0.575 ,  0.37 , 0.3609 ,0.54 ,0.79 , 1.07  , 1.370 , 1.66 , 1.97  , 2.27 , 2.587 , 2.90] # L2
plt.figure(1)
plt.plot( X1000, Y1000, 'k', linewidth = 1, linestyle=':' , label = 'N = 1000')
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.title("Erros no juste de V, independeência de malha")
plt.legend()
plt.xlabel('Valores para o número de Cebeci')
plt.ylabel('Norma L2')
plt.show()