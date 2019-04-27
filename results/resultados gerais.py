import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import rc
from matplotlib import cm
import os
from os import path
import platform

outpath = os.getcwd()    # get current directory


plt.rc('text', usetex=True)
if platform.system() == "Windows":
    plt.rc('font',**{'family':'DejaVu Sans','serif':['Times']})
else :
    plt.rc('font',**{'family':'serif','serif':['Times']})

# plt.rc('font',**{'family':'serif'})
# plt.rc('text', usetex=True)



plt.rcParams.update({'font.size': 14})

##########################################################################################################################################
# Leitura e Plotagem e gráfico tridimensional com análise geral dos resultados do algorítmo para DNS's gerais para o método ajustado (Prandtl e Cebeci).
##########################################################################################################################################

orto = np.loadtxt("ResultadosGeraisOrtodoxos.txt", dtype='float')
C905 = np.loadtxt("ResultadosPrtFixoIdealVs26.txt", dtype='float')
Cc905 = np.loadtxt("ResultadosPrtvariVs26.txt", dtype='float')
RcCm = np.loadtxt("ResultadosPrtFixoIdealVsMod.txt", dtype='float')
Mode = np.loadtxt("ResultadosGeraisModelados.txt", dtype='float')
Gene = np.loadtxt("ResultadosGeraisGenetic.txt", dtype='float')
Gene2p = np.loadtxt("ResultadosGeraisGenetic2temperature.txt", dtype='float')

fig = plt.figure(figsize=(11, 6))                                                                # Determining object image
ax = fig.add_subplot(111, projection='3d')                                                       # Adding 3D axis

xs = orto[: , 0]
ys = orto[: , 1]
zs = orto[: , 3]
ax.scatter(xs, ys, zs , color = "b" , label = 'Prt = 0.71 and A = 26')

xs = C905[: , 0]
ys = C905[: , 1]
zs = C905[: , 3]
ax.scatter(xs, ys, zs , color = 'orange', label = 'Prt = 0.905 and A = 26')

xs = Mode[: , 0]
ys = Mode[: , 1]
zs = Mode[: , 3]
ax.scatter(xs, ys, zs , color = 'black', label = 'Prt(Rey, Pr) and A(Rey)')

xs = Cc905[: , 0]
ys = Cc905[: , 1]
zs = Cc905[: , 3]
ax.scatter(xs, ys, zs , color = 'green', label = 'Prt(Rey) and A = 26')

xs = RcCm[: , 0]
ys = RcCm[: , 1]
zs = RcCm[: , 3]
ax.scatter(xs, ys, zs , color = 'purple', label = 'Prt = 0.905 and A(Ret)')

xs = Gene[: , 0]
ys = Gene[: , 1]
zs = Gene[: , 3]
ax.scatter(xs, ys, zs , color = 'red', label = 'Multivariavel')

xs = Gene2p[: , 0]
ys = Gene2p[: , 1]
zs = Gene2p[: , 3]
ax.scatter(xs, ys, zs , color = 'pink', label = 'Multivariavel com 2 A')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.title('Analise geral dos metodos')
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



tamanho = 10
aspectratio = 0.4
plt.figure(figsize=(tamanho , tamanho * aspectratio))
i = [orto[1 , 0] , orto[3 , 0] , orto[10 , 0] ,  orto[11 , 0]]
ii = [orto[1 , 3] , orto[3 , 3] , orto[10 , 3] ,  orto[11 , 3]]
plt.plot( i, ii, color = 'black' , linestyle='-', label = r'$Pr_t = 0.71$ and $A = 26$')
i = [C905[1 , 0] , C905[3 , 0] , C905[10 , 0] ,  C905[11 , 0]]
ii = [C905[1 , 3] , C905[3 , 3] , C905[10 , 3] ,  C905[11 , 3] ]
plt.plot( i, ii, color = 'black' , linestyle='-.' , label = r'$Pr_t = 0.905$ and $A = 26$')
i = [Cc905[1 , 0] , Cc905[3 , 0] , Cc905[10 , 0] ,  Cc905[11 , 0]]
ii = [Cc905[1 , 3] , Cc905[3 , 3] , Cc905[10 , 3] ,  Cc905[11 , 3]]
plt.plot( i, ii, color = 'black' , linestyle='--' , label = r'$Pr_t(Re_\tau)$ and $A = 26$')
i = [Mode[1 , 0] , Mode[3 , 0] , Mode[10 , 0] ,  Mode[11 , 0]]
ii = [Mode[1 , 3] , Mode[3 , 3] , Mode[10 , 3] ,  Mode[11 , 3]]
plt.plot( i, ii, color = 'black' , linestyle=(0, (1,1)), label = r'$Pr_t(Re_\tau)$ and $A(Re_\tau)$')
# i = [RcCm[1 , 0] , RcCm[3 , 0] , RcCm[10 , 0] ,  RcCm[11 , 0]]
# ii = [RcCm[1 , 3] , RcCm[3 , 3] , RcCm[10 , 3] ,  RcCm[11 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'purple' , label = 'Prt = 0.905 and A(Ret)')
# i = [Gene[1 , 0] , Gene[3 , 0] , Gene[10 , 0] ,  Gene[11 , 0]]
# ii = [Gene[1 , 3] , Gene[3 , 3] , Gene[10 , 3] ,  Gene[11 , 3]]
# plt.plot( i, ii, color = 'black' , linestyle=(0, (1,10)) , label = r'Multi-objective')
i = [Gene2p[1 , 0] , Gene2p[3 , 0] , Gene2p[10 , 0] ,  Gene2p[11 , 0]]
ii = [Gene2p[1 , 3] , Gene2p[3 , 3] , Gene2p[10 , 3] ,  Gene2p[11 , 3]]
plt.plot( i, ii, 'k', linewidth = '1', color = 'black' , label = r'$Pr_t(Re_\tau)$, $A_t(Re_\tau)$ and $A_v(Re_\tau)$ ')
plt.xlabel(r'\textbf{$  Re_\tau $}')
plt.ylabel(r' L2 norm')
plt.legend(fontsize=15 , frameon=False)
plt.grid(color='black', linestyle=':', linewidth=0.5)


# Para Ret = 395


# plt.figure(3)
# i = [orto[2 , 1] , orto[3 , 1] , orto[4 , 1] ,  orto[5 , 1] , orto[6 , 1] , orto[7 , 1] , orto[8 , 1]]
# ii = [orto[2 , 3] , orto[3 , 3] , orto[4 , 3] ,  orto[5 , 3], orto[6 , 3] , orto[7 , 3] , orto[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'b' , label = 'Prt = 0.71 and A = 26')
# i = [Mode[2 , 1] , Mode[3 , 1] , Mode[4 , 1] ,  Mode[5 , 1] , Mode[6 , 1] , Mode[7 , 1] , Mode[8 , 1]]
# ii = [Mode[2 , 3] , Mode[3 , 3] , Mode[4 , 3] ,  Mode[5 , 3], Mode[6 , 3] , Mode[7 , 3] , Mode[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'black' , label = 'Prt(Rey, Pr) and A(Rey)')
# i = [C905[2 , 1] , C905[3 , 1] , C905[4 , 1] ,  C905[5 , 1] , C905[6 , 1] , C905[7 , 1] , C905[8 , 1]]
# ii = [C905[2 , 3] , C905[3 , 3] , C905[4 , 3] ,  C905[5 , 3], C905[6 , 3] , C905[7 , 3] , C905[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'orange' , label = 'Prt = 0.905 and A = 26')
# i = [Cc905[2 , 1] , Cc905[3 , 1] , Cc905[4 , 1] ,  Cc905[5 , 1] , Cc905[6 , 1] , Cc905[7 , 1] , Cc905[8 , 1]]
# ii = [Cc905[2 , 3] , Cc905[3 , 3] , Cc905[4 , 3] ,  Cc905[5 , 3], Cc905[6 , 3] , Cc905[7 , 3] , Cc905[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'green' , label = 'Prt(Rey) and A = 26')
# i = [RcCm[2 , 1] , RcCm[3 , 1] , RcCm[4 , 1] ,  RcCm[5 , 1] , RcCm[6 , 1] , RcCm[7 , 1] , RcCm[8 , 1]]
# ii = [RcCm[2 , 3] , RcCm[3 , 3] , RcCm[4 , 3] ,  RcCm[5 , 3], RcCm[6 , 3] , RcCm[7 , 3] , RcCm[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'purple' , label = 'Prt = 0.905 and A(Ret)')
# i = [Gene[2 , 1] , Gene[3 , 1] , Gene[4 , 1] ,  Gene[5 , 1] , Gene[6 , 1] , Gene[7 , 1] , Gene[8 , 1]]
# ii = [Gene[2 , 3] , Gene[3 , 3] , Gene[4 , 3] ,  Gene[5 , 3], Gene[6 , 3] , Gene[7 , 3] , Gene[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'red' , label = 'Multivariável')
# i = [Genepr[2 , 1] , Genepr[3 , 1] , Genepr[4 , 1] ,  Genepr[5 , 1] , Genepr[6 , 1] , Genepr[7 , 1] , Genepr[8 , 1]]
# ii = [Genepr[2 , 3] , Genepr[3 , 3] , Genepr[4 , 3] ,  Genepr[5 , 3], Genepr[6 , 3] , Genepr[7 , 3] , Genepr[8 , 3]]
# plt.plot( i, ii, 'k', linewidth = '1', color = 'pink' , label = 'Multivariável com Pr')
# plt.legend()
# plt.xlabel('Prandtl')
# plt.ylabel('Norma L2')
# plt.grid(color='black', linestyle=':', linewidth=0.5)

# if platform.system() == "Windows":
#     plt.savefig('graficos\images\gerais.pdf' , bbox_inches='tight')
# else:
#     plt.savefig('graficos/images/gerais.pdf' , bbox_inches='tight')


# plt.subplots_adjust(top=0.95 , bottom=0.19)
plt.show()
##########################################################################################################################################
