import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import h5py
import os
from os import path
import platform

N = 400
Ret = 150
path = os.getcwd().replace('results/graficos/new' , '')

# Space vector
e = np.zeros(N)
for i in range(1 , N):
    e[i] = (i - 0.5) * (1/(N - 0.5) * int(Ret))

with h5py.File(path + "/MFSim/outPut.hdf5", "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    print("Chombo_global: %s" % f['Chombo_global'])
    print("box_id: %s" % f['box_id'])
    print("floatvect_id: %s" % f['floatvect_id'])
    print("intvect_id: %s" % f['intvect_id'])
    print("level_0: %s" % f['level_0'])
    print("level_0/boxes: %s" % f['level_0']['boxes'])
    print("level_0/data:datatype=0: %s" % f['level_0']['data:datatype=0'])
    print("level_0/data:dataoffsets=0: %s" % f['level_0']['data:offsets=0'])
    print("level_0/data:dataattributes: %s" % f['level_0']['data_attributes'])

# Image config
tamanho = 10
aspectratio = 0.4
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams.update({'font.size': 32})
plt.figure(figsize=(tamanho , tamanho * aspectratio))

# Plot temperature
# dados = np.loadtxt(path + "/results/graficos/new/image 150._0.71_400_classico.txt", dtype='float')
# dns = np.loadtxt(path + "/DNS/DNS_RE_150_071.txt", dtype='float')
#
# plt.plot(e , dados , color='black' , linestyle='-', label= r'Present Work')
# plt.plot(Ret - dns[:, 1] , dns[:, 2] , color='black' , linestyle='-.', label= r'DNS')
#
# plt.xlim(0 , Ret)
# plt.ylim(0 , max([max(dados) , max(dns[:, 2])]) * 1.2)
# plt.legend(fontsize=35 , frameon=False)
# ax = plt.gca()
# label = ax.set_xlabel(r' \textbf{ $ \tilde{y} $ }', fontsize=35)
# ax.xaxis.set_label_coords(0.84, -0.025)
# label = ax.set_ylabel(r'\textbf{$ \tilde{\overline{T^\ast}} $}',fontsize=38)
# ax.yaxis.set_label_coords(-0.03, 0.70)
# plt.show()

# Plot velocity
# dados = np.loadtxt(path + "/results/graficos/new/image 150._400_classico.txt", dtype='float')
# dns = np.loadtxt(path + "/DNS/DNS_RE_150.txt", dtype='float')
#
# plt.plot(e , dados , color='black' , linestyle='-', label= r'Present Work')
# plt.plot(Ret - dns[:, 1] , dns[:, 2] , color='black' , linestyle='-.', label= r'DNS')
#
# plt.xlim(0 , Ret)
# plt.ylim(0 , max([max(dados) , max(dns[:, 2])]) * 1.2)
# plt.legend(fontsize=35 , frameon=False)
# ax = plt.gca()
# label = ax.set_xlabel(r' \textbf{ $ \tilde{y} $ }', fontsize=35)
# ax.xaxis.set_label_coords(0.84, -0.025)
# label = ax.set_ylabel(r'\textbf{$ \tilde{\overline{u}} $}',fontsize=38)
# ax.yaxis.set_label_coords(-0.03, 0.70)
# plt.show()
