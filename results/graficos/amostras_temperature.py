import numpy as np
import matplotlib.pyplot as plt
import math
import os
from os import path

path = os.getcwd()
path = path.replace('results\graficos' , '')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#########################
##    controles        ##
#########################


plt.rcParams.update({'font.size': 22})


Ret = " 150"
dnss = "150_0025"
Pr = 0.03
e = np.zeros(400)
metodo = "classico"

tamanho = 10
aspectratio = 0.5


######################### vetor espaço...

for i in range(1 , 400):
    e[i] = (i - 0.5) * (1/(400 - 0.5) * int(Ret))





######################## Gráfico do Prandtl...




prandtl = np.loadtxt(path + "DNS\Prt_RE_640_071.txt", dtype='float')
spaco = prandtl[: , 1]
prt = prandtl[: , 5]

plt.plot(spaco / 640 ,prt , color='black' , linestyle=":" , label=r"$Re_\tau = 640$ , $Pr = 0.71$")

prandtl = np.loadtxt(path + "DNS\Prt_RE_640_0025.txt", dtype='float' )
spaco = prandtl[: , 1]
prt = prandtl[: , 5]

plt.plot(spaco /640,prt , color='black', linestyle="--" , label=r"$Re_\tau = 640$ , $Pr = 0.025$")

prandtl = np.loadtxt(path + "DNS\Prt_RE_395_071.txt", dtype='float')
spaco = prandtl[: , 1]
prt = prandtl[: , 5]

plt.plot(spaco /395,prt, color='black', linestyle="-" , label=r"$Re_\tau = 395$ , $Pr = 0.71$")

prandtl = np.loadtxt(path + "DNS\Prt_RE_395_0025.txt", dtype='float' )
spaco = prandtl[: , 1]
prt = prandtl[: , 5]

plt.plot(spaco /395,prt, color='black', linestyle="-.", label=r"$Re_\tau = 395$ , $Pr = 0.025$")



plt.figure(figsize=(tamanho , tamanho * aspectratio))
plt.xlabel(r'  $ \tilde{y} $ ',fontsize=28)
plt.ylabel(r'$ Pr_t $',fontsize=28)
plt.subplots_adjust(top=0.9 , bottom=0.17)
plt.xlim(0,1)
plt.legend(fontsize=17)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0) , fontsize= 20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0) , fontsize= 30)
plt.savefig('images/DNS_PRt.png')
plt.show(block=False)






####################### Gráficos gerais das temperaturas...

# dados = np.loadtxt("image" + Ret + "._" + str(Pr)  + "_400_" + metodo + ".txt", dtype='float')
# dns = np.loadtxt(path + "DNS/DNS_RE_"+dnss + ".txt", dtype='float')




# plt.figure(figsize=(tamanho , tamanho * aspectratio))

# plt.plot(e , dados)
# plt.plot(- e , dados)

# plt.plot(int(Ret) - dns[:, 1] , dns[:, 2])
# plt.plot(dns[:, 1]- int(Ret) , dns[:, 2])

# plt.xlim(- int(Ret) , int(Ret))
# plt.ylim(0 , max([max(dados) , max(dns[:, 2])]) * 1.2)
# plt.show()




