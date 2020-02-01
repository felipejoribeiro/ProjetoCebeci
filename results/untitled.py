import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import os
from os import path
import platform

plt.rc('text', usetex=True)
if platform.system() == "Windows":
    plt.rc('font',**{'family':'DejaVu Sans','serif':['Times']})
else :
    plt.rc('font',**{'family':'serif','serif':['Times']})


# plt.rc('font',**{'family':'serif'})
# plt.rc('text', usetex=True)


#########################
##    controles        ##
#########################

plt.rcParams.update({'font.size': 15})
tamanho = 10
aspectratio = 0.4

######################## Dados aquisitados...

p = [10, 15 , 20, 25, 30]
e1 = [179.22 , 279.08 , 373.00 , 480.62 , 578.58]
e2 = [41.70 , 63.96 , 86.05 , 109.84 , 134.99]
e3 = [69.91 , 110.00 , 149.99 , 194.90 , 242.25]
e4 = [38.56 , 61.73 , 83.86 , 110.86 , 132.85]
e5 = [146.26 , 226.88 , 302.96 , 398.18 , 470.78]

######################## Estudo de erros... P/ 18

et = [63.64 , 97.14 , 172.20 , 278.44 , 350.04]
ee = [77.21 , 70.58 , 126.00 , 257.31 , 335.43]

erro = [21 , 27.3 , 26.8 , 7.5 , 4]

######################## Gráficos de cada extensometro e linearização...
plt.figure(figsize=(tamanho , tamanho * aspectratio))

linear_e1 = np.polyfit(p , e1 , 1)
linear_e2 = np.polyfit(p , e2 , 1)
linear_e3 = np.polyfit(p , e3 , 1)
linear_e4 = np.polyfit(p , e4 , 1)
linear_e5 = np.polyfit(p , e5 , 1)

P_continum = np.arange( (0) , (30) , 0.1)

plt.plot(P_continum , linear_e1[0] * P_continum + linear_e1[1] , color='black', linewidth = '1' )
plt.plot(P_continum , linear_e2[0] * P_continum + linear_e2[1] , color='red', linewidth = '1' )
plt.plot(P_continum , linear_e3[0] * P_continum + linear_e3[1] , color='green', linewidth = '1' )
plt.plot(P_continum , linear_e4[0] * P_continum + linear_e4[1] , color='blue', linewidth = '1' )
plt.plot(P_continum , linear_e5[0] * P_continum + linear_e5[1] , color='orange', linewidth = '1' )

plt.scatter(p , e1, color='black' ,marker=">"  , label=r"$\epsilon_1$ , $a_1$ = " + str(linear_e1[0]) + r" , $b_1$ = " + str(linear_e1[1]))
plt.scatter(p , e2, color='red' ,marker="2"  , label=r"$\epsilon_2$ , $a_2$ = " + str(linear_e2[0]) + r" , $b_2$ = " + str(linear_e2[1]))
plt.scatter(p , e3, color='green' ,marker= (5, 0) , label=r"$\epsilon_3$, $a_3$ = " + str(linear_e3[0]) + r" , $b_3$ = " + str(linear_e3[1]))
plt.scatter(p , e4, color='blue' ,marker='+'  , label=r"$\epsilon_4$, $a_4$ = " + str(linear_e4[0]) + r" , $b_4$ = " + str(linear_e4[1]))
plt.scatter(p , e5, color='orange' , label=r"$\epsilon_5$, $a_5$ = " + str(linear_e5[0]) + r" , $b_5$ = " + str(linear_e5[1]))

plt.xlabel(r'  $ \Large \epsilon \normalsize \left[  10^{-6} {m}/{m}  \right]$ ',fontsize=20)
plt.ylabel(r'$\Large P \normalsize \left[{kgf}/{cm^2}\right] $',fontsize=20)
plt.legend(fontsize=15 , frameon=False)
plt.subplots_adjust(top=0.9 , bottom=0.19)
plt.title(r'Regre\c{c}\~{a}o linear dos resultados experimentais', fontsize=20)

plt.show()



######################## Gráficos de erros por deformação ...
fig, ax1 = plt.subplots(figsize=(tamanho , tamanho * aspectratio))
ax2 = ax1.twinx()


ax1.scatter([1 , 2 , 3 , 4 , 5], ee, color='black' ,marker=">", label=r" Deforma\c{c}\~{a}o experimental." )
ax1.scatter([1 , 2 , 3 , 4 , 5], et, color='red' ,marker="2"  , label=r" Deforma\c{c}\~{a}o te\'orica.")
ax2.plot([1 , 2 , 3 , 4 , 5] , erro , color='blue', linewidth = '1' )



ax1.set_xlabel(r'  $  sensor_{n}  $ ',fontsize=20)
ax1.set_ylabel(r'$\Large \epsilon \normalsize  \left[  10^{-6} {m}/{m}  \right] $',fontsize=20)
ax2.set_ylabel(r'Erro $\left[\% \right]$', color = 'blue' , fontsize=20)
ax2.tick_params('y', colors='blue')
plt.rc('xtick', labelsize=35)
plt.title(r'Erro de cada sensor', fontsize=20)
plt.subplots_adjust(top=0.9 , bottom=0.19)
leg = ax1.legend(fontsize=15 , frameon=False)

plt.draw()
bb = leg.get_bbox_to_anchor().inverse_transformed(ax1.transAxes)

bb.x0 += 0
bb.x1 += 0
bb.y0 += -0.3
bb.y1 += -0.3
leg.set_bbox_to_anchor(bb, transform = ax1.transAxes)


plt.show()










