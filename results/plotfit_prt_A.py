import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import rc
from matplotlib import cm
import os
from os import path


outpath = os.getcwd()    # get current directory





# # Fit dos valores numéricos do número de Cebeci:


#ajuste Prandtl sem cebeci modelado:
xPrandtl = np.array([ (150.0) , (395.0) , (640.0)  , (1020.0)])
yPrandtl = np.array([(0.94531249999) , (0.8953124999999999) , (0.89531249999999996) , (0.89999999999999991)])
zPrandtl = np.polyfit(xPrandtl , yPrandtl , 3)
print("Modelo para o Prt, cebeci n modelado=" , zPrandtl[0] , "* Ret**3 + " , zPrandtl[1], "* Ret**2 + " , zPrandtl[2] ,  "* Ret + " , zPrandtl[3])
print(zPrandtl)
print("----------------------------------------------")
plt.figure(3)
plt.plot( xPrandtl , yPrandtl ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
xPrandtlvontinum = np.arange( (150) , (1020) , 0.1)
yPrandtlmodelado =  zPrandtl[0] * np.multiply(xPrandtlvontinum, np.multiply(xPrandtlvontinum, xPrandtlvontinum))+ zPrandtl[1] * np.multiply(xPrandtlvontinum, xPrandtlvontinum) + zPrandtl[2] * xPrandtlvontinum + zPrandtl[3]
plt.plot( xPrandtlvontinum , yPrandtlmodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
plt.plot( xPrandtl , (np.zeros(4) + 0.906) ,  'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905')
plt.legend()
plt.title('Ajuste do valor de Prandtl turbulento, cebeci não modelado')
plt.xlabel('Ret')
plt.ylabel('Prt')
plt.show()


#ajuste cebeci:
xcebeci = np.array([ math.log(150.0) , math.log(395.0) , math.log(640.0)  , math.log(1020.0)])
ycebeci = np.array([math.log(28.616180419921875) , math.log(25.673782348632812) , math.log(25.001266479492188) , math.log(25.002136230468750)])
zcebeci = np.polyfit(xcebeci , ycebeci , 2)
plt.figure(4)
plt.plot( (math.exp(1))**(xcebeci) , (math.exp(1)) ** (ycebeci) ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
xcebecivontinum = np.arange( math.log(150) , math.log(1020) , 0.1)
ycebecimodelado =   (math.exp(1))**(0.04510771 * np.multiply(xcebecivontinum, xcebecivontinum) - 0.60942948 * xcebecivontinum + 5.27533404)
plt.plot( (math.exp(1))**xcebecivontinum , ycebecimodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
ycebecimodelado =   (math.exp(1))**( zcebeci[0] * np.multiply(xcebecivontinum, xcebecivontinum) + zcebeci[1] * xcebecivontinum + zcebeci[2])
print("Modelo para o cebeci modelado =" , zcebeci[0], "* Ret**2 + " , zcebeci[1] ,  "* Ret + " , zcebeci[2])
print(zcebeci)
print("----------------------------------------------")
plt.plot( (math.exp(1))**xcebecivontinum , ycebecimodelado ,  'k', linewidth = '1', color = 'orange' , label = 'ajustenovo')
plt.plot( (math.exp(1))**(xcebeci) , (np.zeros(4) + 26) ,  'k', linewidth = '0.8', color = 'blue' , label = 'C = 26')
print()
plt.legend()
plt.title('Ajuste do valor de Cebeci')
plt.xlabel('Ret')
plt.ylabel('A')
plt.show()



#ajuste Prandtl com cebeci modelado:
xPrandtl = np.array([ (150.0) , (395.0) , (640.0)  , (1020.0)])
yPrandtl = np.array([(0.885937) , (0.90156) , (0.910937) , (0.91406)])
zPrandtl = np.polyfit(xPrandtl , yPrandtl , 3)
print("Modelo para o Prt, cebeci modelado =" , zPrandtl[0] , "* Ret**3 + " , zPrandtl[1], "* Ret**2 + " , zPrandtl[2] ,  "* Ret + " , zPrandtl[3])
print(zPrandtl)
print("---------------------------------------------- ")
plt.figure(3)
plt.plot( xPrandtl , yPrandtl ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
xPrandtlvontinum = np.arange( (150) , (1020) , 0.1)
yPrandtlmodelado =  zPrandtl[0] * np.multiply(xPrandtlvontinum, np.multiply(xPrandtlvontinum, xPrandtlvontinum))+ zPrandtl[1] * np.multiply(xPrandtlvontinum, xPrandtlvontinum) + zPrandtl[2] * xPrandtlvontinum + zPrandtl[3]
plt.plot( xPrandtlvontinum , yPrandtlmodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
plt.plot( xPrandtl , (np.zeros(4) + 0.906) ,  'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905')
plt.legend()
plt.title('Ajuste do valor de Prandtl turbulento, cebeci modelado')
plt.xlabel('Ret')
plt.ylabel('Prt')
plt.show()




#ajuste com o algorítmo genético diferencial do Prt e Cebeci:
xPrandtl = np.array([ (150.0) , (395.0) , (640.0)  , (1020.0)])
yPrandtl = np.array([(0.60035962262979192) , (0.75938280042931627) , ( 0.83398599512545113 ) , (0.84876244186418570)])
zPrandtl = np.polyfit(xPrandtl , yPrandtl , 3)
print("Modelo para o Prt, genetic=" , zPrandtl[0] , "* Ret**3 + " , zPrandtl[1], "* Ret**2 + " , zPrandtl[2] ,  "* Ret + " , zPrandtl[3])
print(zPrandtl)
print("---------------------------------------------- ")
plt.figure(3)
plt.plot( xPrandtl , yPrandtl ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
xPrandtlvontinum = np.arange( (150) , (1020) , 0.1)
yPrandtlmodelado =  zPrandtl[0] * np.multiply(xPrandtlvontinum, np.multiply(xPrandtlvontinum, xPrandtlvontinum))+ zPrandtl[1] * np.multiply(xPrandtlvontinum, xPrandtlvontinum) + zPrandtl[2] * xPrandtlvontinum + zPrandtl[3]
plt.plot( xPrandtlvontinum , yPrandtlmodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
plt.plot( xPrandtl , (np.zeros(4) + 0.906) ,  'k', linewidth = '1', color = 'blue' , label = 'Prt = 0.905')
plt.legend()
plt.title('Ajuste do valor de Prandtl turbulento, genetic')
plt.xlabel('Ret')
plt.ylabel('Prt')
plt.show()

#cebeci:
xcebeci = np.array([ math.log(150.0) , math.log(395.0) , math.log(640.0)  , math.log(1020.0)])
ycebeci = np.array([math.log(46.730289505179861) , math.log(34.938766563537811) , math.log(30.114934486932917) , math.log(29.921388094701264)])
zcebeci = np.polyfit(xcebeci , ycebeci , 3)
plt.figure(4)
plt.plot( (math.exp(1))**(xcebeci) , (math.exp(1)) ** (ycebeci) ,  'k', linewidth = '1', color = 'black' , label = 'Dados numéricos')
xcebecivontinum = np.arange( math.log(150) , math.log(1020) , 0.1)
ycebecimodelado =   math.exp(1)**(zcebeci[0] * np.multiply(xcebecivontinum, np.multiply(xcebecivontinum, xcebecivontinum))+ zcebeci[1] * np.multiply(xcebecivontinum, xcebecivontinum) + zcebeci[2] * xcebecivontinum + zcebeci[3])
plt.plot( (math.exp(1))**xcebecivontinum , ycebecimodelado ,  'k', linewidth = '1', color = 'red' , label = 'ajuste')
print("Modelo para o cebeci, genetic=  exp(" , zcebeci[0] , "* Ret**3 + " , zcebeci[1], "* Ret**2 + " , zcebeci[2] ,  "* Ret + " , zcebeci[3] , ")")
print(zcebeci)
print("----------------------------------------------")
plt.plot( (math.exp(1))**(xcebeci) , (np.zeros(4) + 26) ,  'k', linewidth = '0.8', color = 'blue' , label = 'C = 26')
print()
plt.legend()
plt.title('Ajuste do valor de Cebeci, genetic')
plt.xlabel('Ret')
plt.ylabel('A')
plt.show()
##########################################################################################################################################
