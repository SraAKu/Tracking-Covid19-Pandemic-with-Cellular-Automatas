import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

L = 500                 #Dimensiones
NCeldas = L*L           #Numero total celdas

Densidad = 'Constante=1'
Conectividad = 'Constante=1'
accion = 'nula'

miu = 7e-1             #Coeficiente de infección a celdas vecinas  
v = 1e-2               #Virulencia de la epidemica/Coefiiente de infección dentro de la celda
ro = 1/NCeldas        #Densidad de población de la red/Constante por ahora

compendio = 'AC_D:%s_C:%s_A:%s_L=%d,.txt' %(Densidad,Conectividad,accion,L)
file1 = open(compendio, 'r')
Y = np.loadtxt(file1)
nsteps = len(Y[:,0])
NewI = np.zeros(nsteps)
NewI[0] = Y[0,1]
for i in range(nsteps-1):
    NewI[i+1] = Y[i+1,1] - Y[i,1]

fig = plt.figure(figsize=(16,9))
plt.xlim(0,max(Y[:,0]))
plt.ylim(0,NCeldas)
plt.grid()
plt.plot(Y[:,0],Y[:,1],label='infectados',color='red')
plt.plot(Y[:,0],Y[:,2],label='susceptibles',color='green')
plt.plot(Y[:,0],Y[:,3],label='NuevosInfectads',color='blue')
plt.plot(Y[:,0],Y[:,4],label='Recuperados',color='orange')
plt.plot(Y[:,0],Y[:,5],label='Infectados_t',color='purple')
plt.legend()
plt.savefig('AC_D:%s_C:%s_A:%s_L=%d,.png' %(Densidad,Conectividad,accion,L))

fig2 = plt.figure(figsize=(16,9))
plt.xlim(0,max(Y[:,0]))
plt.ylim(0,max(Y[:,3]))
plt.grid()
plt.plot(Y[:,0],Y[:,3],label='NuevosInfectads',color='blue')
plt.legend()
plt.savefig('NuevosInfectados:%s_C:%s_A:%s_L=%d,.png' %(Densidad,Conectividad,accion,L))

fig3 = plt.figure(figsize=(16,9))
plt.xlim(0,max(Y[:,1]))
plt.ylim(0,max(Y[:,3]))
plt.grid()
plt.plot(Y[:,1],Y[:,3],label='NuevosInfectads',color='blue')
plt.legend()
plt.savefig('NuevosInfectados_Infectados:%s_C:%s_A:%s_L=%d,.png' %(Densidad,Conectividad,accion,L))

######################## GRAFICADOR DE SIMULACIÓN VS DATOS REALES ###########################

file2 = pd.read_csv('/home/juan_aku/Escritorio/AC-Covid/Covid_Data/time_series_covid19_confirmed_global.csv')
data = file2.to_numpy()
dataPais = data[137][34:]         #PRIMEROS 100 DATOS
pasos_diaT = 10

#Densidad = 'Normal_Sigma=L|2,L|4.47'
Densidad = 'Normal_Sigma=L|4,L|13.63'
Conectividad = 'Constante=0.15'
accion = 'lineal'
ubicacion = 'D_out Loc=(3L|4 - 70,L|4 + 150)'
compendio = 'AC_D:%s_C:%s_A:%s_U:%s_L=%d,.txt' %(Densidad,Conectividad,accion,ubicacion,L)
filed = open(compendio, 'r')
Y1 = np.loadtxt(filed)
Is1 = Y1[:,1]
X1 = np.arange(len(Is1))

'''
Densidad = 'Normal_Sigma=L|4,L|13.63'
Conectividad = 'Constante=0.15'
accion = 'nula'
ubicacion = 'D|2'
compendio = 'AC_D:%s_C:%s_A:%s_U:%s_L=%d,.txt' %(Densidad,Conectividad,accion,ubicacion,L)
filed = open(compendio, 'r')
Y2 = np.loadtxt(filed)
Is2 = Y2[:,1]
X2 = np.arange(len(Is2))


Densidad = 'Normal_Sigma=L|4,L|13.63'
Conectividad = 'Constante=0.15'
accion = 'nula'
ubicacion = 'D|2 - 150'
compendio = 'AC_D:%s_C:%s_A:%s_U:%s_L=%d,.txt' %(Densidad,Conectividad,accion,ubicacion,L)
filed = open(compendio, 'r')
Y3 = np.loadtxt(filed)
Is3 = Y3[:,1]
X3 = np.arange(len(Is3))


Densidad = 'Normal_Sigma=L|4,L|13.63'
Conectividad = 'Constante=0.15'
accion = 'nula'
ubicacion = 'D|2 + 150'
compendio = 'AC_D:%s_C:%s_A:%s_U:%s_L=%d,.txt' %(Densidad,Conectividad,accion,ubicacion,L)
filed = open(compendio, 'r')
Y4 = np.loadtxt(filed)
Is4 = Y4[:,1]
X4 = np.arange(len(Is4))


Densidad = 'Normal_Sigma=L|4,L|13.63'
Conectividad = 'Constante=0.15'
accion = 'nula'
ubicacion = 'D_out Loc=(3L|4 - 70,L|4 + 150)'
compendio = 'AC_D:%s_C:%s_A:%s_U:%s_L=%d,.txt' %(Densidad,Conectividad,accion,ubicacion,L)
filed = open(compendio, 'r')
Y5 = np.loadtxt(filed)
Is5 = Y5[:,1]
X5 = np.arange(len(Is5))


#Densidad = 'Normal_Sigma=L|5'
#Conectividad = 'Constante=0.25'
#accion = 'nula'
#compendio = 'AC_D:%s_C:%s_A:%s_L=%d,.txt' %(Densidad,Conectividad,accion,L)
#filed = open(compendio, 'r')
#Y6 = np.loadtxt(filed)
#Is6 = Y6[:,1]
#X6 = np.arange(len(Is6))

'''

print(len(dataPais))
pasos_tiempo = len(Is1)/len(dataPais)    #Correspondencia entera entre pasos de simulación y tiempo
print(pasos_tiempo)
#pasos_tiempo = len(Is2)/len(dataPais)
#print(pasos_tiempo)
#pasos_tiempo = len(Is3)/len(dataPais)
#print(pasos_tiempo)
#pasos_tiempo = len(Is4)/len(dataPais)
#print(pasos_tiempo)
#pasos_tiempo = len(Is5)/len(dataPais)
#print(pasos_tiempo)
#pasos_tiempo = len(Is6)/len(dataPais)
#print(pasos_tiempo)
X0 = np.arange(len(dataPais))

conect0 = 0.15
conect1 = -0.05
dia0 = 0
dias_accion = 15
deltac = (conect0 - conect1)/(dia0 - (dias_accion*pasos_diaT))
conects = list()
conects.append(conect0)

for i in range(dias_accion*pasos_diaT):
    conects.append(conect0 + (deltac*(i+1)))

conects = np.array(conects)
XX = np.arange(len(conects))
XX = XX/pasos_diaT + 36
XX0 = np.arange(0,36)
XX1 = np.arange(36+dias_accion, len(dataPais))
conects0 = np.full(36,conect0)
conects1 = np.full(len(dataPais)-36-dias_accion,conect1)

fig = plt.figure(figsize=(14,9))

plt.subplot(211)
plt.plot(XX0,conects0,color='blue')
plt.plot(XX,conects,color='blue')
plt.plot(XX1,conects1,color='blue')
plt.title('Connectivity over time')
plt.xlim(0,len(dataPais))
plt.ylabel('Connectivity $S$')
plt.legend()
plt.grid(True)

plt.subplot(212)
plt.scatter(X0,dataPais,color='red',label='Confirmed Cases')
plt.scatter(X1/pasos_diaT,Is1,s=6.5,label='Normal Densities\nOutbreak Location = $(L/4 + 150, 3L/4 - 70)$')
#plt.scatter(X2/pasos_diaT,Is2,s=8,label='Normal Densities\nOutbreak Location = $D/2$')
#plt.scatter(X3/pasos_diaT,Is3,s=8,label='Normal Densities\nOutbreak Location = $D/2 - 150$')
#plt.scatter(X4/pasos_diaT,Is4,s=8,label='Normal Densities\nOutbreak Location = $D/2 + 150$')
#plt.scatter(X5/pasos_diaT,Is5,s=8,label='Normal Densities\nOutbreak Location = $(L/4 + 150, 3L/4 - 50)$')
#plt.scatter(X6/pasos_diaT,Is1,s=8,color='green',label='Pop Normaly distributed\n$\sigma=L|5$')
#plt.scatter(X6/pasos_diaT,Is6,s=2,label='Pop Normaly distributed\n$\sigma=L/5$')
#plt.scatter(X5/pasos_diaT,Is5,s=2,label='Pop Normaly distributed\n$\sigma=L/6$')
#plt.scatter(X2/pasos_diaT,Is2,s=2,label='Pop Normaly distributed\n$\sigma=L/9$')
#plt.scatter(X3/pasos_diaT,Is3,s=2,label='Pop Normaly distributed\n$\sigma=L/8$')
#plt.scatter(X4/pasos_diaT,Is4,s=2,label='Pop Normaly distributed\n$\sigma=L/7$')
plt.title('Covid-19 simulated cases for Italy')
plt.xlim(0,len(X0))
plt.ylabel('Cumulative confirmed cases')
#plt.yscale('log')
plt.xlabel('Days from outbreak report')
plt.legend()
plt.grid(True)

plt.savefig('SimuladosvsReportadosLineal(Normal_Sigma=L|2,L|4.47)_accion_%s_conectividad_%s_ubicacion_%s.png' %(accion, Conectividad, ubicacion))
