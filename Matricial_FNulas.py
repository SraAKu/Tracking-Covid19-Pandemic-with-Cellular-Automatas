###################################################################################################
###              METODO DE AUTOMATAS CELULARES APLICADO A LA PROPAGACIÓN EPIDEMIOLOGICA         ###
###                                   MULTIPLICACION DE MATRICES                                ###
###    Juan Diego Alzate Giraldo - Universidad Nacioal de Colombia - judalzategi@unal.edu.co    ###
###################################################################################################

import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd

animate = True

####  IMPORTACIÓN DATOS EPIDEMIOLOGICOS COVID-19  #####
file = pd.read_csv('/home/juan_aku/Escritorio/AC-Covid/Covid_Data/time_series_covid19_confirmed_global.csv')
data = file.to_numpy()
####  DATOS POR PAIS:  ITALIA  ####
dataPais = data[137][13:]         #PRIMEROS 100 DATOS      

# Modelo epidemiologico SIR con autómatas celulares # 
R0 = 2.68                               #Obtenido de la literatura
#R0 = 3.0
Poblacion = 60.54e6                      #Población Italia
primeros_infectados = dataPais[0]        #Primer reporte de infectados
Irepresentativo = 1/Poblacion            #Fraccioón representativa de un ciudadano infectado
I0 = primeros_infectados/Poblacion  
S0 = 1 - I0    
brote_infectados = dataPais[21]
I21 = brote_infectados/Poblacion
S21 = 1 - I21
accionbrote_infectados = dataPais[36]      #Infectados despues de 15 días del primer reporte
I36 = accionbrote_infectados/Poblacion
pasos_diaT = 10                         #Pasos de simulación por dia 
dias_accion = 60                        #Días de acción gubernamental

#alfa = 0.5                             #Proporción de acción para mitigar la pandemia   

L = 500                      #Dimensiones
NCeldas = L*L                #Numero total celdas

I = np.zeros((L,L))          #Definición de Infectados
S = np.ones((L,L))           #Definición de Susceptibles 
R = np.zeros((L,L))          #Definición de Recuperados
It = np.zeros((L,L))         #Infectados en tiempo t

Densidad = 'Normal_Sigma=L|2,L|4.47'
#Densidad = 'Constante=1'
N = np.zeros((L,L))

sigma = L/2
miu = int(L/4)
normal = lambda coordenada : (1/np.sqrt(2*np.pi*sigma))*np.exp(-((coordenada-miu)**2)/(2*sigma**2))
for i in range(L):
    for j in range(L):
        N[i][j] = normal(i)*normal(j)*5/6 + N[i][j]

sigma = L/(np.sqrt(20))
miu = int(3*L/4)
for i in range(L):
    for j in range(L):
        N[i][j] = normal(i)*normal(j)/6 + N[i][j]

#N = np.full((L,L),1)

N = np.array(N)
N = abs(N)
NT = np.sum(N)
N = N/NT

NP = N * Poblacion             #Matriz de Poblaciones

Conectividad = 'Constante=0.15'
conect0 = 0.15
dia0 = 0
rec = (1/(15*pasos_diaT*(R0-1)))*(np.log((S21*I36)/(I21*(1-I36))))    #Coeficiente de recuperación / Calculado
v =  (rec*R0)/(1+4*conect0)                #Virulencia de la epidemia / Calculada a partir de R0 y rec
print((rec*R0)/(1+4*conect0))
#v =( np.sqrt((1+conect0-rec*R0-rec)**2 + 4*rec*R0*(1-rec)) - (1+conect0-rec*R0-rec))/2
print(v)

accion = 'nula'
conect0 = 0.15
conect1 = 0.15
deltac = (conect0 - conect1)/(dia0 - (dias_accion*pasos_diaT))

#accion = 'sigmoide'
#dc = 15
#alfa = alfa = (1/dc)*(np.log((conect0-conect1)/(1-conect0)))
#sigmoide = lambda d : (1-conect1)/(1+np.exp(-alfa*(d-dc)))

conect = conect0
con = np.full((L,L),conect)      #Coeficiente de infección de celdas vecinas

nsteps = 0                     #Numero de pasos

########  Localización primer contagio  ######
#centro = int(L/2)               #Indice del centro de la red
centro = int(L/4)

I0 = I0 / N[centro][centro]     #Calculo de proporción representativa del infectado en la celda
I[centro][centro] = I0          #Asignación de individuo infectado
S[centro][centro] = 1 - I0
Infectados = I0 * (Poblacion*N[centro][centro])
Susceptibles = Poblacion - Infectados
Recuperados = 0

k = 0
directorio = 'animacion_L=' + str(L) + '_Normal_Sigma=L|2,L|4.47'
if not(os.path.isdir(directorio)):
    os.mkdir(directorio)

archivo         = 'AC_D:%s_C:%s_A:%s_L=%d,.txt' %(Densidad,Conectividad,accion,L)
observables = open(archivo, 'w')
observables.write('%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n' %(nsteps, Infectados, Susceptibles, I0, Recuperados, I0))

NInfectados = list()               #Lista de nuevos infectados
NInfectados.append(I0)
ultimo = False

Dif = con * v
Difusion = np.multiply(N,Dif)

Is0 = list()                       #Lista de infectados en la celda central
Is0.append(I[centro][centro])

Is = list()                        #Lista de infectados acumulados
Is.append(primeros_infectados)                     

bandera1 = True
bandera2 = True

while(not(ultimo)):
    #Matrices de infectados y suceptibles temporal#

    Derecha = np.zeros((L,L))
    Arriba = np.zeros((L,L))
    Izquierda = np.zeros((L,L))
    Abajo =  np.zeros((L,L))
    #Recorrido de red por cada paso temporal#

    C_0 = np.multiply(Difusion,I)

    Derecha[:,0:L-1] = C_0[:,1:L]
    Arriba[1:L,:] = C_0[0:L-1,:]
    Izquierda[:,1:L] = C_0[:,0:L-1]
    Abajo[0:L-1,:] = C_0[1:L,:]

    C_1 = Derecha + Arriba + Izquierda + Abajo
    C = np.divide(C_1,N)

    InfeccionVecinos = np.multiply(S,C)
    InfeccionInterna_0 = np.multiply(I,S)
    InfeccionInterna = InfeccionInterna_0 * v
    NuevosInfectados = InfeccionVecinos + InfeccionInterna
    NuevosInfectados = abs(NuevosInfectados)

    R = R + I * rec                  #Recuperados acumulados
    I = I + NuevosInfectados         #Infectados acumulados
    S = S - NuevosInfectados         #Susceptibles
    It = I - R

    #print(I)
    #print("*****************")
    #print(R)

    NInfectados.append(np.sum(NuevosInfectados))

    for j in range(L):
        for i in range(L):

            if(I[i][j] > 1.0):
                I[i][j] = 1.0
                S[i][j] = 0.0

    I_celda = np.multiply(I,NP)
    Infectados = np.sum(I_celda)
    Susceptibles = Poblacion - Infectados

    if ( bandera2 and (nsteps/pasos_diaT >= 21)):
        coordL = int(3*L/4)    #Ubicación infectados de Lombardía
        I[coordL][coordL] = I[coordL][coordL] + (16*Irepresentativo/N[coordL][coordL])
        bandera2 = False        

    Is.append(Infectados)
    print(Infectados)
    if( Infectados >= dataPais[-1]):
        ultimo = True

    nsteps = nsteps + 1
    Is0.append(I[centro][centro])

    observables.write('%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n' %(nsteps, Infectados, Susceptibles, NInfectados[nsteps], np.sum(R), np.sum(It)))

    #############   SECCIÓN DE CONTROL GUBERNAMENTAL   #########

    if( bandera1 and (nsteps/pasos_diaT >= 36) and (deltac != 0)):
        #conect = (alfa*rec)/(4*v) - 0.25
        dia0 = dia0 + 1
        conect = conect0 + (deltac*dia0)
        #conect = 1 - sigmoide(dia0)
        print("******************")
        print(nsteps/pasos_diaT)
        print('Conectividad = %.4f' %conect)
        con = np.full((L,L),conect)
        Dif = con * v
        Difusion = np.multiply(N,Dif)
        if(dia0 >= dias_accion*pasos_diaT):
            bandera1 = False

    if animate:
        if(not(nsteps%10)):
            #Mostrar la malla
            plt.clf()
            plt.imshow(I*10e4, cmap='OrRd', clip_on=True )
            plt.xticks([])
            plt.yticks([])
            plt.title('%d x %d Covid-AC model, paso = %d' %(L,L,nsteps))
            plt.pause(0.01)
            filename = directorio+'/'+'step'+str(k)+'.png'
            plt.savefig(filename)
            k=k+1

print(conect)
pasos_tiempo = len(Is)/len(dataPais)    #Correspondencia entera entre pasos de simulación y tiempo