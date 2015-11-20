import numpy as np
import matplotlib.pyplot as plt
import sys

#Sacamos los datos creados en C
datos = np.genfromtxt("datos.dat")
a = datos[:,0]
b = datos[:,1]
c = datos[:,2]
d = datos[:,3]
l = datos[:,4]

#Sacamos los datos del archivo
observado = np.genfromtxt("monthrg.dat")
n = 1995-1800+1
x_obs = np.linspace(0, n, n)
y_obs = np.zeros(n)

for i in range(n):
    sum = 0
    for j in range(12):
        sum = observado[2280+(i*12+j),3] + sum
    y_obs[i] = sum

#Encontramos los mejores parametros
masPosible = np.argmax(l)
mejor_a = a[masPosible]
mejor_b = b[masPosible]
mejor_c = c[masPosible]
mejor_d = d[masPosible]

#Calculamos la aproximacion
y_calulado = mejor_a*np.cos((2*np.pi/mejor_d)*x_obs + mejor_b) + mejor_c

#Creamos un arcivo de texto que denota el tiempo en el make
papitas = open('marcador.txt', "w")
papitas.write("Hola que mas, no dormi mucho...\n")
papitas.write("voy a dormir")
papitas.close()
#Graficamos la aproximacion con los datos, y tambien todos los pares de parametros.
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(6,5,21)
ax1.plot(x_obs, y_obs)
ax1.plot(x_obs, y_calulado)
ax2 = plt.subplot(6,5,1)
ax2.scatter(a, l)
ax3 = plt.subplot(6,5,6)
ax3.scatter(b,l)
ax4 = plt.subplot(6,5,7)
ax4.scatter(b,a)
ax5 = plt.subplot(6,5,11)
ax5.scatter(c,l)
ax6 = plt.subplot(6,5,12)
ax6.scatter(c,a)
ax7 = plt.subplot(6,5,13)
ax7.scatter(c,b)
ax8 = plt.subplot(6,5,16)
ax8.scatter(d,l)
ax9 = plt.subplot(6,5,17)
ax9.scatter(d,a)
ax10 = plt.subplot(6,5,18)
ax10.scatter(d,b)
ax11 = plt.subplot(6,5,19)
ax11.scatter(d,c)

plt.show()
