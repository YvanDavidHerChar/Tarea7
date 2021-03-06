import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Sacamos los datos creados en C
datos = np.genfromtxt("datos.dat")
a = datos[:,0]
b = datos[:,1]
c = datos[:,2]
d = datos[:,3]
l = datos[:,4]

#Graficamos todos los pares de parametros.
fig = plt.figure(figsize=(10,10))
ax2 = plt.subplot(6,5,6)
ax2.scatter(a, l)
ax3 = plt.subplot(6,5,11)
ax3.scatter(b,l)
ax4 = plt.subplot(6,5,12)
ax4.scatter(b,a)
ax5 = plt.subplot(6,5,16)
ax5.scatter(c,l)
ax6 = plt.subplot(6,5,17)
ax6.scatter(c,a)
ax7 = plt.subplot(6,5,18)
ax7.scatter(c,b)
ax8 = plt.subplot(6,5,21)
ax8.scatter(d,l)
ax9 = plt.subplot(6,5,22)
ax9.scatter(d,a)
ax10 = plt.subplot(6,5,23)
ax10.scatter(d,b)
ax11 = plt.subplot(6,5,24)
ax11.scatter(d,c)

pp = PdfPages('TengoSuenho.pdf')
pp.savefig(fig)
pp.close()
