import numpy as np
import math
import matplotlib.pyplot as plt
etot=0.05
kd=2

stotlogA_=np.linspace(-3,3,100)
stotA_=np.power(10.,stotlogA_)

bA_=-(etot+stotA_+kd)
esA_=(-bA_-np.sqrt(np.power(bA_,2.)-4.*1.0*etot*stotA_))/(2.)
#print esA_
#fig,ax=plt.subplots()
plt.scatter(stotA_,esA_)
#plt.xscale('log')
#plt.yscale('log')
plt.show()
#scon
#r
