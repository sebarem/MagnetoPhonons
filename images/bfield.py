from numpy import *
from scipy import *
from pylab import *
import MagnetoPhononC as M
import scipy.optimize as s
import scipy.io as sio
import sys
from pylab import *
reload(M)
fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches

fig_width = 5 # width in inches
fig_height = 5      # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 16,'font.size': 16,'legend.fontsize':16,'xtick.labelsize':16,'ytick.labelsize':16, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)


#B-field sweeps for different n
#Set the parameters
delta = 0.03
l = 4.5e-3
shift0 = 1582.5
cutoff = 100
Brange =  arange(2.4,4.0,0.01)#arange(1.0,5.5,0.02)#
e0 = M.convert(shift0,2)
n = 0.61e12
savefile = "C:/Users/Owner/research/publications/sr-papers/2013-ElectrostaticControlOfMagnetophononResonanceInGraphene/images/"

Brange = Brange**2

fig,ax1 = subplots()
ax2=ax1.twinx()

#----------------------------------------------------------------
SolPlus = []
SolMinus = []

for B in Brange:
    x1 = M.FPlus(e0,n,l,e0,B,delta,cutoff)
    x2 = M.FMinus(e0,n,l,e0,B,delta,cutoff)
        
    SolPlus.append(x1)
    SolMinus.append(x2)
    
SolPlus = array(SolPlus)
SolMinus = array(SolMinus)

#-------------------------------------------------------------------
ax1.plot(Brange,M.convert(SolPlus.real + e0,1),'-',color='red',lw=3)
ax1.plot(Brange,M.convert(SolMinus.real + e0,1),'-',color='blue',lw=3)
ax1.vlines((6.5,12.6),1580,1600,linestyle='--',color='black',lw=3)
ax2.plot(Brange,M.FillingFactor(n,Brange),color='orange',lw=3)
ax2.set_ylabel('Filling Factor',color='orange')
ax1.set_xlim((6,15))#((1.0,5.5))#
ax1.set_ylim((1580,1600))#((1505,1655))#
ax1.set_xlabel('B(T)')#('$\sqrt{B}$')#
ax1.set_ylabel('Ramanshift(cm$^{-1})$')
ax1.xaxis.set_ticks([6.0,8.0,10.0,12.0,14.0])
ax1.legend(['    n=+1','    n=-1'],loc='upper left')
fig.tight_layout()
savefig('bfield2'+'.png',dpi=1000)
show()
    
    
 