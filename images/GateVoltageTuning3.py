from numpy import *
from scipy import *
from pylab import *
import MagnetoPhononD as M
import scipy.optimize as s
import sys
import transport as t
reload(M)

fig_size = [3.0,2.5]
#fig_size =[5.0,3.75]
params = {'backend': 'WXAgg','axes.labelsize': 8,'font.size': 8,'legend.fontsize': 8,'xtick.labelsize': 8,'ytick.labelsize': 8, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)
"""
fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 8,'font.size': 8,'legend.fontsize': 8,'xtick.labelsize': 8,'ytick.labelsize': 8, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)
"""
#rapidly plot simulation over simulated gatesweep data
try:
    delta = double(sys.argv[1])
    l = double(sys.argv[2])
    shift0 = double(sys.argv[3])
    cutoff = int(sys.argv[4])
    nrange = arange(double(sys.argv[5]),double(sys.argv[6]),double(sys.argv[7]))
    e0 = M.convert(shift0,2)
    B = double(sys.argv[8])
    vF = double(sys.argv[9])
    l2 = double(sys.argv[10])
except:
    print 'Usage: ', sys.argv[0], ' delta l shift0 cutoff n-low n-high n-step B vF'
    sys.exit(1)
width = 0.001
savedir = './'

#----------------------------------------------------------------
FitDataPath='C:/Users/Owner/research/graphene/graphene_hf/results/2009_MagnetoPhonons_121/EXP20090309_121_raman-12_6T-Fit_DoubleLorentz_Intraband.txt'
FitData = loadtxt(FitDataPath)

mask1 = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,20,21,24,25,26,27,30,31,32,33,36,37,38,39,42,43,44,47,48,49,52,53,54,55,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79])

filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
transportpath = '20090309_121_raman/20090310_121_06h56.txt'
TransportData = loadtxt(filepath+transportpath)

Voltage = TransportData[:,0]*100
my_temp = zeros(len(mask1))
for i in range(len(mask1)): my_temp[i] = Voltage[mask1[i]]
Voltage = my_temp
#----------------------------------------------------------------
SolPlus = M.FPlus(e0,nrange,l,e0,B,delta,cutoff,vF,l2)
SolMinus = M.FMinus(e0,nrange,l,e0,B,delta,cutoff,vF,l2)

SolPlusAllowed = M.FPlus(e0,nrange,l,e0,B,delta,cutoff,vF,0.0)
SolMinusAllowed = M.FMinus(e0,nrange,l,e0,B,delta,cutoff,vF,0.0)

#-------------------------------------------------------------------
f,ax=subplots()

ax.plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,1],'o',color='0.4',ms=6)
ax.plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,4],'o',color='0.4',ms=6)
ax.plot(M.FillingFactor(nrange,B),M.convert(SolPlus.real+e0,1),linestyle='-',color='red',lw=3)
ax.plot(M.FillingFactor(nrange,B),M.convert(SolMinus.real+e0,1),'-',color='blue',lw=3)



ax.set_xlim((-10,10))
ax.set_xticks([-10,-6,-2,2,6,10])

ax.vlines([-6,-2,0,2,6],1582,1600,linestyle='-',color='orange',lw=3)

ax.set_ylim((1582,1600))
ax.set_xlabel('FillingFactor')
ax.set_ylabel('Ramanshift(cm$^{-1}$)')
f.tight_layout()
f.savefig(savedir+'GateVoltageTuning-POSITION.png',dpi=300)

ax.fill_between(M.FillingFactor(nrange,B),M.convert(SolPlus.real+e0,1),M.convert(SolPlusAllowed.real+e0,1),facecolor='red',alpha=0.5)
ax.fill_between(M.FillingFactor(nrange,B),M.convert(SolMinus.real+e0,1),M.convert(SolMinusAllowed.real+e0,1),facecolor='blue',alpha=0.5)
f.savefig(savedir+'GateVoltageTuning-POSITION-SymmetricTransitions.png',dpi=300)

f2,ax2=subplots()

ax2.plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,1],'o',color='0.4',ms=6)
ax2.plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,4],'o',color='0.4',ms=6)
ax2.plot(M.FillingFactor(nrange,B),M.convert(SolPlusAllowed.real+e0,1),linestyle='-',color='red',lw=3)
ax2.plot(M.FillingFactor(nrange,B),M.convert(SolMinusAllowed.real+e0,1),'-',color='blue',lw=3)

ax2.plot(M.FillingFactor(nrange,B),M.convert(SolPlus.real+e0,1),linestyle='--',color='red',lw=2)
ax2.plot(M.FillingFactor(nrange,B),M.convert(SolMinus.real+e0,1),'--',color='blue',lw=2)


ax2.set_xlim((-10,10))
ax2.set_xticks([-10,-6,-2,2,6,10])

ax2.vlines([-6,-2,0,2,6],1582,1600,linestyle='-',color='orange',lw=3)

ax2.set_ylim((1582,1600))
ax2.set_xlabel('FillingFactor')
ax2.set_ylabel('Ramanshift(cm$^{-1}$)')
f2.tight_layout()
f2.savefig(savedir+'GateVoltageTuning-POSITION-NO-Symmetric-transitions.png',dpi=300)

"""
figure(2)
plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,2],'o',color='0.4',ms=4)
plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,5],'o',color='0.4',ms=4)
plot(M.FillingFactor(nrange,B),(2*abs(M.convert(SolPlus.imag,1)))+7,'-',color='red',lw=4)
plot(M.FillingFactor(nrange,B),(2*abs(M.convert(SolMinus.imag,1)))+7,'-',color='blue',lw=4)
vlines([-6,-2,0,2,6],6,18,linestyle='-',color='orange',lw=3)
xlim((-10,10))
ylim((6,18))
xlabel('FillingFactor')
ylabel('FWHM(cm$^{-1}$)')
tight_layout()
"""
show()