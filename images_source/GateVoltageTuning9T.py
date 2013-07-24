from numpy import *
from scipy import *
from pylab import *
import MagnetoPhononC as M
import scipy.optimize as s
import sys
import transport as t
reload(M)
"""
fig_size = [8.0,6.0]
params = {'backend': 'WXAgg','axes.labelsize': 20,'font.size': 15,'legend.fontsize': 15,'xtick.labelsize': 15,'ytick.labelsize': 15, 'text.usetex': False,'figure.figsize': fig_size}
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

#rapidly plot simulation over simulated gatesweep data
try:
    delta = double(sys.argv[1])
    l = double(sys.argv[2])
    shift0 = double(sys.argv[3])
    cutoff = int(sys.argv[4])
    nrange = arange(double(sys.argv[5]),double(sys.argv[6]),double(sys.argv[7]))
    e0 = M.convert(shift0,2)
    B = double(sys.argv[8])
    savefile = sys.argv[9]
except:
    print 'Usage: ', sys.argv[0], ' delta lambda energy0 cutoff n-low n-high n-step n'
    sys.exit(1)
width = 0.001

#----------------------------------------------------------------
FitDataPath='C:/Users/Owner/research/graphene/graphene_hf/results/2009_MagnetoPhonons_121/B9T-Field-Results-AllFixed.txt'
FitData = loadtxt(FitDataPath)



filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
transportpath = '20090304_121_raman/20090305_121_01h52.txt'
TransportData = loadtxt(filepath+transportpath)

Voltage = TransportData[:,0]*100

#----------------------------------------------------------------
SolPlus = []
SolMinus = []
for n in nrange:
    #x1 = s.fsolve(M.FPlus,[e0,width],args=(n,l,e0,B,delta,cutoff))
    #x2 = s.fsolve(M.FMinus,[e0,width],args=(n,l,e0,B,delta,cutoff))
    x1 = M.FPlus(e0,n,l,e0,B,delta,cutoff)
    x2 = M.FMinus(e0,n,l,e0,B,delta,cutoff)
    print n,'|',x1,'|',x2
    SolPlus.append(x1)
    SolMinus.append(x2)
SolPlus = array(SolPlus)
SolMinus = array(SolMinus)
#-------------------------------------------------------------------
figure(1)

plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,1],'o',color='0.4',ms=4)
plot(M.FillingFactor(t.Gate2Density(Voltage-2.5,300),B),FitData[:,4],'o',color='0.4',ms=4)
plot(M.FillingFactor(nrange,B),M.convert(SolPlus.real+e0,1),linestyle='-',color='red',lw=3)
plot(M.FillingFactor(nrange,B),M.convert(SolMinus.real+e0,1),'-',color='blue',lw=3)
xlim((-10,10))

vlines([-6,-2,0,2,6],1582,1600,linestyle='-',color='orange',lw=3)

ylim((1582,1590))
g=gca()
g.yaxis.set_major_formatter(FormatStrFormatter('%d'))

xlabel('FillingFactor')
ylabel('Ramanshift(cm$^{-1}$)')
tight_layout()
#savefig('GateVoltageTuning2.png',dpi=1000)

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
show()