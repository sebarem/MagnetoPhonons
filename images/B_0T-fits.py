from scipy import *
from numpy import *
import scipy.io as sio
import transport as t
from scipy.optimize import leastsq
from pylab import *
import raman as r
import time

fig_width_pt = 246.0  # Get this from LaTeX using /showthe/columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 3.4*0.5  # width in inches
fig_height =fig_width #*2 # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 10,'font.size': 10,'legend.fontsize': 10,'xtick.labelsize': 10,'ytick.labelsize': 10, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)

filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
spectrapath = '20090305_121_raman/clean_data.txt'
transportpath = '20090305_121_raman/20090305_121_14h52.txt'
data=loadtxt(filepath+spectrapath)
TransportData = loadtxt(filepath+transportpath)


voltage = TransportData[:,0]*100
position = loadtxt('C:/Users/Owner/research/graphene/graphene_hf/results/2009_MagnetoPhonons_121/Combined-Positions.txt')
width = loadtxt('C:/Users/Owner/research/graphene/graphene_hf/results/2009_MagnetoPhonons_121/Combined-Width.txt')


#Now plot the model
alpha = 4.8e-3 #4.5e-3
w0 = 1582.0 #1581.4 #undisturbed Ramanshift in cm^-1 1586
diff = 0.0
e45delta = 0.010
e38delta = 0.05+diff#0.021 #0.051+2*diff#0.050 #0.0455
delta0 = 5.5
vF=1.1e6 #1.10e6
dn=0.3
energy = t.Density2Energy(t.Gate2Density(voltage+1,300),vF)


Energy = arange(-0.7,0.7,0.01) #Energy range to be plotted in eV
#n = t.Energy2Density(Energy)
#-----------------------------------------------------
e45Ando = r.SelfEnergy(Energy,alpha,w0,e45delta)
e45AndoG = delta0-2*r.eV2Raman(e45Ando.imag)
e45AndoW = w0+r.eV2Raman(e45Ando.real)
#-----------------------------------------------------
e38Ando = r.SelfEnergy(Energy,alpha,w0,e38delta)
e38AndoG = delta0-2*r.eV2Raman(e38Ando.imag)
e38AndoW = w0+r.eV2Raman(e38Ando.real)

def build_peak(energies,my_range,alpha,delta,w0,delta0):
    p1=[] #array to hold distribution of Gammas
    p2=[] #array to hold distibution of Omegas
    my_peak = zeros(len(my_range)) #range of Ramanshifts
    for m in energies:
        Andol = r.SelfEnergy(m,alpha,w0,delta)
        #print Ando
        AndoGl = delta0-2*w0*Andol.imag/0.196
        AndoWl = w0+w0*Andol.real/0.196
        p1.append(AndoGl)
        p2.append(AndoWl)
        my_peak = my_peak + r.lorentz(my_range,[2/(pi*AndoGl),AndoWl,AndoGl,0,0])
    
    p1 = array(p1)
    
    p2=array(p2)
    return my_peak


my_range=arange(1550,1620,0.1)
p0 = array([7000.0,1585.0,10.0,0.0,0.0])

n=arange(0e12,8e12,0.1e12)

p=[]
StartTime = time.clock()

for nl in n:
    print 'Building Spectrum!',time.clock()-StartTime
    my_n = nl/1e12
    N = (dn*randn(1000)+my_n)
    temp_n = N*1e12 #holds the n distribution
    temp_E = t.Density2Energy(temp_n,vF) #holds the energy distribution calculated from n distribution
    my_peak =build_peak(temp_E,my_range,alpha,e45delta,w0,delta0)
    
    plsq = leastsq(r.dlorentz,p0,args=(my_peak,my_range))
    p.append(plsq[0])
    

e45p=array(p)

f,ax = subplots()   
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
vF=1.13e6
energy = t.Density2Energy(t.Gate2Density(voltage+1,300),vF)
ax.plot(voltage+1,position,'o',color='0.4',markersize=8)
ax.set_xlabel('Gatevoltage(V)')
ax.set_ylabel('Ramanshift(cm$^{-1}$)')
ax.set_xticks([-40,-20,0.0,20,40])
#ax.vlines((-0.1962/2,0.1962/2),1580,1594,color='grey')
ax.set_xlim((-40,40))
ax.set_ylim((1580,1602))
ax.set_yticks([1580,1584,1588,1592,1596,1600])

f.tight_layout()

f.savefig('Combined-Fits-Shift-DataOnly.png',dpi=300)


ax.plot((-1)*t.Density2Gate(n,300),e45p[:,1],'--',color='red',lw=3)
ax.plot(t.Density2Gate(n,300),e45p[:,1],'--',color='red',lw=3)


f.savefig('l-48_d-10_dn-03_vf-110_w0-15820_d0-55_DLCombi-SHIFT.png',dpi=300)

show()