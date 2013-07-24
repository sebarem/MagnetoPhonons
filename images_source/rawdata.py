from scipy import *
from numpy import *
from pylab import *
from scipy.optimize import leastsq
import transport as t
import raman as r
reload(r)

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 0.5*fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*2      # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 12,'font.size': 12,'legend.fontsize': 10,'xtick.labelsize': 12,'ytick.labelsize': 12, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)


#-----------------------------------------------------
#DEFINITION OF DOUBLE LORENTIZAN FIT FUNCTION
#-----------------------------------------------------
p = array([100.0,1582.0,10.0,100.0,1588.0, 10.0, 0.0,200.0]) #initial conditions
def g(x,a,w):
    value = 1-abs(a*sin(x/w)*exp(-x**2/w**2))-abs(0.1e-13*x)
    return value
def dlorentz2mod(p,y,x,n):
    p[2]=g(n,0.5,1.0e12)*p[5]
    temp1=max(p[2],p[5])
    temp2=min(p[2],p[5])
    #p[2]=temp2
    #p[5]=temp1
    #if p[1]>=p[4]:
    #    temp=p[4]
    #    p[4]=p[1]
    #    p[1]=temp
    value=r.dlorentz2(p,y,x)
    return value
def fitting(p,x,y,voltage):
    n=t.Gate2Density(voltage)
    plsq=leastsq(dlorentz2mod,p,args=(y,x,n))
    print plsq[0]
    return plsq[0]
    
#Create a nice greyscale plot of the G-Band splitting as function of the filling Factor

#Import Data
#Sample: 121
#Exp 20090309_121_raman
#cleaned for cosmics
filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
spectrapath = '20090309_121_raman/clean_data.txt'
data=loadtxt(filepath+spectrapath)


RangeLow = 1530
RangeHigh = 1600
low=argmin(abs(data[:,0]-RangeLow))
high=argmin(abs(data[:,0]-RangeHigh))

x = data[low:high,0]+26
xfit = arange(RangeLow+26, RangeHigh+26,0.1)
y1 = data[low:high,20] #spectrum at Vgate = -20V
y2 = data[low:high,32] #spectrum at Vgate = -8V

p1=fitting(p,x,y1,-20.0)
p2=fitting(p,x,y2,-8.0)

f,ax=subplots()

ax.plot(x,y1,'o',mfc='black',markersize=5)
ax.plot(x,y2,'o',mfc='red',mec='black',markersize=5)

r.PlotDoubleLorentz(ax,xfit,p1,'2','black')
r.PlotDoubleLorentz(ax,xfit,p2,'2')


ax.set_xlabel('Raman shift(cm$^{-1}$)')
ax.set_ylabel('Intensity(a.u.)')
ax.xaxis.set_ticks([1560,1580,1600,1620])

f.tight_layout()
f.savefig('rawdata.png',dpi=1000)
show()