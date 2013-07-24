from scipy import *
from numpy import *
from pylab import *

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 3.4*0.5  # width in inches
fig_height =fig_width #*2 # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 10,'font.size': 10,'legend.fontsize': 10,'xtick.labelsize': 10,'ytick.labelsize': 10, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)


#Create a nice greyscale plot of the G-Band splitting as function of the filling Factor

#Import Data
#Sample: 121
#Exp 20090309_121_raman
#cleaned for cosmics
filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
spectrapath = '20090309_121_raman/clean_data.txt'
data=loadtxt(filepath+spectrapath)

#Import Data from double lorentzian fits to substract offsets. 
#mask1 was already applied to this data so it doesn't have to be filtered again
FitData = loadtxt('C:/Users/Owner/research/graphene/graphene_hf/results/2009_MagnetoPhonons_121/EXP20090309_121_raman-12_6T-Fit_DoubleLorentz.txt')
offset = 1540*FitData[:,6]+FitData[:,7]

RangeLow = 1548#1530
RangeHigh = 1580#1600
low=argmin(abs(data[:,0]-RangeLow))
high=argmin(abs(data[:,0]-RangeHigh))

mask1 = array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,19,20,21,24,25,26,27,30,
               31,32,33,36,37,38,39,42,43,44,47,48,49,52,53,54,55,58,59,60,61,62,
               63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79])
data1=data[low:high,mask1+1]
for i in range(len(offset)):
    data1[:,i] = data1[:,i]-offset[i]
f,ax=subplots()
ax.imshow(data1,cmap=cm.Greys_r,extent=[-40,40,1570,1602],origin='lower',aspect='auto')
for im in ax.get_images():
    im.set_clim(0,450)
ax.xaxis.set_ticks([-40,-20,0,20,40])
ax.yaxis.set_ticks([1570,1580,1590,1600])
ax.set_xlabel('Gate voltage(V)')
ax.set_ylabel('Raman shift(cm$^{-1}$)')
f.tight_layout()
f.savefig('split.png',dpi=1000)
show()