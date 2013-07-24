from scipy import *
from numpy import *
from pylab import *

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 3.4/246.0               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 3.4*0.5  # width in inches
fig_height =fig_width  #*2     # height in inches
fig_size = [fig_width,fig_height]
params = {'backend': 'GTKAgg','axes.labelsize': 10,'font.size': 10,'legend.fontsize': 10,'xtick.labelsize': 10,'ytick.labelsize': 10, 'text.usetex': False,'figure.figsize': fig_size}
rcParams.update(params)

#Create a nice greyscale plot of the G-Band splitting as function of the filling Factor

#Import Data
#Sample: 121
#Exp 20090309_121_raman
#cleaned for cosmics
filepath = 'C:/Users/Owner/research/graphene/data/samples/b_field_samples/121_200902/'
spectrapath = '20090305_121_raman/clean_data.txt'
data=loadtxt(filepath+spectrapath)



RangeLow = 1570#1556
RangeHigh = 1602#1626
low=argmin(abs(data[:,0]-RangeLow))
high=argmin(abs(data[:,0]-RangeHigh))


data1=data[low:high,2:40]
print data1.shape
#for i in range(len(offset)):
#    data1[:,i] = data1[:,i]-offset[i]




f,ax=subplots()
ax.imshow(data1,cmap=cm.Greys_r,extent=[-40,40,RangeLow,RangeHigh],origin='lower',aspect='auto')
for im in ax.get_images():
    im.set_clim(190,800)
ax.xaxis.set_ticks([-40,-20,0,20,40])
ax.yaxis.set_ticks([1570,1580,1590,1600])
ax.set_xlabel('Gate voltage(V)')
ax.set_ylabel('Raman shift(cm$^{-1}$)')
f.tight_layout()
f.savefig('greyscale.png',dpi=1000)
show()
show()