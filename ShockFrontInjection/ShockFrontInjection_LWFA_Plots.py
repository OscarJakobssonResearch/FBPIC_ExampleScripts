import os
import math
import os.path
import numpy as np
from os.path import basename
import matplotlib.pyplot as plt
from opmd_viewer import OpenPMDTimeSeries
from opmd_viewer.addons import LpaDiagnostics
from scipy.ndimage.filters import gaussian_filter

ts_laser = LpaDiagnostics('./diags/hdf5/',check_all_files=False)
ts_particles = OpenPMDTimeSeries('./diags/hdf5/',check_all_files=False)

fields={'cmap':'coolwarm'}
particles={'cmap':'jet'}


SimulationName='PulseTrainSimulation'

#Make directories for plots
Path='Plots'+SimulationName
if not os.path.exists(Path):
    os.makedirs(Path)
if not os.path.exists(Path+'/EzLineout'):
    os.makedirs(Path+'/EzLineout')


if not os.path.exists(Path+'/LaserEnvelope'):
    os.makedirs(Path+'/LaserEnvelope')   
if not os.path.exists(Path+'/LaserEnvelope/Lineout'):
    os.makedirs(Path+'/LaserEnvelope/Lineout')
if not os.path.exists(Path+'/rho'):
    os.makedirs(Path+'/rho')
      
if not os.path.exists('TextFiles'):
    os.makedirs('TextFiles')
if not os.path.exists(Path+'/Ez'):
    os.makedirs(Path+'/Ez')    
if not os.path.exists(Path+'/uz'):
    os.makedirs(Path+'/uz') 
# The main loop will run for every i*factor,
# from 0 to ts.iterations.size/factor
# Setting factor=2 will process every other file and so on.
factor = 1
a0PerSnapshot=[]
PulsewaistPerSnapshot=[]
PulselengthPerSnapshot=[]
EzMinPerSnapshot=[]

for i in range(0,int(ts_laser.iterations.size/factor)):    

    iter = ts_laser.iterations[factor*i]
    print('Iteration:%d' %iter)

    Elaser=ts_laser.get_laser_envelope(iteration=iter, pol='x',m='all',index='all', plot=True); #freq_filter=0
    plt.savefig('%s/LaserEvelope_Iteration%09i.png' %(Path+'/LaserEnvelope', iter), bbox_inches='tight')
    plt.close()

    #1D Outline of laser envelope
    ElaserSlice=ts_laser.get_laser_envelope(iteration=iter, pol='x',slicing_dir='y',theta=0, plot=True);
    plt.savefig('%s/LineoutLaserEvelope_Iteration%09i.png' %(Path+'/LaserEnvelope/Lineout', iter), bbox_inches='tight')
    plt.close()


    Ez, info_Ez = ts_particles.get_field(iteration=iter, field='E', coord='z', plot=True,**particles ) 
    #plt.clim(-2e11,1e11)
    plt.clim(-5e10,5e10)
    plt.savefig('%s/Ez_Iteration%09i.png' %(Path+'/Ez', iter), bbox_inches='tight')
    plt.close()

    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=20, axis='y')
    plt.locator_params(nbins=20, axis='x')
    plt.plot(Ez[100])
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    #plt.xticks(rotation=45)
    
    Erestricted=[]
    for i in range(int(1.0*len(Ez[100]))):
        Erestricted.append(Ez[100][i+int(0.0*len(Ez[100]))])
    
    EzMinPerSnapshot.append(np.amax(Erestricted))
    #plt.ylim(-6e9,6e9)
    plt.savefig('%sEzLineout_%09i2.png' %(Path+'/EzLineout/', iter), bbox_inches='tight')
    plt.close()

    rho, info_rho = ts_particles.get_field(iteration=iter, field='rho', plot=True,**particles)
    plt.savefig('%s/rho_%09i.png' %(Path+'/rho', iter), bbox_inches='tight')
    plt.close()

    a0=ts_laser.get_a0(iteration=iter, pol='x')
    #print(a0)
    a0PerSnapshot.append(a0)

    waist=ts_laser.get_laser_waist(iteration=iter, pol='x')
    #print(waist)
    PulsewaistPerSnapshot.append(waist)

    length=ts_laser.get_ctau(iteration=iter, pol='x')
    #print(length)
    PulselengthPerSnapshot.append(length)
    
    ts_particles.get_particle( ['z', 'uz'], species='electrons', iteration=iter );
    z_selected, uz_selected = ts_particles.get_particle( ['z', 'uz'], species='electrons', iteration=iter, select={'uz':[1, None]} )
    plt.plot(z_selected,uz_selected, 'b.',markersize=0.1)
    plt.savefig('%s/uz_%09i.png' %(Path+'/uz', iter), bbox_inches='tight')
    plt.close()
  
np.savetxt('TextFiles/'+'a0PerSnapshot'+'SimulationName'+'.txt',a0PerSnapshot)
np.savetxt('TextFiles/'+'PulsewaistPerSnapshot'+'SimulationName'+'.txt',PulsewaistPerSnapshot)
np.savetxt('TextFiles/'+'PulselengthPerSnapshot'+'SimulationName'+'.txt',PulselengthPerSnapshot)



np.savetxt('TextFiles/'+'EzLineout'+'SimulationName'+'.txt',EzMinPerSnapshot)
EzLineout = []
with open('TextFiles/'+'EzLineout'+'SimulationName'+'.txt') as data: 
    datalines = (line.rstrip('\r\n') for line in data)
    for line in datalines:
        EzLineout.append(float(line))


fig2, axis2 = plt.subplots()
plt.locator_params(nbins=10, axis='y')
plt.locator_params(nbins=10, axis='x')
plt.plot(EzLineout, color='steelblue')
axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
#plt.xticks(rotation=45)
#plt.xlim(0,260)
#plt.ylim(0,6e9)
plt.savefig('TextFiles/'+'EzLineout.pdf')
plt.close()

with open('TextFiles/'+'a0PerSnapshot'+'SimulationName'+'.txt') as data: 
    datalines = (line.rstrip('\r\n') for line in data)
    for line in datalines:
        a0PerSnapshot.append(float(line))
PulsewaistPerSnapshot = []
with open('TextFiles/'+'PulsewaistPerSnapshot'+'SimulationName'+'.txt') as data: 
    datalines = (line.rstrip('\r\n') for line in data)
    for line in datalines:
        PulsewaistPerSnapshot.append(float(line))
PulselengthPerSnapshot = []
with open('TextFiles/'+'PulselengthPerSnapshot'+'SimulationName'+'.txt') as data: 
    datalines = (line.rstrip('\r\n') for line in data)
    for line in datalines:
        PulselengthPerSnapshot.append(float(line))

plt.plot(a0PerSnapshot, color='steelblue')
plt.savefig('TextFiles/'+'a0PerSnapshot.pdf')
plt.close()
plt.plot(PulsewaistPerSnapshot, color='steelblue')
plt.savefig('TextFiles/'+'PulsewaistPerSnapshot.pdf')
plt.close()
plt.plot(PulselengthPerSnapshot, color='steelblue')
plt.savefig('TextFiles/'+'PulselengthPerSnapshot.pdf')
plt.close()

names = ['/EzLineout','/Ez', '/rho','/LaserEnvelope', '/LaserEnvelope/Lineout','/uz']


for i in names:
    print(names)
    print(i)
    os.system("ffmpeg -r 10 -y -pattern_type glob -i '%s/*.png' -c:v libx264 \
        -pix_fmt yuv420p -movflags +faststart -vf 'pad=ceil(iw/2)*2:ceil(ih/2)\
        *2' -crf 5 %s.mp4" % (Path+i,Path+i))

