import os
import math
import os.path
import numpy as np
from os.path import basename
import matplotlib.pyplot as plt
from opmd_viewer import OpenPMDTimeSeries
from opmd_viewer.addons import LpaDiagnostics
from scipy.ndimage.filters import gaussian_filter
from scipy import fftpack
from scipy.signal import savgol_filter

from mpl_toolkits import mplot3d
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

ts_laser = LpaDiagnostics('./diags/hdf5/',check_all_files=False)
ts_particles = OpenPMDTimeSeries('./diags/hdf5/',check_all_files=False)

fields={'cmap':'coolwarm'}
particles={'cmap':'jet'}


SimulationName='TrailingPulse_PulseTrain'

#Make directories for plots
Path='Plots_'+SimulationName
if not os.path.exists(Path):
    os.makedirs(Path)
if not os.path.exists(Path+'/EzLineout'):
    os.makedirs(Path+'/EzLineout')
if not os.path.exists(Path+'/ExLineout'):
    os.makedirs(Path+'/ExLineout')
if not os.path.exists(Path+'/ExLineoutTrailingPulse'):
    os.makedirs(Path+'/ExLineoutTrailingPulse')
if not os.path.exists(Path+'/LaserEnvelope'):
    os.makedirs(Path+'/LaserEnvelope')   
if not os.path.exists(Path+'/LaserFrequency'):
    os.makedirs(Path+'/LaserFrequency')     
if not os.path.exists(Path+'/LaserFrequencyTrailingPulse'):
    os.makedirs(Path+'/LaserFrequencyTrailingPulse')
if not os.path.exists(Path+'/LaserEnvelope/Lineout'):
    os.makedirs(Path+'/LaserEnvelope/Lineout')
if not os.path.exists(Path+'/rho'):
    os.makedirs(Path+'/rho')     
if not os.path.exists('TextFiles'):
    os.makedirs('TextFiles')
if not os.path.exists(Path+'/Ez'):
    os.makedirs(Path+'/Ez')  
if not os.path.exists(Path+'/Ex'):
    os.makedirs(Path+'/Ex')       
if not os.path.exists(Path+'/Er'):
    os.makedirs(Path+'/Er')  
if not os.path.exists(Path+'/gamma_z'):
    os.makedirs(Path+'/gamma_z')  
if not os.path.exists(Path+'/y_z'):
    os.makedirs(Path+'/y_z')  

# The main loop will run for every i*factor,
# from 0 to ts_particles.iterations.size/factor
# Setting factor=2 will process every other file and so on.
factor = 1
a0PerSnapshot=[]
PulsewaistPerSnapshot=[]
PulselengthPerSnapshot=[]
EzMinPerSnapshot=[]
MaxFrequencyTrailingPulse=[]


#for i in range(0,int(ts_laser.iterations.size/factor)):    
for i in range(10800,10801): 
    iter =10800
    #iter = ts_laser.iterations[factor*i]
    print('Iteration:%d' %iter)
    
    
    Ez, info_Ez = ts_particles.get_field(iteration=iter, field='E', coord='z', plot=True,**particles ) 
    #plt.clim(-5e8,5e8)
    plt.savefig('%s/Ez_Iteration%09i.png' %(Path+'/Ez', iter), bbox_inches='tight')
    plt.close()


    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=20, axis='y')
    plt.locator_params(nbins=20, axis='x')
    plt.plot(Ez[600])
    #axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    #axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    plt.savefig('EzLineout_%09i.png' %(iter), bbox_inches='tight')
    plt.close()

    EzTimesteps = np.linspace(info_Ez.z[0]-info_Ez.z[0], info_Ez.z[len(info_Ez.z)-1]-info_Ez.z[0], len(info_Ez.z))
    EzTimesteps = [ x/299792458 for x in EzTimesteps]

    np.savetxt('EzLineout_%09i.txt'%(iter),Ez[600])
    np.savetxt('EzLineout_%09i_Timesteps.txt'%(iter),EzTimesteps)  
    plt.plot(EzTimesteps,Ez[600])

    plt.savefig('EzLineoutReversed_%09i.png' %(iter), bbox_inches='tight')    

    plt.close()

    rho, info_rho = ts_particles.get_field(iteration=iter, field='rho')
    plt.plot(rho[600])
    plt.savefig('%s/rho_%09i.png' %(Path+'/rho', iter), bbox_inches='tight')
    plt.close()
