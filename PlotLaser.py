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


SimulationName='BeamLoading_PulseTrain'

#Make directories for plots
Path='Plots_'+SimulationName
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


'''
x, y, z, gamma = ts_particles.get_particle(var_list=['x', 'y', 'z', 'gamma'],
                                 iteration=ts_particles.iterations[0], species='beam')
U0 = gamma.sum()
gamma_ave = gamma.mean()
sig_gamma = gamma.std()
gamma0 = gamma[0]
delta_gamma = gamma0 - gamma_ave
'''

for i in range(10800,10801):    

    iter =10800 #ts_laser.iterations[factor*i]
    print('Iteration:%d' %iter)

    #fig, ax = plt.subplots()

    Elaser, ElaserInfo=ts_laser.get_laser_envelope(iteration=iter, pol='x',m='all',index='all',plot=True); #freq_filter=0
    ElaserSlice, ElaserSliceInfo=ts_laser.get_laser_envelope(iteration=iter, pol='x',slicing_dir='y',theta=0);
    

    #fig=plt.imshow(Elaser)
    
    
    #ax.imshow(Elaser, aspect='auto',zorder=0)
    plt.savefig('test2%09i.png' %(iter), bbox_inches='tight', dpi=100)
    #plt.plot(ElaserSlice)
    #plt.close()
    ax = plt.gca()
    ElaserInfo.z = [ x*1e6 for x in ElaserInfo.z]
    ElaserSlice = [ ((1e6*np.max(ElaserInfo.r)/2)*x/(np.max(ElaserSlice))-1e6*np.max(ElaserInfo.r)) for x in ElaserSlice]
    ax.plot(ElaserInfo.z,ElaserSlice, color='white',linewidth=1.3)


    plt.savefig('test%09i.png' %(iter), bbox_inches='tight', dpi=100)

    # Elaser[:,int(len(ElaserInfo.z)/2)] gives the transverse lineout at the longitudinal center of the beam
    Elaser[:,int(len(ElaserInfo.z)/2)] = [ ElaserInfo.z[0]+(np.max(ElaserInfo.z)/6)*x/(np.max(Elaser[:,int(len(ElaserInfo.z)/2)])) for x in Elaser[:,int(len(ElaserInfo.z)/2)]]

    ElaserInfo.r = [ x*1e6 for x in ElaserInfo.r]
    ax.plot(Elaser[:,int(len(ElaserInfo.z)/2)],ElaserInfo.r,color='white',linewidth=1.3)
    plt.savefig('transverse4.pdf')
    plt.close()         
