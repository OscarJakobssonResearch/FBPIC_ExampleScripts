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
factor = 2

a0PerSnapshot=[]
PulsewaistPerSnapshot=[]
PulselengthPerSnapshot=[]
EzMinPerSnapshot=[]



x, y, z, gamma = ts_particles.get_particle(var_list=['x', 'y', 'z', 'gamma'],
                                 iteration=ts_particles.iterations[0], species='beam')

U0 = gamma.sum()
gamma_ave = gamma.mean()
sig_gamma = gamma.std()
gamma0 = gamma[0]
delta_gamma = gamma0 - gamma_ave


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
    #plt.clim(-5e8,5e8)
    plt.savefig('%s/Ez_Iteration%09i.png' %(Path+'/Ez', iter), bbox_inches='tight')
    plt.close()

    Er, info_Er = ts_particles.get_field(iteration=iter, field='E', coord='r', plot=True)
    plt.savefig('%s/Er_Iteration%09i.png' %(Path+'/Er', iter), bbox_inches='tight')
    plt.close()


    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=20, axis='y')
    plt.locator_params(nbins=20, axis='x')
    plt.plot(Ez[200])
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    plt.savefig('%sEzLineout_%09i2.png' %(Path+'/EzLineout/', iter), bbox_inches='tight')
    plt.close()

    rho, info_rho = ts_particles.get_field(iteration=iter, field='rho', plot=True,**particles)
    plt.savefig('%s/rho_%09i.png' %(Path+'/rho', iter), bbox_inches='tight')
    plt.close()

    x, y, z, gamma = ts_particles.get_particle(var_list=['x', 'y', 'z', 'gamma'],
                                     iteration=iter, species='beam')
    r_sim = np.sqrt(x**2.+y**2.)

    plt.plot(z,gamma/gamma0, linestyle=' ', marker='o')
    plt.savefig('%s/gamma_z_Iteration%09i.png' %(Path+'/gamma_z', iter), bbox_inches='tight')
    plt.close()     

    plt.plot(z,y, linestyle=' ', marker='o')
    plt.savefig('%s/y_z_Iteration%09i.png' %(Path+'/y_z', iter), bbox_inches='tight')
    plt.close() 

    '''
    ts_particles.get_particle( ['z', 'uz'], species='electrons', iteration=iter );
    z_selected, uz_selected = ts_particles.get_particle( ['y', 'uy'], species='electrons', iteration=iter, select={'uy':[0, None]} )
    plt.plot(z_selected,uz_selected, 'b.',markersize=0.1)
    plt.savefig('%s/uy_%09i.png' %(Path+'/uy', iter), bbox_inches='tight')
    plt.close()    
    '''


names = ['/EzLineout','/Ez','/Er', '/rho', '/gamma_z', '/y_z','/LaserEnvelope', '/LaserEnvelope/Lineout']


for i in names:
    print(names)
    print(i)
    os.system("ffmpeg -r 10 -y -pattern_type glob -i '%s/*.png' -c:v libx264 \
        -pix_fmt yuv420p -movflags +faststart -vf 'pad=ceil(iw/2)*2:ceil(ih/2)\
        *2' -crf 5 %s.mp4" % (Path+i,Path+i))

