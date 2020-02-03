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
    
    Elaser=ts_laser.get_laser_envelope(iteration=iter, pol='x',m='all',index='all', plot=True); #freq_filter=0
    plt.savefig('%s/LaserEvelope_Iteration%09i.png' %(Path+'/LaserEnvelope', iter), bbox_inches='tight')
    plt.close()

    #1D Outline of laser envelope
    ElaserSlice=ts_laser.get_laser_envelope(iteration=iter, pol='x',slicing_dir='y',theta=0, plot=True);
    plt.savefig('%s/LineoutLaserEvelope_Iteration%09i.png' %(Path+'/LaserEnvelope/Lineout', iter), bbox_inches='tight')
    plt.close()
 


    Laser3D,Laser3DInfo=ts_laser.get_laser_envelope(iteration=iter, pol='x',m='all',index='all', plot=True);

    print(Laser3D)
    print(Laser3DInfo.z)
    print(Laser3DInfo.z[0])
    print(Laser3DInfo.z[1])        
    print(Laser3DInfo.z[1] - Laser3DInfo.z[0]) 
    

    y = np.linspace(Laser3DInfo.z[0], Laser3DInfo.z[len(Laser3DInfo.z)-1], len(Laser3DInfo.z))
    x = np.linspace(Laser3DInfo.r[0], Laser3DInfo.r[len(Laser3DInfo.r)-1], len(Laser3DInfo.r))

    X, Y = np.meshgrid(y, x)
    Z = Laser3D

    #black background color:  https://stackoverflow.com/questions/51107968/change-3d-background-to-black-in-matplotlib
    plt.style.use('dark_background')

    fig = plt.figure(figsize=(20,5))  #figsize=(40,5)
    ax = fig.gca(projection='3d')   #ax = plt.axes(projection='3d')
    #ax.set_aspect('equal')
    #ax.contour3D(X, Y, Z, 50, cmap='binary')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
    #plt.xlim(700e-6,900e-6)
    #plt.savefig('test.pdf')
    ax.view_init(30, -90)

    
    
    #fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    ax.plot_surface(X, Y, Z,cmap='plasma', edgecolor='none',rcount=100, ccount=100)
    fig.savefig('test1.png', bbox_inches='tight')
    

    
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min()]).max()    
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    for xb, yb in zip(Xb, Yb):
       ax.plot([xb], [yb], 'w')
    plt.grid()
    #plt.show()
    
    #remove axes   https://stackoverflow.com/questions/45148704/how-to-hide-axes-and-gridlines-in-matplotlib-python
    plt.axis('off')
    plt.grid(b=None)

    ax.set_xlim(np.array([Laser3DInfo.z[0]+(-Laser3DInfo.z[0]+Laser3DInfo.z[len(Laser3DInfo.z)-1])/4 ,Laser3DInfo.z[len(Laser3DInfo.z)-1]-(-Laser3DInfo.z[0]+Laser3DInfo.z[len(Laser3DInfo.z)-1])/4]))

    #bbox = fig.bbox_inches.from_bounds(1, 1, 8, 6)
    fig.savefig('test2.png', bbox_inches='tight')
    
    
    
    def rotate(angle):
        ax.view_init(30,azim=angle)
        ax.set_xlim(np.array([Laser3DInfo.z[0]+(-Laser3DInfo.z[0]+Laser3DInfo.z[len(Laser3DInfo.z)-1])/4 ,Laser3DInfo.z[len(Laser3DInfo.z)-1]-(-Laser3DInfo.z[0]+Laser3DInfo.z[len(Laser3DInfo.z)-1])/4]))

    rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2),interval=100)

    rot_animation.save('rotation.gif', dpi=80, writer='imagemagick')
    '''
    

    FrequencySpectrum=ts_laser.get_spectrum(iteration=iter, pol='x', plot=True);
    plt.xlim(2.1e15,2.5e15)
    plt.savefig('%s/LaserFrequency_Iteration%09i.png' %(Path+'/LaserFrequency', iter), bbox_inches='tight')
    plt.close()
    

    
    Ez, info_Ez = ts_particles.get_field(iteration=iter, field='E', coord='z', plot=True,**particles ) 
    #plt.clim(-5e8,5e8)
    plt.savefig('%s/Ez_Iteration%09i.png' %(Path+'/Ez', iter), bbox_inches='tight')
    plt.close()

    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=10, axis='y')
    plt.locator_params(nbins=10, axis='x')
    plt.plot(Ez[200])
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    plt.savefig('%sEzLineout_%09i2.png' %(Path+'/EzLineout/', iter), bbox_inches='tight')
    plt.close()

    Er, info_Er = ts_particles.get_field(iteration=iter, field='E', coord='r', plot=True)
    plt.savefig('%s/Er_Iteration%09i.png' %(Path+'/Er', iter), bbox_inches='tight')
    plt.close()


    Ex, info_Ex = ts_particles.get_field(iteration=iter, field='E', coord='x', plot=True,**particles ) 
    #plt.clim(-5e8,5e8)
    plt.savefig('%s/Ex_Iteration%09i.png' %(Path+'/Ex', iter), bbox_inches='tight')
    plt.close()

    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=10, axis='y')
    plt.locator_params(nbins=10, axis='x')
    plt.plot(Ex[200])
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    plt.savefig('%sExLineout_%09i2.png' %(Path+'/ExLineout/', iter), bbox_inches='tight')
    plt.close()


    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=10, axis='y')
    plt.locator_params(nbins=10, axis='x')
    plt.plot(Ex[200][860:1250])
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    plt.savefig('%sExLineoutTrailingPulse_%09i2.png' %(Path+'/ExLineoutTrailingPulse/', iter), bbox_inches='tight')
    plt.close()



    fig2, axis2 = plt.subplots()
    plt.locator_params(nbins=10, axis='y')
    plt.locator_params(nbins=10, axis='x')

    #plt.xlim(860,1250)
    #plt.ylim(-1.1e12,1.1e12)
    plt.xticks(rotation=90)
    axis2.xaxis.grid(linestyle='solid',color='gray',linewidth=0.1)
    axis2.yaxis.grid(linestyle='solid',color='gray',linewidth=0.1)

    Erestricted=np.zeros((len(Ez[200])))  #remove all laser pulses
    #for i in range(int(1.0*len(Ez[200]))):
    for i in range(860,1250):              # add back only the last one. This gives a smoother curve for the fft of the last pulse for some reason, comapred to doing fft only on the cropped last pulse
        Erestricted[i]=(Ex[200][i])
    #for i in range(860,1250):    
    #    Erestricted.append(Ez[200][i])
    
    X = fftpack.fft(Erestricted)
    f_s= 2*math.pi/((info_Ez.z[1] - info_Ez.z[0])/299792458)
    freqs = fftpack.fftfreq(len(Erestricted)) * f_s
    fig, ax = plt.subplots()
    plt.xlim(2.1e15,2.6e15)
    plt.plot(freqs, np.abs(X))
  
    MaxFrequencyTrailingPulse.append(freqs[np.argmax(np.abs(X))])

    #ax.set_xlabel('Frequency in Hertz [Hz]')
    #ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
    #ax.set_xlim(-f_s / 2, f_s / 2)
    #ax.set_ylim(-5, 110)

    plt.savefig('%sAngularFrequencyTrailingPulse_%09i2.png' %(Path+'/LaserFrequencyTrailingPulse/', iter), bbox_inches='tight')
    plt.close()

    rho, info_rho = ts_particles.get_field(iteration=iter, field='rho', plot=True,**particles)
    plt.savefig('%s/rho_%09i.png' %(Path+'/rho', iter), bbox_inches='tight')
    plt.close()




    y = np.linspace(info_Ez.z[0], info_Ez.z[len(info_Ez.z)-1], len(info_Ez.z))
    x = np.linspace(info_Ez.r[0], info_Ez.r[len(info_Ez.r)-1], len(info_Ez.r))

    X, Y = np.meshgrid(y, x)
    Z = Ez

    fig = plt.figure(figsize=(12,5))
    ax = plt.axes(projection='3d')
    #ax.contour3D(X, Y, Z, 50, cmap='binary')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
    #plt.xlim(700e-6,900e-6)
    #plt.savefig('test.pdf')
    ax.view_init(30, -90)
    ax.plot_surface(X, Y, Z,cmap='jet', edgecolor='none',rcount=100, ccount=100,antialiased=True)
    fig.savefig('testEz.pdf', bbox_inches='tight')
    '''


    '''
    x, y, z, gamma = ts_particles.get_particle(var_list=['x', 'y', 'z', 'gamma'],
                                     iteration=iter, species='beam')
    r_sim = np.sqrt(x**2.+y**2.)

    plt.plot(z,gamma/gamma0, linestyle=' ', marker='o')
    plt.savefig('%s/gamma_z_Iteration%09i.png' %(Path+'/gamma_z', iter), bbox_inches='tight')
    plt.close()     

    plt.plot(z,y, linestyle=' ', marker='o')
    plt.savefig('%s/y_z_Iteration%09i.png' %(Path+'/y_z', iter), bbox_inches='tight')
    plt.close() 

    
    ts_particles.get_particle( ['z', 'uz'], species='electrons', iteration=iter );
    z_selected, uz_selected = ts_particles.get_particle( ['y', 'uy'], species='electrons', iteration=iter, select={'uy':[0, None]} )
    plt.plot(z_selected,uz_selected, 'b.',markersize=0.1)
    plt.savefig('%s/uy_%09i.png' %(Path+'/uy', iter), bbox_inches='tight')
    plt.close()    
    '''
'''
plt.plot(MaxFrequencyTrailingPulse,linestyle=' ', marker='o')
plt.ylim(2.1e15,2.6e15)
#print(MaxFrequencyTrailingPulse)
plt.savefig('MaxFrequencyTrailingPulse.pdf')


names = ['/EzLineout','/ExLineout','/ExLineoutTrailingPulse','/Ez','/Er','/Ex', '/rho', '/gamma_z', '/y_z','/LaserEnvelope','/LaserFrequency', '/LaserFrequencyTrailingPulse' , '/LaserEnvelope/Lineout']


for i in names:
    print(names)
    print(i)
    os.system("ffmpeg -r 10 -y -pattern_type glob -i '%s/*.png' -c:v libx264 \
        -pix_fmt yuv420p -movflags +faststart -vf 'pad=ceil(iw/2)*2:ceil(ih/2)\
        *2' -crf 5 %s.mp4" % (Path+i,Path+i))
'''
