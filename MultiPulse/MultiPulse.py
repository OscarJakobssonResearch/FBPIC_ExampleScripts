# TODO:  Add diagnostics such that as much as possible is saved
# Add option for ion motion

"""
This is an input script to run a simulation of multi-pulse laser-wakefield acceleration using FBPIC.
The simulation window is fixed in space and the lasers pass through the window from left to right

Usage
-----
- Pulse trains with uniform pulse spacing can be simulated by using the 
UniformPulseTrainN(NumberOfPulses,UniformPulseSpacing) function inside add_laser_pulse.
E.g to simulate a pulse train with 50 pulses uniformly spaced by the plasma wavelength use:
add_laser_pulse(sim,UniformPulseTrainN(50,lambdap),method="antenna",z0_antenna=z0Antenna,v_antenna=0) 

- Pulse trains with non-uniform pulse spacing can be simulated by using the 
NonUniformPulseTrainN(NumberOfPulses) function inside add_laser_pulse. The non-uniform pulse spacing 
is defined inside the function by PulseSpacing(), which sets the central position of each laser in the train. 
This can for example be used to read in an array of pulse positions or setting a linearly increasing spacing.
E.g to simulate a pulse train with 50 pulses with a linearly increasing distance between the pulses:
add_laser_pulse(sim,NonUniformPulseTrainN(50),method="antenna",z0_antenna=z0Antenna,v_antenna=0) 
and 
 TODO: def PulseSpacing(i):
    SpacingInTime=
    return(c*SpacingInTime)

- Lasers are initialized using an antenna. This is located at the left-hand side of the simulation window 
at z=0. 

- Data is written to disk each plasma period.
This means that if the pulses are uniformly spaced by the plasma wavelength the output data will have the 
pulses at the same position in the simulation window every time.

"""

# -------
# Imports
# -------
import math
import numpy as np
from scipy.constants import c, e, m_e
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser_pulse
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser
from fbpic.lpa_utils.bunch import add_particle_bunch_gaussian
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic, \
	 set_periodic_checkpoint, restart_from_checkpoint

# ----------
# Parameters
# ----------

# Whether to use the GPU
use_cuda = True

# Order of the stencil for z derivatives in the Maxwell solver.
# Use -1 for infinite order, i.e. for exact dispersion relation in
# all direction (adviced for single-GPU/single-CPU simulation).
# Use a positive number (and multiple of 2) for a finite-order stencil
# (required for multi-GPU/multi-CPU with MPI). A large `n_order` leads
# to more overhead in MPI communications, but also to a more accurate
# dispersion relation for electromagnetic waves. (Typically,
# `n_order = 32` is a good trade-off.)
# See https://arxiv.org/abs/1611.05712 for more information.
n_order = -1

#Plasma density
n_e = 3e17*1.e6 	# Density (electrons.meters^-3)

epsilon0 = 8.854187817e-12;
omegap= math.sqrt(e*e*n_e /(m_e*epsilon0));
lambdap=2*math.pi*c/omegap

# The lasers in the pulse train
NumberOfPulses=10
UniformPulseSpacing=lambdap
a0 = 0.10		  
w0 = 40.e-6		  	# Laser waist
ctau = math.sqrt(2)*c/omegap  	# Laser duration
z0 = -lambdap 	# Laser centroid # laser will be at 179.26e-6 when snapshot is taken, this corresponds to iteration 5697 of the single pulse test 
#Position of antenna
z0Antenna=-1.e-6




# The simulation box
Nz = int(1200*((NumberOfPulses+2))/4)			# Number of gridpoints along z
zmax = 0.e-6		# Right end of the simulation box (meters)
zmin = -(NumberOfPulses+3)*lambdap 		# Left end of the simulation box (meters)
Nr = 200			# Number of gridpoints along r
rmax = 100.e-6		# Length of the box along r (meters)
Nm = 2				# Number of modes used

# The simulation timestep
dt = (zmax-zmin)/Nz/c	# Timestep (seconds)

# The particles
p_zmin = 0.e-6  	# Position of the beginning of the plasma (meters)
p_zmax = 10.e-3 	# Position of the end of the plasma (meters)
p_rmax = 120.e-6 	# Maximal radial position of the plasma (meters)
p_nz = 2			# Number of particles per cell along z
p_nr = 2		 	# Number of particles per cell along r
p_nt = 8			# Number of particles per cell along theta


# The moving window
v_window = c	   # Speed of the window

# The diagnostics and the checkpoints/restarts
dt_period = 100e-15		  # Period of the diagnostics in number of timesteps
save_checkpoints = False # Whether to write checkpoint files
checkpoint_period = 30000	 # Period for writing the checkpoints
use_restart = False		 # Whether to restart from a previous checkpoint
track_electrons = False	 # Whether to track and write particle ids


w_matched = w0
rel_delta_n_over_w2 = 1./( np.pi * 2.81e-15 * w_matched**4 * n_e )

def dens_func( z, r ) :
	n = np.ones_like(z)
	n = np.where( z >= 0., n * ( 1. + rel_delta_n_over_w2 * r**2), n )
	#n=np.where(z>0,0,0)
	return(n)

# The interaction length of the simulation (meters)
L_interact = 10e-3 # increase to simulate longer distance!
# Interaction time (seconds) (to calculate number of PIC iterations)
T_interact = ( L_interact + (zmax-zmin) ) / c
# (i.e. the time it takes for the moving window to slide across the plasma)







#Function to generate laser profile of pulse train
	# To be used in add_laser_pulse(sim,UniformPulseTrainN(N,spacing),...)
def UniformPulseTrainN(N,spacing):
	#First Pulse
	PulseTrain_profile=GaussianLaser( a0, w0, ctau/c, z0, lambda0=0.8e-6,theta_pol=0, cep_phase=0.,zf=0e-6 )
	#Successive pulses equally spaced by 'spacing'
	for i in range(1,N):
		PulseTrain_profile += GaussianLaser( a0, w0, ctau/c, z0-i*spacing,lambda0=0.8e-6,theta_pol=0, cep_phase=0.,zf=0e-6 )
	return(PulseTrain_profile)





# The electron beam
L0 = 0.e-6      # Position at which the beam should be "unfreezed"
z0Beam = -(NumberOfPulses+1)*lambdap-lambdap/4      # Initial centroid of the beam
sig_z = 10.e-6    # RMS size
sig_r = 10.e-6    # RMS size
n_emit = 0.e-6   # Normalized emittance in m.rad
gamma0 = 10.
sig_gamma = 0.
Qtot = 1.e-12   # Charge in Coulomb





# ---------------------------
# Carrying out the simulation
# ---------------------------

# NB: The code below is only executed when running the script,
# (`python lwfa_script.py`), but not when importing it (`import lwfa_script`).
if __name__ == '__main__':

	# Initialize the simulation object
	sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt, zmin=zmin,
		boundaries={'z':'open', 'r':'open'}, n_order=n_order, use_cuda=use_cuda )

	# Create the plasma electrons
	elec = sim.add_new_species( q=-e, m=m_e, n=n_e,
		dens_func=dens_func, p_zmin=p_zmin, p_zmax=p_zmax, p_rmax=p_rmax,
		p_nz=p_nz, p_nr=p_nr, p_nt=p_nt )


	#beam = add_particle_bunch_gaussian( sim, q=-e, m=m_e,
	#		sig_r=sig_r, sig_z=sig_z, n_emit=0,  #initial emittance is 0
	#		gamma0=gamma0, sig_gamma=sig_gamma,
	#		n_physical_particles=Qtot/e, n_macroparticles=5000,
	#		zf=L0, tf=(L0-z0Beam)/c, z_injection_plane=L0 )


	# Load initial fields
	# Add a laser to the fields of the simulation
	#add_laser( sim, a0, w0, ctau, z0 )
	#add_laser_pulse(sim,PulseTrain_profile,method="antenna",z0_antenna=z0Antenna,v_antenna=0) #fw_propagating=True
	add_laser_pulse(sim,UniformPulseTrainN(NumberOfPulses,UniformPulseSpacing),method="antenna",z0_antenna=z0Antenna,v_antenna=0) #fw_propagating=True

	if use_restart is False:
		# Track electrons if required (species 0 correspond to the electrons)
		if track_electrons:
			elec.track( sim.comm )
	else:
		# Load the fields and particles from the latest checkpoint file
		restart_from_checkpoint( sim )

	# Configure the moving window
	sim.set_moving_window( v=v_window )

	# Add diagnostics
	sim.diags = [ FieldDiagnostic(dt_period/sim.dt, sim.fld, comm=sim.comm),
				  ParticleDiagnostic(dt_period/sim.dt, {"electrons" : elec},
					select={"ux" : [0, None ],"uy" : [0, None ]}, comm=sim.comm )]
	# Add checkpoints
	if save_checkpoints:
		set_periodic_checkpoint( sim, 50000 )

	# Number of iterations to perform
	N_step = int(T_interact/sim.dt)

	### Run the simulation
	sim.step( N_step )
	print('')
