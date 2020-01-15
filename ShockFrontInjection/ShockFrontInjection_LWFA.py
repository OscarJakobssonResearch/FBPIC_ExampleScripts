# -------
# Imports
# -------
import numpy as np
from scipy.constants import c, e, m_e
import matplotlib.pyplot as plt
# Import the relevant structures in FBPIC
from fbpic.main import Simulation
from fbpic.lpa_utils.laser import add_laser
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

# The simulation box
Nz = 2400		  # Number of gridpoints along z
zmax = 0.	   # Right end of the simulation box (meters)
zmin = -120.e-6	  # Left end of the simulation box (meters)
Nr = 100			 # Number of gridpoints along r
rmax = 110.e-6	  # Length of the box along r (meters)
Nm = 2			 # Number of modes used

# The simulation timestep
dt = (zmax-zmin)/Nz/c	# Timestep (seconds)

# The particles
p_zmin = 0.e-6  # Position of the beginning of the plasma (meters)
p_zmax = 1000.e-6 # Position of the end of the plasma (meters)
p_rmax = 120.e-6 # Maximal radial position of the plasma (meters)
n_e = 3.0e23 # Density (electrons.meters^-3)
p_nz = 2		 # Number of particles per cell along z
p_nr = 2		 # Number of particles per cell along r
p_nt = 8		 # Number of particles per cell along theta


a0 = 2.3         # Laser amplitude
w0 = 40.e-6      # Laser waist
ctau = 13.5e-6	   # Laser duration
z0 = -30.e-6	# Laser centroid
zfoc = 0.e-6
lambda0=0.8e-6



# The moving window
v_window = c	   # Speed of the window

# The diagnostics and the checkpoints/restarts
dt_period = p_zmax/c/100		  # Period of the diagnostics in number of timesteps
save_checkpoints = False # Whether to write checkpoint files
use_restart = False		 # Whether to restart from a previous checkpoint
track_electrons = False	 # Whether to track and write particle ids

w_matched = 40.e-6
rel_delta_n_over_w2 = 1./( np.pi * 2.81e-15 * w_matched**4 * n_e )

ramp_up_start = 0.e-6   #+ 35.e-6 added such that ramp starts at zero
ramp_up_length = 0.5e-3
ramp_down_start = 0.5e-3


ramp_down_length = 0.05e-3
jet_dens_fac = 1


m_down = -1/ramp_down_length
c_down = (-ramp_down_start-ramp_down_length)*m_down

#Density function
def dens_func( z, r ) :
    n = np.ones_like(z)
    # Transverse guiding parabolic profile
    n = np.where( z >= p_zmin, n * ( 1. + rel_delta_n_over_w2 * r**2), n )
    # Make linear ramp
    # Make ramp up
    n = np.where( (z< ramp_up_start+ramp_up_length) & (z >= ramp_up_start), n+jet_dens_fac*(z-ramp_up_start)/ramp_up_length, n )
    # n = np.where( z<ramp_up_start+ramp_up_length, (z-ramp_up_start)/ramp_up_length, n )
    # Make ramp down
    n = np.where( (z >= ramp_down_start) & \
                  (z < ramp_down_start+ramp_down_length), n + jet_dens_fac*(m_down*z+c_down), n )
    return(n)

DensArray=[]
DistArray=[]
for i in range(-100,10000):
    Dist=i*1e-6
    DistArray.append(Dist)
    DensArray.append(dens_func(Dist,0))
plt.plot(DistArray,DensArray)
plt.savefig('Density_6a_shifted.pdf')
np.savetxt('DensityArray.txt',DensArray)
np.savetxt('DistanceArray.txt',DistArray)

# The interaction length of the simulation (meters)
L_interact = (p_zmax-p_zmin) # increase to simulate longer distance!
# Interaction time (seconds) (to calculate number of PIC iterations)
T_interact = ( L_interact + (zmax-zmin) ) / v_window
# (i.e. the time it takes for the moving window to slide across the plasma)

# ---------------------------
# Carrying out the simulation
# ---------------------------

# NB: The code below is only executed when running the script,
# (`python lwfa_script.py`), but not when importing it (`import lwfa_script`).
if __name__ == '__main__':

	# Initialize the simulation object
	sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt, zmin=zmin,
		boundaries='open', n_order=n_order, use_cuda=use_cuda )

	# Create the plasma electrons
	elec = sim.add_new_species( q=-e, m=m_e, n=n_e,
		dens_func=dens_func, p_zmin=p_zmin, p_zmax=p_zmax, p_rmax=p_rmax,
		p_nz=p_nz, p_nr=p_nr, p_nt=p_nt )

	# Load initial fields
	# Add a laser to the fields of the simulation
	add_laser( sim, a0, w0, ctau, z0, lambda0=lambda0, zf=zfoc) 

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
					select={"uz" : [1., None ]}, comm=sim.comm ) ]
	# Add checkpoints
	if save_checkpoints:
		set_periodic_checkpoint( sim, 500000 )

	# Number of iterations to perform
	N_step = int(T_interact/sim.dt)

	### Run the simulation
	sim.step( N_step )
	print('')
