import numpy as np
import pickle
import sys
import os
import time

# Command line arguments
w = float(sys.argv[1]) # Anderson disorder strength
Ncells = 3*int(sys.argv[2]) # Number of unit cells along one axis
n_kgrid = int(sys.argv[3]) # Number of k-points along each axis
nrep = int(sys.argv[4]) # Number of disorder configurations
Ntot = Ncells*Ncells # Total number of unit cells
Natom = 2*Ntot # Total number of atoms
u_list = np.linspace(0, 1, n_kgrid+1)[:n_kgrid] # List of periodic boundary condition phases in b1 direction, i.e. k dot a1 / 2 pi
v_list = np.linspace(0, 1, n_kgrid+1)[:n_kgrid] # List of periodic boundary condition phases in b2 direction, i.e. k dot a2 / 2 pi
# For the purposes of convention, round both u and v to be in [-0.5, 0.5] in case they are not
u_list = ((u_list+0.5)%1)-0.5
v_list = ((v_list+0.5)%1)-0.5
print("u list: %s" % u_list)
print("v list: %s" % v_list)
# This is mathematically unnecessary due to periodicity, but normally the k-points are listed this way
cchop = 1 # carbon-carbon hopping
aCC = 1 # carbon-carbon distance
# Set the periodic boundary condition hopping to zero if we are not summing over k-points and 1 otherwise
if n_kgrid==1:
	pbchop = 0
else:
	pbchop = 1

ev_path = "en_v2"
if not os.path.exists(ev_path):
	os.makedirs(ev_path)

"""
Indexing: u=2(Ncells)m+2n+T
T = 0 for A, 1 for B
n,m = 0, 1, ... Ncells-1
"""
def ind(T,n,m):
	return 2*Ncells*m + 2*n + T

# Pristine Hamiltonian

def make_Hbulk(Ncells, cchop, aCC):
	#Hamiltonian
	Ht = np.zeros( ( 2*Ncells*Ncells , 2*Ncells*Ncells ), dtype=np.complex128 )
	for n in np.arange(Ncells-1):
		for m in np.arange(Ncells-1):
			Ht[ ind(0,n,m) , ind(1,n,m) ] = cchop
			Ht[ ind(1,n,m) , ind(0,n,m) ] = cchop
			Ht[ ind(0,n,m) , ind(1,n+1,m) ] = cchop
			Ht[ ind(1,n+1,m) , ind(0,n,m) ] = cchop
			Ht[ ind(0,n,m) , ind(1,n,m+1) ] = cchop
			Ht[ ind(1,n,m+1) , ind(0,n,m) ] = cchop
	# "Lower right" boundary
	for n in np.arange(Ncells-1):
		Ht[ ind(0,n,Ncells-1) , ind(1,n,Ncells-1) ] = cchop
		Ht[ ind(1,n,Ncells-1) , ind(0,n,Ncells-1) ] = cchop
		Ht[ ind(0,n,Ncells-1) , ind(1,n+1,Ncells-1) ] = cchop
		Ht[ ind(1,n+1,Ncells-1) , ind(0,n,Ncells-1) ] = cchop
	# "Upper right" boundary
	for m in np.arange(Ncells-1):
		Ht[ ind(0,Ncells-1,m) , ind(1,Ncells-1,m) ] = cchop
		Ht[ ind(1,Ncells-1,m) , ind(0,Ncells-1,m) ] = cchop
		Ht[ ind(0,Ncells-1,m) , ind(1,Ncells-1,m+1) ] = cchop
		Ht[ ind(1,Ncells-1,m+1) , ind(0,Ncells-1,m) ] = cchop
	# Corner
	Ht[ ind(0,Ncells-1,Ncells-1) , ind(1,Ncells-1,Ncells-1) ] = cchop
	Ht[ ind(1,Ncells-1,Ncells-1) , ind(0,Ncells-1,Ncells-1) ] = cchop
	
	return Ht

def make_Vbulk(Ncells, cchop, aCC):
	# Velocity operator
	Vx = np.zeros( ( 2*Ncells*Ncells , 2*Ncells*Ncells ), dtype=np.complex128 )
	Vxconst = 1j*cchop*aCC*np.sqrt(3)/2
	for n in np.arange(Ncells-1):
		for m in np.arange(Ncells-1):
			Vx[ ind(0,n,m) , ind(1,n+1,m) ] = Vxconst
			Vx[ ind(1,n+1,m) , ind(0,n,m) ] = -Vxconst
			Vx[ ind(0,n,m) , ind(1,n,m+1) ] = -Vxconst
			Vx[ ind(1,n,m+1) , ind(0,n,m) ] = Vxconst
	# "Lower right" boundary
	for n in np.arange(Ncells-1):
		Vx[ ind(0,n,Ncells-1) , ind(1,n+1,Ncells-1) ] = Vxconst
		Vx[ ind(1,n+1,Ncells-1) , ind(0,n,Ncells-1) ] = -Vxconst
	# "Upper right" boundary
	for m in np.arange(Ncells-1):
		Vx[ ind(0,Ncells-1,m) , ind(1,Ncells-1,m+1) ] = -Vxconst
		Vx[ ind(1,Ncells-1,m+1) , ind(0,Ncells-1,m) ] = Vxconst
	'''
	Vy = np.zeros( ( 2*Ncells*Ncells , 2*Ncells*Ncells ), dtype=np.complex128 )
	Vyconst = 1j*cchop*aCC
	halfVyconst = 0.5*Vyconst
	for n in np.arange(Ncells-1):
		for m in np.arange(Ncells-1):
			Vy[ ind(0,n,m) , ind(1,n,m) ] = -Vyconst
			Vy[ ind(1,n,m) , ind(0,n,m) ] = Vyconst
			Vy[ ind(0,n,m) , ind(1,n+1,m) ] = halfVyconst
			Vy[ ind(1,n+1,m) , ind(0,n,m) ] = -halfVyconst
			Vy[ ind(0,n,m) , ind(1,n,m+1) ] = halfVyconst
			Vy[ ind(1,n,m+1) , ind(0,n,m) ] = -halfVyconst
	# "Lower right" boundary
	for n in np.arange(Ncells-1):
		Vy[ ind(0,n,Ncells-1) , ind(1,n,Ncells-1) ] = -Vyconst
		Vy[ ind(1,n,Ncells-1) , ind(0,n,Ncells-1) ] = Vyconst
		Vy[ ind(0,n,Ncells-1) , ind(1,n+1,Ncells-1) ] = halfVyconst
		Vy[ ind(1,n+1,Ncells-1) , ind(0,n,Ncells-1) ] = -halfVyconst
	# "Upper right" boundary
	for m in np.arange(Ncells-1):
		Vy[ ind(0,Ncells-1,m) , ind(1,Ncells-1,m) ] = -Vyconst
		Vy[ ind(1,Ncells-1,m) , ind(0,Ncells-1,m) ] = Vyconst
		Vy[ ind(0,Ncells-1,m) , ind(1,Ncells-1,m+1) ] = halfVyconst
		Vy[ ind(1,Ncells-1,m+1) , ind(0,Ncells-1,m) ] = -halfVyconst
	# Corner
	Vy[ ind(0,Ncells-1,Ncells-1) , ind(1,Ncells-1,Ncells-1) ] = -Vyconst
	Vy[ ind(1,Ncells-1,Ncells-1) , ind(0,Ncells-1,Ncells-1) ] = Vyconst
	'''
	# Done
	# return [Vx, Vy]
	return Vx

def apply_Hpbc(Ht, u, v, Ncells, pbchop, aCC):
	# Adds PBCs to the bulk Hamiltonian, input as Ht
	for m in np.arange(Ncells-1):
		Ht[ ind(1,0,m) , ind(0,Ncells-1,m) ] = pbchop*np.exp(1j*2*np.pi*u)
		Ht[ ind(0,Ncells-1,m) , ind(1,0,m) ] = pbchop*np.exp(-1j*2*np.pi*u)
	for n in np.arange(Ncells-1):
		Ht[ ind(1,n,0) , ind(0,n,Ncells-1) ] = pbchop*np.exp(1j*2*np.pi*v)
		Ht[ ind(0,n,Ncells-1) , ind(1,n,0) ] = pbchop*np.exp(-1j*2*np.pi*v)
	# Done
	return Ht

def apply_Vpbc(Vx, u, v, Ncells, pbchop, aCC):
	# Adds PBCs to the bulk velocity operator, input as Vx
	Vxconst = 1j*pbchop*aCC*np.sqrt(3)/2
	for m in np.arange(Ncells-1):
		Vx[ ind(0,Ncells-1,m) , ind(1,0,m) ] = Vxconst*np.exp(-1j*2*np.pi*u)
		Vx[ ind(1,0,m) , ind(0,Ncells-1,m) ] = -Vxconst*np.exp(1j*2*np.pi*u)
	for n in np.arange(Ncells-1):
		Vx[ ind(0,n,Ncells-1) , ind(1,n,0) ] = -Vxconst*np.exp(-1j*2*np.pi*v)
		Vx[ ind(1,n,0) , ind(0,n,Ncells-1) ] = Vxconst*np.exp(1j*2*np.pi*v)
	'''
	Vy = np.zeros( ( 2*Ncells*Ncells , 2*Ncells*Ncells ), dtype=np.complex128 )
	halfVyconst = 1j*pbchop*aCC*0.5
	for m in np.arange(Ncells-1):
		Vy[ ind(0,Ncells-1,m) , ind(1,0,m) ] = halfVyconst*np.exp(-1j*2*np.pi*u)
		Vy[ ind(1,0,m) , ind(0,Ncells-1,m) ] = -halfVyconst*np.exp(1j*2*np.pi*u)
	for n in np.arange(Ncells-1):
		Vy[ ind(0,n,Ncells-1) , ind(1,n,0) ] = halfVyconst*np.exp(-1j*2*np.pi*v)
		Vy[ ind(1,n,0) , ind(0,n,Ncells-1) ] = -halfVyconst*np.exp(1j*2*np.pi*v)
	'''
	# Done
	return Vx # return [Vx, Vy]

# Make the bulk Hamiltonian and velocity operator
Hbulk = make_Hbulk(Ncells, cchop, aCC)
Vxbulk = make_Vbulk(Ncells, cchop, aCC)
print("Finished constructing bulk Hamiltonian and velocity")

# Calculate and save analytic bounds on the eigenvalues
pristine_energy_bound = np.abs(apply_Hpbc(Hbulk,0,0,Ncells,pbchop,aCC)).sum(axis=0).max() # Row sum over the Hamiltonian with PBCs if applicable
energy_bound = pristine_energy_bound + w # Add disorder to row sum
# Save the minimum and maximum energy over all k and all disorder configurations
with open(ev_path+"/e-min-max-w%s-n%s.pkl" % (w, Ntot), "wb") as outfile:
	pickle.dump([-energy_bound, energy_bound], outfile, pickle.HIGHEST_PROTOCOL)

n_done = 0
n_total_runs = len(u_list)*len(v_list)*nrep
for rep in range(nrep):
	# Introduce disorder
	site_disorder = np.random.uniform(-1,1,Natom)
	Ht = Hbulk + w*np.diag(site_disorder)
	# Sum over the k-point mesh for fixed disorder
	for u in u_list:
		for v in v_list:
			# Get PBC parts of the Hamiltonian and modify the bulk part with it
			Ht = apply_Hpbc(Ht, u, v, Ncells, pbchop, aCC)
			# Set types and precision
			e_vecs = np.zeros((Natom,Natom), dtype=np.complex128)
			e_vals = np.zeros((Natom), dtype=np.float64)
			# Compute eigenpairs and matrix elements
			t0 = time.time()
			e_vals, e_vecs = np.linalg.eigh(Ht)
			t = time.time()-t0
			ts = t % 60
			tm = int((t-ts)/60 % 60)
			th = int((t-tm*60-ts)/3600 % 60)
			print("Finished diagonalizing in time %d:%d:%.1f" % (th, tm, ts))
			# Set types and precision
			Vxnm = np.zeros((Natom,Natom), dtype=np.complex128)
			Vxnm2 = np.zeros((Natom,Natom), dtype=np.float64)
			# Get PBC parts of Vx and modify the bulk part with it
			# [Vx, Vy] = add_Vpbc(Vx, u, v, Ncells, pbchop, aCC)
			Vx = apply_Vpbc(Vxbulk, u, v, Ncells, pbchop, aCC)
			t0 = time.time()
			Vxnm = e_vecs.conjugate().transpose().dot( Vx.dot( e_vecs ) )
			Vxnm2 = np.real( np.multiply( Vxnm.conjugate() , Vxnm ) )
			t = time.time()-t0
			ts = t % 60
			tm = int((t-ts)/60 % 60)
			th = int((t-tm*60-ts)/3600 % 60)
			print("Finished calculating velocity matrix elements in time %d:%d:%.1f" % (th, tm, ts))
			# Save the data
			data = [e_vals, Vxnm2, Ntot]
			with open(ev_path+"/e-v-w%s-n%s-pbc%su%.2fv%.2f-rep%s.pkl" % (w, Ntot, pbchop, u, v, rep), "wb") as outfile:
				pickle.dump(data, outfile, pickle.HIGHEST_PROTOCOL)
			# Status report
			n_done += 1
			print("Finished %d of %d" % (n_done, n_total_runs ) )
			# Update maximum and minimum energies
