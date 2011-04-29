import numpy as np
import math as m

def wavevector_energy(kx, ky, valence=False):
	'''Determine the energy of a graphene lattice for each wave vector in the given
	array. The x coordinate of the wave vectors is given by kx and the y coordinate
	is given by ky. The shape of kx and ky must match. If `valence` is true, this
	gives the energy of the valence band instead of the conductive band.'''
	
	# Get wave vector array in unit cell vector basis
	a = np.matrix([
		[m.cos(m.pi/6.0), m.sin(m.pi/6.0)],
		[m.cos(m.pi/6.0), -m.sin(m.pi/6.0)]
	])
	k1 = (a[0,0]*kx+a[0,1]*ky) / (2.0*np.pi)
	k2 = (a[1,0]*kx+a[1,1]*ky) / (2.0*np.pi)
	
	# Helper functions
	u = lambda k1, k2: 2*(np.cos(2*np.pi*k1) + np.cos(2*np.pi*k2) + np.cos(2*np.pi*(k1-k2)))
	f = lambda k1, k2: 3 + u(k1, k2)
	
	# Compute energies
	e = np.sqrt(f(k1, k2))
	if valence:
		e *= -1
	return np.array(e)
