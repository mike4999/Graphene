import random as rand
import numpy as np

def hamiltonian(xcells, ycells, disorder=1.0):
	'''Generate a hamiltonian matrix using the Anderson model for graphene. The lattice
	is specified by a width and height in number of unit cells.'''
	width = 4*xcells
	height = 2*ycells
	n = width*height # Number of sites
	ham = np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			if i == j:
				ham[i,j] = rand.uniform(-disorder/2.0, disorder/2.0)
			else:
				i_col, j_col = i//height, j//height
				i_row, j_row = i%height, j%height
				ham[i,j] = 0
				if j_row == i_row:
					if (j_col+1)%width == i_col or j_col == (i_col+1)%width:
						ham[i,j] = 1
						continue
				if (j_row+1)%height == i_row:
					if i_col%4 == 1 and (j_col+1)%width == i_col:
						ham[i,j] = 1
						continue
					elif i_col%4 == 2 and j_col == (i_col+1)%width:
						ham[i,j] = 1
						continue
				if j_row == (i_row+1)%height:
					if i_col%4 == 0 and j_col == (i_col+1)%width:
						ham[i,j] = 1
						continue
					elif i_col%4 == 3 and (j_col+1)%width == i_col:
						ham[i,j] = 1
						continue
	return ham

def ipr(vecs):
	n = vecs.shape[0]
	numer = np.sum(np.sum(vecs**2, axis=0)**2)
	denom = n * np.sum(np.sum(vecs**4, axis=0))
	return numer/denom
