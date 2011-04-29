#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

// Generate a random number uniformly distributed between low and high.
double uniform(double low, double high) {
	double norm = rand()/(double)RAND_MAX;
	return low + (high-low)*norm;
}

// Generate a hamiltonian matrix using the Anderson model. Width must be a multiple of
// 4 and height must be a multiple of 2 (to match up with unit cells). Returns 0 on
// success, nonzero on failure.
gsl_matrix *hamiltonian_alloc(int xcells, int ycells, double disorder) {
	int width = 4*xcells, height = 2*ycells;
	int n = width*height;
	gsl_matrix *ham = gsl_matrix_alloc(n, n);
	int i, j, i_col, j_col, i_row, j_row;
	for(i = 0; i < n; i += 1) {
		for(j = 0; j < n; j += 1) {
			if(i == j) {
				gsl_matrix_set(ham, i, j, uniform(-disorder/2.0, disorder/2.0));
			} else {
				i_col = i/height;
				j_col = j/height;
				i_row = i%height;
				j_row = j%height;
				gsl_matrix_set(ham, i, j, 0);
				if(j_row == i_row) {
					if((j_col+1)%width == i_col || j_col == (i_col+1)%width) {
						gsl_matrix_set(ham, i, j, 1);
						continue;
					}
				}
				if((j_row+1)%height == i_row) {
					if(i_col%4 == 1 && (j_col+1)%width == i_col) {
						gsl_matrix_set(ham, i, j, 1);
						continue;
					} else if(i_col%4 == 2 && j_col == (i_col+1)%width) {
						gsl_matrix_set(ham, i, j, 1);
						continue;
					}
				}
				if(j_row == (i_row+1)%height) {
					if(i_col%4 == 0 && j_col == (i_col+1)%width) {
						gsl_matrix_set(ham, i, j, 1);
						continue;
					} else if(i_col%4 == 3 && (j_col+1)%width == i_col) {
						gsl_matrix_set(ham, i, j, 1);
						continue;
					}
				}
			}
		}
	}
	return ham;
}

// Given an input array, put the values into n bins (n determined by freq->size/values->size)
// whose frequencies are stored to freq and midpoints are stored to centers.
void histogram(gsl_vector *input, gsl_vector *freq, gsl_vector *centers) {
	int i, j, n = freq->size;
	double min = gsl_vector_min(input), max = gsl_vector_max(input);
	double low, high, val;
	double step = (max-min)/n;
	for(i = 0; i < n; i += 1) {
		gsl_vector_set(centers, i, min + i*step + step/2.0);
		gsl_vector_set(freq, i, 0);
	}
	for(i = 0; i < input->size; i += 1) {
		val = gsl_vector_get(input, i);
		for(j = 0; j < n; j += 1) {
			low = min + j*step;
			high = low + step;
			if(val >= low) {
				if(val < high || (j == n-1 && val <= high)) {
					gsl_vector_set(freq, j, gsl_vector_get(freq, j)+1);
				}
			}
		}
	}
}

// Normalize a vector so its sum is 1
void normalize(gsl_vector *v) {
	int i;
	double sum = 0;
	for(i = 0; i < v->size; i += 1) {
		sum += gsl_vector_get(v, i);
	}
	gsl_vector_scale(v, 1.0/sum);
}

int main(int argc, char **argv) {
	int i, j, n;
	gsl_matrix *ham, *vecs;
	gsl_vector *vals;
	
	if(argc <= 1) {
		fprintf(stderr, "Usage: %s [--states | --ipr | --spacing]\n", argv[0]);
		return 0;
	}
	
	srand(time(NULL));
	
	// Generate hamiltonian
	ham = hamiltonian_alloc(10, 10, 1.0);
	n = ham->size1;
	
	// Find eigenvalues and eigenvectors of hamiltonian
	vals = gsl_vector_alloc(n);
	vecs = gsl_matrix_alloc(n, n);
	{
		gsl_eigen_symmv_workspace *wk = gsl_eigen_symmv_alloc(ham->size1);
		gsl_eigen_symmv(ham, vals, vecs, wk);
		gsl_eigen_symmv_free(wk);
		gsl_eigen_symmv_sort(vals, vecs, GSL_EIGEN_SORT_VAL_ASC);
	}
	
	if(argc >= 2) {
		if(strcmp(argv[1], "--states") == 0) {
			// Generate a histogram of the energies
			gsl_vector *bins = gsl_vector_alloc(50);
			gsl_vector *freqs = gsl_vector_alloc(bins->size);
			histogram(vals, freqs, bins);
			normalize(freqs);
			for(i = 0; i < bins->size; i += 1) {
				printf("%12g %12g\n", gsl_vector_get(bins, i), gsl_vector_get(freqs, i));
			}
			
		} else if(strcmp(argv[1], "--ipr") == 0) {
			// Find the IPR
			double numer = 0, numer_col;
			double denom = 0;
			double x;
			for(i = 0; i < n; i += 1) {
				numer_col = 0;
				for(j = 0; j < n; j += 1) {
					x = gsl_matrix_get(vecs, j, i);
					numer_col += pow(x, 2.0);
					denom += pow(x, 4.0);
				}
				numer += pow(numer_col, 2.0);
			}
			denom *= n;
			printf("%12g\n", numer/denom);
			
		} else if(strcmp(argv[1], "--spacing") == 0) {
			// Find the spacing
			gsl_vector *diffs = gsl_vector_alloc(n-1);
			for(i = 0; i < n-1; i += 1) {
				double a = gsl_vector_get(vals, i);
				double b = gsl_vector_get(vals, i+1);
				gsl_vector_set(diffs, i, b-a);
			}
			
			// Generate a histogram
			gsl_vector *bins = gsl_vector_alloc(100);
			gsl_vector *freqs = gsl_vector_alloc(bins->size);
			histogram(diffs, freqs, bins);
			normalize(freqs);
			for(i = 0; i < bins->size; i += 1) {
				printf("%12g %12g\n", gsl_vector_get(bins,i), gsl_vector_get(freqs,i));
			}
			
		} else {
			fprintf(stderr, "Unrecognized argument\n");
			return -1;
		}
	}
	
	return 0;
}
