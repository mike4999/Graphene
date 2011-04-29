import os

def save3d(f, xx, yy, values):
	'''Save a set of coordinates associated with some values into the format used by gnuplot. xx
	and yy should be meshes in the format provided by numpy.meshgrid. xx, yy, and all entries in
	values must be the same shape.'''
	if type(values) != tuple:
		values = (values,)
	fmt = (len(values)+2)*' %12g' + '\n'
	for i, xrow in enumerate(xx):
		x = xrow[i]
		for j, y in enumerate(yy[i]):
			v = [val[i,j] for val in values]
			f.write(fmt % tuple([x, y]+v))
		f.write('\n')
