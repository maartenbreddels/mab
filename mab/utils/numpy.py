from __future__ import absolute_import
import numpy
import mmap

def mmapzeros(shape, dtype=numpy.float64):
	"""Creates a numpy array similar to numpy.zeros, but with memory mapped data, so any forked process shares the data.
	
		:param shape: shape of the numpy array
		:param dtype: type of numpy array
	
	This can avoid moving large numpy arrays around. 
	Can be usefull in combination with mab.parallelize.parallelize::
		
		y = mmapzeros((10000, 10))
		indices = range(10000)
		@parallellize(cores=10)
		def f(i):
			y[i] = arange(10) * i
			
		f(indices)
		
	now y is filled with the data, if no mmapped array was used, the data will not be filled!
			
	"""
	sizeof_type = numpy.zeros(1, dtype=dtype).nbytes
	f = mmap.mmap(-1, sizeof_type * numpy.prod(shape))
	ar = numpy.frombuffer(f, dtype=dtype)
	ar = ar.reshape(shape)
	return ar 
	
def range_and_borders(N, xmin, xmax):
	x = (numpy.arange(0, N)+0.5) / (N) * (xmax-xmin) + xmin
	x_borders = (numpy.arange(0, N+1)) / (N+0.0) * (xmax-xmin) + xmin
	return x, x_borders
	