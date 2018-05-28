import numpy

def bingrid(borders, x, *rest):
	i = 0
	while i < len(borders)-1:
		xmin, xmax = borders[i], borders[i+1]
		if i == len(borders) - 2:
			indices = (x >= xmin) & (x <= xmax)
		else:
			indices = (x >= xmin) & (x < xmax)
		#values = tuple([x[indices[i1:i2]] + [k[indices[i1:i2]] for k in rest])
		#if sum(indices) > 0:
		#	print mean(x[indices]), min(x[indices]), max(x[indices]), xmin, xmax
		#else:
		#	print "-------", xmin, xmax 
			 
		values = tuple([x[indices]] + [k[indices] for k in rest])
		yield values
		i += 1

def binrange(N, x, *rest):
	i1 = 0
	#if sort:
	indices = numpy.argsort(x)
	while i1 < len(x):
		i2 = min(len(x), i1+N)
		#values = tuple([x[indices[i1:i2]] + [k[indices[i1:i2]] for k in rest])
		values = tuple([x[indices[i1:i2]]] + [k[indices[i1:i2]] for k in rest])
		yield (i2-i1,) + values
		i1 += N