import numpy


def optimize(f, xbest0, Nmax=500, vertices=None):
	dimx = len(xbest0)
	#y = f(xbest0)
	numpy.random.seed(0)
	if vertices is None:
		vertices = numpy.zeros((dimx+1, dimx))
		if 1:
			for i in range(dimx+1):
				vertices[i] = numpy.random.random(dimx)
				vertices[i] /= sum(vertices[i])
		else:
			for i in range(1,dimx+1):
				vertices[i][i-1] = 1
			vertices[i] = 1
			
		
	
	fy = numpy.array([f(x) for x in vertices])
	
	alpha = 1.
	gamma = 2.
	rho = 0.5
	sigma = 0.5
	
	done = False
	n = 0
	while not done:
		indices = numpy.argsort(fy)
		vertices = vertices[indices]
		#print `indices`, `vertices`, `fy`
		fy = fy[indices]
		
		# calculate center of gravity
		x0 = sum(vertices[:-1], 0)/dimx
		
		# reflected point
		xr = x0 + alpha * (x0 - vertices[-1])
		fr = f(xr)
		#if n == 0:
		#print "-" * 20, n
		#print vertices
		#print fy
		#print x0, xr, fr
		
		if (fr < fy[-2]): # better than worst
			#if debug:
			#	print "reflected point better than worst"
			if (fr >= fy[0]): # not better than the best
				# ok, simply copy vertex
				vertices[-1] = xr
				fy[-1] = fr
				print ".",
				#print xr, fr
			else: # best point so far
				# expanded point
				#print "*",
				xe = x0 + gamma * (x0 - vertices[-1])
				mask = xe < 0
				fe = f(xe)
				if fe < fr:
					# expanded better than reflected, copy expanded
					print "e",
					vertices[-1] = xe
					fy[-1] = fe
				else: # else copy reflected
					print "r",
					vertices[-1] = xr
					fy[-1] = fr
		else: # f(xr) <= f(worst)
			# compute contracted point
			xc = vertices[-1] + rho * (x0 - vertices[-1])
			fc = f(xc)
			fworst = fy[-2]
			if fc <= fworst: # contracted better than worst
				# copy contracted
				print "c", 
				#print xc, fc
				vertices[-1] = xc
				fy[-1] = fc
			else: # reduction
				print "red", 
				#print xc, fc
				for i in range(1,dimx+1): # for all but the best point
					vertices[i] = vertices[0] + sigma * (vertices[i] - vertices[0])
					#mask = vertices[i] < 0
					#vertices[i][mask] = xbest[mask] + (random.random(dim)[mask]-0.5) * 1e-6
					#vertices[i] = maximum(0, vertices[i])
					#vertices[i] /= sum(vertices[i])
					fy[i] = f(vertices[i])
					# - log10(random.random()) * Temp
		#if 0 and n == 50:
		#	i = numpy.argmin(fy)
		#	best = vertices[i]
		#	for i in range(dimx+1):
		#		vertices[i] = best + numpy.random.random(dimx) * 0.1-0.1/2
			
			
		n += 1
		if n == Nmax:
			done = True
	i = numpy.argmin(fy)
	return vertices[i]

