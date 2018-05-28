# -*- coding: utf-8 -*-
from numpy import *
#import numpy

def markmax(image, x, y, id=1):
	markers = (image * 0).astype(int8)
	markers[x,y] = id
	
	Nx, Ny = markers.shape
	#maxiters = 1000
	#iteration = 0

	moved = True
	while moved:
		moved = False
		for i in range(Nx)[::-1]:
			#for j in range(Ny):
			for j in range(Ny)[::-1]:
				value = image[i,j]
				marker = markers[i,j]
				#print i, j, marker, value
				if markers[i,j] > 0:
					if i > 0: # left
						if (markers[i-1,j] == 0) and (image[i-1,j] > value):
							markers[i-1,j] = marker
							markers[i,j] = 0
							#print "left"
							moved = True
					if i < Nx-1: # right
						if (markers[i+1,j] == 0) and (image[i+1,j] > value):
							markers[i+1,j] = marker
							markers[i,j] = 0
							#print "right"
							moved = True
					if j > 0: # above
						if (markers[i,j-1] == 0) and (image[i,j-1] > value):
							markers[i,j-1] = marker
							markers[i,j] = 0
							#print "above"
							moved = True
					if j < Ny -1: # below
						if (markers[i,j+1] == 0) and (image[i,j+1] > value):
							markers[i,j+1] = marker
							markers[i,j] = 0
							#print "below"
							moved = True
	return markers

def watershed(image, markers):
	moved = True
	Nx, Ny = markers.shape
	while (sum(markers == 0) > 0) and moved: # if contains zeros.. continue
		moved = False
		#print "#" * 70
		#print markers
		for i in range(Nx)[::-1]:
			for j in range(Ny)[::-1]:
				value = image[i,j]
				marker = markers[i,j]
				#print i, j, marker, value
				if markers[i,j] > 0:
					if i > 0: # left
						if (markers[i-1,j] == 0) and (image[i-1,j] < value):
							markers[i-1,j] = marker
							#print "left"
							moved = True
						if j > 0: # left+above
							if (markers[i-1,j-1] == 0) and (image[i-1,j-1] < value):
								markers[i-1,j-1] = marker
								#print "left+above"
								moved = True
						if j < Ny -1: # left+below
							if (markers[i-1,j+1] == 0) and (image[i-1,j+1] < value):
								markers[i-1,j+1] = marker
								#print "left+below"
								moved = True
					if i < Nx-1: # right
						if (markers[i+1,j] == 0) and (image[i+1,j] < value):
							markers[i+1,j] = marker
							#print "right"
							moved = True
						if j > 0: # right+above
							if (markers[i+1,j-1] == 0) and (image[i+1,j-1] < value):
								markers[i+1,j-1] = marker
								#print "right+above"
								moved = True
						if j < Ny -1: # right+below
							if (markers[i+1,j+1] == 0) and (image[i+1,j+1] < value):
								markers[i+1,j+1] = marker
								#print "right+below"
								moved = True
					if j > 0: # above
						if (markers[i,j-1] == 0) and (image[i,j-1] < value):
							markers[i,j-1] = marker
							#print "above"
							moved = True
					if j < Ny -1: # below
						if (markers[i,j+1] == 0) and (image[i,j+1] < value):
							markers[i,j+1] = marker
							#print "below"
							moved = True
