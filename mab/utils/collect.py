# -*- coding: utf-8 -*-
import os
from numpy import *


class dSph(object):
	def __init__(self, objects, output):
		self.objects = objects
		self.output = output
	
	def run(self, args, opts, scope):
		evidences_list = []
		for object in self.objects:
			modelpath = "models/%s" % object
			evidences = []
			#for name in "nfw core_g1".split(): # core_g2 einasto_a0.2 einasto_a0.4".split():
			for name in "nfw einasto_a0.2 einasto_a0.4 core_g1_3 core_g2_3 core_g1_4 core_g2_4".split():
				filename = os.path.join(modelpath, "schw", name, "evidence.txt")
				value = eval(file(filename).read())
				evidences.append(value)
			evidences_list.append(evidences)
			
		combined = numpy.prod(evidences_list, axis=0)
		if 1:
			evidences_list.append(combined)
		evidences_list = numpy.array(evidences_list).T
		evidences_list = evidences_list * 1/evidences_list[0,:]
		evidences_list = evidences_list.tolist()
		print evidences_list
			#globals = "
