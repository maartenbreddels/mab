# -*- coding: utf-8 -*-
import os
import numpy
import sys
from numpy import *
from mab import logging

logger = logging.getLogger("mab.grid")

class HyperQuadTree(object):
	def __init__(self, binary, dirname, dimension, max_newpoints, max_level, initial_subdivides=2, basename="quadtree.%03d"):
		self.binary = binary
		self.dimension = dimension
		self.initial_subdivides = initial_subdivides
		self.dirname = dirname
		self.basename = basename
		self.max_newpoints = max_newpoints
		self.max_level = max_level
		
		self.points = []
		self.newpoints = []
		self.unknown_points = []
		self.values = {}
		if not os.path.exists(dirname):
			logger.info("making dir %s" % dirname)
			os.makedirs(dirname)
			
	def normalization(self, f, iteration=0):
		maxvalue = max([f(*point) for point in self.points])
		return exp(maxvalue)
		
	def lognormalization(self, f, iteration=0):
		maxvalue = max([f(*point) for point in self.points])
		return maxvalue
		
	def update(self, f, iteration=0):
		maxvalue = max([f(*point) for point in self.points])
		#print self.points
		#print maxvalue
		#dsa
		for point in self.points:
			self.values[point] = (f(*point))-maxvalue
		self.write_solved(iteration=iteration)
		self.rewrite_known(iteration=iteration)
		
	def eval(self, f, iteration=0):
		unknown = os.path.join(self.dirname, self.format_name(iteration)+".unknown")
		known =  os.path.join(self.dirname, self.format_name(iteration)+".known")
		lines = file(unknown).readlines()
		lines.extend(file(known).readlines())
		for line in lines:
			point = tuple([float(k) for k in line.split()])
			#print point
			f(*point[:-1])
		
	def write_solved(self, iteration=0):
		unknown = os.path.join(self.dirname, self.format_name(iteration)+".unknown")
		solved =  os.path.join(self.dirname, self.format_name(iteration)+".solved")
		logger.debug("writing file %s" % solved)
		lines = file(unknown).readlines()
		file_solved = open(solved, "w")
		assert len(lines) == len(self.newpoints), "%d %d"% (len(lines), len(self.newpoints))
		#values = []
		#for line in lines:
		#	point = tuple([float(k) for k in line.split()])
		#	values.append(f(*point[:-1]))
		#print lines
		for line in lines:
			point = tuple([float(k) for k in line.split()])[:-1]
			value = self.values[point] - max(self.values.values())
			#value = self.values[point]#+20
			#print value
			#print value
			#value = self.values[point]
			print >>file_solved, line.strip(), value
		file_solved.close()
		
	def rewrite_known(self, iteration=0):
		known = os.path.join(self.dirname, self.format_name(iteration)+".known")
		f = file(known)
		lines = f.readlines()
		f.close()
		
		assert len(lines) == len(self.points) - len(self.newpoints)
		                                                    
		f = file(known, "w")
		logger.debug("(re)writing file %s" % known)
		for line in lines:
			parts = line.split()
			point = tuple([float(k) for k in parts[:-2]]) # clip of value and volume
			value = self.values[point]#+20
			value = self.values[point] - max(self.values.values())
			print >>f, " ".join(parts[:-1]), value
		f.close()
		
	def optimize(self, iteration=0):
		#self.write_solved(f, iteration=iteration)
		inputname = os.path.join(self.dirname, self.format_name(iteration))
		outputname = os.path.join(self.dirname, self.format_name(iteration+1))
		rest = "-- 0 0 0 %d %d" % (self.max_newpoints[iteration], self.max_level)
		template = "%s -d %d --output=%s --input=%s --limit-level-difference=1 --optimize " + rest
		cmd = template % (self.binary, self.dimension, outputname, inputname)
		print cmd
		err = os.system(cmd)
		if err:
			sys.exit(err)
	
	def grid(self, size=100, iteration=0):
		#for p in self.unknown_points:
		inputname = os.path.join(self.dirname, self.format_name(iteration))
		cmd = "%s -d %d --grid --input=%s " % (self.binary, self.dimension, inputname)
		for i in range(self.dimension):
			cmd += " %d " % size
		print cmd
		if os.system(cmd) != 0:
			sys.exit(-1)
		
		data = file(inputname+".grid").read()
		gridlist = eval(data)
		grid = numpy.array(gridlist)
		import pdb;
		#pdb.set_trace()
		return log(grid)
		
			
	def read(self, iteration=0):
		def readdata(name):
			grid = []
			inputname = self.format_name(iteration)
			filename = os.path.join(self.dirname, inputname + name)
			lines = file(filename).readlines()
			for line in lines:
				values = [float(k) for k in line.split()]
				grid.append(tuple(values))
			return numpy.array(grid)
			
		grid = readdata(".unknown")
		points = grid[:,:-1]
		for point in points:
			p = tuple(point)
			assert len(p) == self.dimension
			#print p
			self.points.append(p)
			self.newpoints.append(p)
		grid = readdata(".known")
		if len(grid) > 0:
			points = grid[:,:-2]
			for point in points:
				p = tuple(point)
				assert len(p) == self.dimension
				self.points.append(p)
		
		
	def run(self, args, opts, scope):
		self.init()
	
	def format_name(self, iteration):
		return self.basename % iteration
	
	def init(self):
		output_filename = os.path.join(self.dirname, self.format_name(0))
		cmd = "%s -d %d -i -o %s -s %d " % (self.binary, self.dimension, output_filename, self.initial_subdivides)
		cmd += " -- "
		#for a, b in self.ranges_org:
		#	cmd += str(a) + " "
		#for a, b in self.ranges_org:
		#	cmd += str(b) + " "
		for i in range(self.dimension):
			cmd += "0 "
		for i in range(self.dimension):
			cmd += "1 "
		print cmd
		if os.system(cmd) != 0:
			return

	def iterate(self, iteration, optimize=True):
		if optimize:
			inputname = os.path.join(self.dirname, self.format_name(iteration))
			outputname = os.path.join(self.dirname, self.format_name(iteration+1))
			rest = "-- 0 0 0 %d %d" % (self.max_newpoints[iteration], self.max_level)
			template = "%s -d %d --output=%s --input=%s --limit-level-difference=1 --optimize " + rest
			cmd = template % (self.binary, self.dimension, outputname, inputname)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
		else:
			outputname_known = os.path.join(self.dirname, self.format_name(iteration+1))+ ".known"
			
			inputname_known = os.path.join(self.dirname, self.format_name(iteration))+ ".known"
			cmd = "cp %s %s" % (inputname_known, outputname_known)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
			
			inputname_solved = os.path.join(self.dirname, self.format_name(iteration))+ ".solved"
			cmd = "cat %s >> %s" % (inputname_solved, outputname_known)
			print cmd
			err = os.system(cmd)
			if err:
				sys.exit(err)
			
		
