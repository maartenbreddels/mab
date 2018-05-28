# -*- coding: utf-8 -*-
from scipy.optimize import fsolve, fmin
from numpy import *
import mab
import os

class CommandLinkBestFit(object):
	def __init__(self, parameterset, bestfit_dirname):
		self.parameterset = parameterset
		self.bestfit_dirname = bestfit_dirname
		
		
	def run(self, args, opts, scope):
		self.parameterset.load()
		map = self.parameterset.valuemap
		keys = map.keys()
		values = map.values()
		
		logps = []
		for parametervalue in self.parameterset.parameter_values:
			logps.append(self.parameterset.valuemap[tuple(parametervalue.values)])
		index = argmax(logps)
			
		#print keys[index]
		#print values[index]
		#print self.parameterset.parameter_values
		#print len(self.parameterset.parameter_values)
		#print index
		#print len(keys)
		assert len(keys) == len(self.parameterset.parameter_values), "wrong iteration: %d %d" % (len(keys), len(self.parameterset.parameter_values))
		bestfit = self.parameterset.parameter_values[index]
		fullpath = os.path.join(scope["model_set_path"], bestfit.name)
		fullpath = os.path.abspath(fullpath)
		#print fullpath
		cmd = "ln -s --force %s %s" % (fullpath, self.bestfit_dirname)
		print "command:", cmd
		if os.path.exists(self.bestfit_dirname):
			os.system("rm %s" % self.bestfit_dirname)
		os.system(cmd)
		#list = self.parameterset.parameter_range_list
		#globals = [(list[i].name, keys[index][i]) for i in range(len(list))]
		#print globals
		#includes = [os.path.join(
		#subscope = mab.utils.iniscope.SubScope(scope, globals=globals, includes=includes)
		#subscope.load()
		#print subscope.scope["model_id_path"]
		
		#print keys
		#import pdb
		#pdb.set_trace()

