# -*- coding: utf-8 -*-
from optparse import OptionParser
import sys
import numpy
import os
from numpy import *
from mab.random import *
from kaplot import *
from scipy.optimize import fsolve, fmin
import mab.gd.schw.grid
import mab.utils.progressbar
import mab.gd.logging as logging
import mab.gd.configure
import mab.gd.schw.configure
import mab.gd.jeans.configure
import mab.gd.configure
import mab.gd.logging as logging
#import mab.utils.numpy
#import itertools
#import mab.parallelize
#import bz2

class Extra(object):
	pass

	def load(self, name):
		return mab.gd.configure.createobj_fromini(self.ini, name, self.scope)

def program(objectname, kaplot=False, schw=False, schwparams=False, jeans=False, **kwargs):
	#logger = logging.getLogger("gd.schw.ic")
	#logger2 = logging.getLogger("gd.schw.bla")
	from optparse import OptionParser
	if kaplot:
		parser = defaultparser()
	else:
		parser = OptionParser("")
	
	#parser = OptionParser("")
	parser.add_option("--modelpath")
	parser.add_option("--%s" % objectname, default=objectname)
	if schwparams:
		parser.add_option("--schwsetname")
		parser.add_option("--schwmodelname", default="single")
	parser.add_option("--aliases")
	parser.add_option("--globals")
	parser.add_option("--jeans", default=jeans)
		
	#parser.add_option("--cores", default=1, type=int)
	#parser.add_option("--%s" % objectname, default=objectname)
	for name, value in kwargs.items():
		if (value is True) or (value is False):
			parser.add_option("--" + name,      dest=name, action="store_true", default=value)
			parser.add_option("--no-" +name,   dest=name, action="store_false")
		elif isinstance(value, int):
			parser.add_option("--" + name, default=value, type=int)
		else:
			parser.add_option("--%s" % name, default=value)
		
	logging.add_options(parser)
	print "lala"
	args = sys.argv[1:]
	if kaplot:
		opts, args = parseargs(parser=parser)
	else:
		opts, args = parser.parse_args(args)
	logging.configure(opts)
	jeans = bool(opts.jeans)
	
	aliases = []
	if opts.aliases:
		print "aliases:", opts.aliases
		for alias in opts.aliases.split(";"):
			name, value = alias.split("=")
			aliases.append((name, value))
	globals = []
	if opts.globals:
		print "globals:", opts.globals
		for global_value in opts.globals.split(";"):
			name, value = global_value.split("=")
			globals.append((name, value))
	
	modelpath = opts.modelpath
	if modelpath is None:
		#modelpath = file("default_galaxy_modelpath").read().strip()
		lines = file("default_galaxy_modelpath").readlines()
		lines = [k.strip() for k in lines]
		lines = [k for k in lines if not k.startswith("#")]
		assert len(lines) >= 1, "at least 1 non-commented line required"
		modelpath = lines[0]
	print "modelpath:", modelpath
	
	objectname = getattr(opts, objectname)
	if schw and not jeans:
		if schwparams:
			ini, scope = mab.gd.schw.configure.loadini(modelpath, opts.schwsetname, opts.schwmodelname, objectname=objectname, aliases=aliases, globals=globals)
		else:
			ini, scope = mab.gd.schw.configure.loadini(modelpath, objectname=objectname, aliases=aliases, globals=globals)
	elif jeans:
		ini, scope = mab.gd.jeans.configure.loadini(modelpath, objectname=objectname, aliases=aliases, globals=globals)
	else:
		ini, scope = mab.gd.configure.loadini(modelpath, objectname=objectname, aliases=aliases, globals=globals)

	object = scope[objectname]
	extra = Extra()
	extra.args = args
	extra.opts = opts
	extra.ini = ini
	extra.scope = scope
	extra.modelpath = modelpath
	extra.objectname = objectname
	if schwparams:
		extra.schwsetname = opts.schwsetname
		extra.schwmodelname = opts.schwmodelname

	#parametersweep = mab.gd.configure.createobj_fromini(ini, opts.parametersweep, scope)
	return object, extra

