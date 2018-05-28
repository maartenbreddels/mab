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
import mab.logging as logging_new
import mab.utils.iniscope
#import mab.utils.numpy
#import itertools
#import mab.parallelize
#import bz2

logger = logging.getLogger("utils.program")

class Dummy(object):
	pass

class Extra(object):
	pass

	def load(self, name):
		return mab.gd.configure.createobj_fromini(self.ini, name, self.scope)

from mab.utils.iniscope import IniScope

def program2(kaplot=False, **kwargs):
	from optparse import OptionParser
	if kaplot:
		parser = defaultparser()
	else:
		parser = OptionParser("")
		
	for name, value in kwargs.items():
		if (value is True) or (value is False):
			parser.add_option("--" + name,      dest=name, action="store_true", default=value)
			parser.add_option("--no-" +name,   dest=name, action="store_false")
		elif isinstance(value, int):
			parser.add_option("--" + name, default=value, type=int)
		else:
			parser.add_option("--%s" % name, default=value)
		
	logging.add_options(parser)
	parser.add_option("--list", action="store_true", default=False)
	parser.add_option("--iteration", default=0, type=int)
	parser.add_option("--aliases")
	parser.add_option("--options")
	parser.add_option("--configurations", default="")
	parser.add_option("--globals")
	parser.add_option("--include")
	parser.add_option("--modelpath", default=None)
	parser.add_option("--root")
	parser.add_option("--type", default=None)
	
	args = sys.argv[1:]
	if kaplot:
		opts, args = parseargs(parser=parser)
	else:
		opts, args = parser.parse_args(args)
	
	logging.configure(opts)
	logging_new.configure(opts)
	
	aliases = []
	configurations = []
	if opts.configurations:
		for conf in opts.configurations.split(";"):
			configurations.append(conf)
	if opts.aliases:
		#print "aliases:", opts.aliases
		for alias in opts.aliases.split(";"):
			name, value = alias.split("=")
			aliases.append((name, value))
	options = Dummy()
	if opts.options:
		#print "options:", opts.options
		for option in opts.options.split(";"):
			name, value = option.split("=")
			setattr(options, name, value)
	globals = [("iteration", repr(opts.iteration)), ("cores", repr(opts.cores))]
	if opts.modelpath:
		globals.append(("modelpath", repr(opts.modelpath)))
	if opts.globals:
		logger.info("globals: %s" % opts.globals)
		for global_value in opts.globals.split(";"):
			logger.debug(" global: %s" % global_value)
			#print "--->", global_value
			name, value = global_value.split("=")
			if name[0] == "\"":
				name = name[1:]
			if value[-1] == "\"":
				value = value[:-1]
			try:
				globals.append((name, (value)))
				#print "-->", `global_value`, ",", `name`, "=", `value`
			except:
				print "Error evaluating global '%s' (expr: %s)" % (name, value)
				raise
	
		
	include_files = []
	if opts.type:
		include_files.append(opts.type+".ini")
	if opts.include:
		include_files.extend(opts.include.split(";"))
		
	#print aliases
	if opts.root:
		filename = os.path.join(opts.root, "mab.ini")
	else:
		filename = "mab.ini"
		
	scope = IniScope(filename, aliases=aliases, options=options, globals=globals, configurations=configurations)
	scope.init()
	scope.readfiles(*include_files)
	#if opts.type:
	#	scope.readfiles(opts.type+".ini")
	#else:
	#	scope.readfiles()
	#print opts.globals
	scope["opts"] = opts
	scope["args"] = args 
	return scope

def program(objectname, dynamic=False, kaplot=False,  **kwargs):
	#logger2 = logging.getLogger("gd.schw.bla")
	from optparse import OptionParser
	if kaplot:
		parser = defaultparser()
	else:
		parser = OptionParser("")
	
	#parser = OptionParser("")
	parser.add_option("--modelpath")
	#if schwparams:
	#	parser.add_option("--schwsetname")
	#	parser.add_option("--schwmodelname", default="single")
	parser.add_option("--aliases")
	parser.add_option("--options")
	parser.add_option("--globals")
	parser.add_option("--type", default=None)
	if dynamic:
		parser.add_option("--%s" % objectname, default=None)
	else:
		parser.add_option("--%s" % objectname, default=objectname)
	
	parser.add_option("--help-commands",      dest="help_commands", action="store_true", default=False)
		
		
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
	args = sys.argv[1:]
	if kaplot:
		opts, args = parseargs(parser=parser)
	else:
		opts, args = parser.parse_args(args)
	
	logging.configure(opts)
	
	aliases = []
	if opts.aliases:
		print "aliases:", opts.aliases
		for alias in opts.aliases.split(";"):
			name, value = alias.split("=")
			aliases.append((name, value))
	options = Dummy()
	if opts.options:
		print "options:", opts.options
		for option in opts.options.split(";"):
			name, value = option.split("=")
			setattr(options, name, value)
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
	
	if dynamic and not opts.help_commands:
		objectname = "_".join([objectname, getattr(opts, objectname)])
	else:
		objectname = getattr(opts, objectname)
	filenames = []
	if opts.type:
		filenames.append(os.path.join(modelpath, opts.type+".ini"))
	if hasattr(options, "modelset") and hasattr(options, "modelname"):
		filename = os.path.join(modelpath, opts.type, options.modelset, options.modelname, opts.type+".ini")
		print filename
		filenames.append(filename)
	print filenames
	#filenames.append(os.path.join(modelpath, opts.name+".ini"))
	extra = Extra()
	extra.args = args
	extra.opts = opts
	scope = {"arguments":extra, "options":options}
	ini, scope = mab.gd.configure.loadini(modelpath, scope=scope, objectname=objectname, aliases=aliases, filenames=filenames, globals=globals)
	#ini, scope = mab.utils.iniscope.loadini(modelpath, scope=scope, objectname=objectname, aliases=aliases, filenames=filenames, globals=globals)
	
	if opts.help_commands:
		sections = ini.sections()
		print "Available commands:"
		for section in sections:
			if section.startswith("command"):
				print "\t", section[len("command_"):]
		sys.exit(0)
	
	extra.ini = ini
	extra.scope = scope
	extra.modelpath = modelpath
	extra.objectname = objectname
	object = scope[objectname]
	#if schwparams:
	#	extra.schwsetname = opts.schwsetname
	#	extra.schwmodelname = opts.schwmodelname

	#parametersweep = mab.gd.configure.createobj_fromini(ini, opts.parametersweep, scope)
	return object, extra

