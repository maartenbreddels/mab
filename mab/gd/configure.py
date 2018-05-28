import ConfigParser as configparser
import sys
import os
import numpy
from math import *
import mab.gd.logging as logging
import math
import collections
import mab.ordereddict
import glob

logger = logging.getLogger("gd.configure")

enable_logging = True
enable_imports = True

module_exts = [".py"]
imported = []

def recursive_import(basedirname, dirname):
	if not enable_imports:
		return
	if dirname in imported:
		return
	imported.append(dirname)
	filenames = glob.glob(os.path.join(dirname, "*"))
	for filename in filenames:
		if  os.path.isdir(filename):
			#print "DIR", filename
			recursive_import(basedirname, filename)
		else:
			name, ext = os.path.splitext(filename)
			if ext in module_exts:
				modulename = name[len(basedirname)+1:].replace("/", ".")
				#print "MOD", modulename
				if enable_logging:
					logger.debug("basedir: %s; name: %s; modulename: %s" % (basedirname, name, modulename))
				exec "import %s" % modulename
				
				

class AutoloadScope(object):
	def __init__(self, filenames, aliases, globals, **defaults):
		self.filenames = filenames
		self.aliases = aliases
		self.globals = globals
		self.defaults = defaults
		self.reset()
		
	def __setitem__(self, name, value):
		#print "set:", name
		self.dict[name] = value
	def __getitem__(self,  name):
		#print "get:", name
		if name not in self.dict:
			if name in self.ini.sections():
				self.dict[name] = createobj_fromini(self.ini, name, self)
		return self.dict[name]
	
	def __contains__(self, name):
		return self.dict.__contains__(name)
	
	def flush(self):
		self.dict = dict()
	
	def __delitem__(self, name):
		del self.dict[name]
		
	def reset(self):
		self.dict = dict(self.defaults)
		self.ini = configparser.ConfigParser(dict_type=mab.ordereddict.OrderedDict)
		self.ini.optionxform = str
		
	def readfiles(self, *extrafilenames):
		allfilenames = self.filenames + list(extrafilenames)
		#print "all filenames", allfilenames 
		for filename in allfilenames:
			self.ini.read(filename)
		
	def init(self):
		ini = self.ini
		aliases = self.aliases
		if "imports" in ini.sections():
			for name, value in ini.items("imports"):
				logger.debug("import: %s (%s)" % (name, value))
				exec "import %s" % name
				basename = name.split(".")[0]
				#print "   basename", basename
				basemodule = eval(basename)
				module = eval(name)
				self[basename] = basemodule
				if value == "recursive":
					dirname = os.path.dirname(module.__file__)
					#print module, dirname
					basedirname = os.path.split(os.path.dirname(basemodule.__file__))[0]
					recursive_import(basedirname, dirname)
					#print dirname
					#sys.exit(0)
		if "globals" in ini.sections():
			self["configuration"] = self
			for name, value in ini.items("globals"):
				try:
					v = eval(value, globals(), self)
				except:
					logger.error("error evaluating global %s=%s" % (name, value))
					raise
				logger.debug("global: %s=%r (expr: %s)" % (name, v, value))
				self[name] = v
		if "parameters" in ini.sections():
			for name, value in ini.items("parameters"):
				v = eval(value, globals(), self)
				logger.debug("parameter: %s=%r (expr: %s)" % (name, v, value))
				self[name] = v
				#print "param:", name, value
		#print aliases
		for name, value in aliases:
			ini.set(name, "alias", value)
			logger.debug("setting alias: %s=%s " % (name, value))
		for name, value in self.globals:
			#self[name] = value
			logger.debug("(overriding) setting global: %s=%s " % (name, value))
			#print name, value
			self[name] = eval(value, globals(), self)


def createobj_fromini(ini, name, scope):
	#def get(name):
	#	print "$$$$$$$$$$$$", name
	#	return createobj_fromini(ini, name, scope)
	if enable_logging:
		logger.debug("create: %s" % name)
	items = ini.items(name)
	kwargs = {}
	args = None
	classname = None
	constructor = None
	#print len(items), "alias" in dict(items), dict(items)
	#scope["get"] = get
	if "alias" in dict(items):
		value = dict(items)["alias"]
		obj = eval(value, globals(), scope)
		if enable_logging:
			logger.debug("alias : %s -> %r (expr: %s)" % (name, obj, value))
		#obj = createobj_fromini(ini, dict(items)["alias"], scope)
		#scope["get"] = get
	elif "alias_condition" in dict(items):
		if enable_logging:
			logger.debug("alias (true) %s -> %s" % (name, dict(items)["alias_true"]))
			logger.debug("alias (false) %s -> %s" % (name, dict(items)["alias_false"]))
		expr = dict(items)["alias_condition"]
		result = eval(expr, globals(), scope)
		if enable_logging:
			logger.debug("condition: %s == %r" % (expr, result))
		if result:
			obj = createobj_fromini(ini, dict(items)["alias_true"], scope)
		else:
			obj = createobj_fromini(ini, dict(items)["alias_false"], scope)
		#scope["get"] = get
	else:
		for arg, value in items:
			#print name, ":", arg
			if arg == "class":
				classname = value
			elif arg == "constructor":
				constructor = value
			elif arg == "args":
				args = []
				items = value.split(",")
				for item in items:
					if item in scope:
						item = scope[item]
					elif ini.has_section(item):
						item = createobj_fromini(ini, item, scope)
						#scope["get"] = get
					else:
						try:
							item = eval(item, globals(), scope)
						except:
							if enable_logging:
								logger.error("error evaluating key %s in section %s (value is %r)" % (arg, name, item))
							raise
					args.append(item)
			else:
				if value in scope:
					value = scope[value]
				elif ini.has_section(value):
					#print "create obj"
					value = createobj_fromini(ini, value, scope)
					#scope["get"] = get
				elif value == "parameter":
					nestedname = "%s.%s" % (name, arg)
					if nestedname not in scope:
						raise Exception, "missing parameter: %s" % nestedname
					value = scope[nestedname]
				else:
					try:
						value = eval(value, globals(), scope)
					except:
						logger.error("error evaluating key %s in section %s (value is %r)" % (arg, name, value))
						raise
				kwargs[arg] = value
		assert classname is not None, "class attribute not given for section: %s" % name
		fullclassname = classname
		if 0:
			modulename, classname = fullclassname.rsplit(".", 1)
			mod = __import__(modulename)
			obj = mod
			for attrname in fullclassname.split(".")[1:]:
				obj = getattr(obj, attrname)
			#print obj
		try:
			obj = eval(fullclassname, globals(), scope)
		except:
			logger.error("error evaluating classname: %r" % fullclassname)
			raise
		try:
			if constructor:
				args = []
				for param in constructor.split(","):
					args.append(kwargs[param])
				obj = obj(*args)
			else:
				if args:
					obj = obj(*args, **kwargs)
				else:
					obj = obj(**kwargs)
		except:
			logger.error("error constructing object of class: %s (section %s)" % (fullclassname, name))
			raise
		#print name
	scope[name] = obj
	return obj

if 0:
	for section in ini.sections():
		print "[%s]" % section
		for name, value in ini.items(section):
			print "%s=%s" % (name, value)
			

#$scope = Scope()
def creategalaxy():
	parameters = {}
	for name, value in ini.items("parameters"):
		parameters[name] = eval(value)
	galaxy = createobj_fromini(ini, "galaxy", parameters)
	return galaxy
	
def loadini(modelpath, scope=None, objectname="galaxy", filenames=[], aliases={}, globals=[]):
	global enable_imports
	#assert scope is None
	if scope is None:
		scope = dict()
	#filename = sys.argv[1]
	scope["modelpath"] = modelpath
	#filenames = []
	
	filename = os.path.join(modelpath, "galaxymodel.ini")
	filenames.insert(0, filename)
	filename = os.path.join(modelpath, "galaxy.ini")
	filenames.insert(0, filename)
	
	allfilenames = []
	for filename in filenames:
		ini = configparser.ConfigParser()
		ini.optionxform = str
		if os.path.exists(filename):
			ini.read(filename)
			if "options" in ini.sections():
				values = dict(ini.items("options"))
				if "include" in values:
					dirname = os.path.dirname(filename)
					f = os.path.abspath(os.path.join(dirname, values["include"])) 
					allfilenames.append(f)
					#print values["include"], dirname
					assert os.path.exists(f), "filename %s doesn't exist" % f
			allfilenames.append(filename)
	
	logger.debug("used ini files: %r" % allfilenames)
	#assert os.path.exists(filename), "%s not found" % filename 
	if isinstance(scope, AutoloadScope):
		scope = scope.dict 
	scope = AutoloadScope(allfilenames, aliases, globals, **scope)
	scope.readfiles()
	ini = scope.ini
	#assert os.path.exists(filename), "%s not found" % filename
	#if os.path.exists(filename):
	#	ini.read(filename)
	#for filename in filenames:
	#	ini.read(filename)
	scope.init()
	#createobj_fromini(ini, objectname, scope)
	#print scope.keys()
	#enable_imports = False
	return ini, scope
	
#print creategalaxy()