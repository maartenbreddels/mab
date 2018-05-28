from .. import configure
import os

import mab.gd.logging as logging
logger = logging.getLogger("jeans.configure")

def loadini(modelpath, scope=None, objectname="jeansmodel", aliases=[], globals=[]):
	#if schwsetname:
	#	scope["schwsetname"] = schwsetname
	#if schwmodelname:
	#	scope["schwmodelname"] = schwmodelname
	#ini, scope = configure.loadini(modelpath, scope, objectname=objectname)
	filenames = [os.path.join(modelpath, "jeans.ini")]
	#if schwsetname is not None:
	#	filenames.append(os.path.join(modelpath, "schw", schwsetname, "schwmodel.ini"))
	#if schwmodelname is not None:
	#	filenames.append(os.path.join(modelpath, "schw", schwsetname, schwmodelname, "galaxy.ini"))
	#if schwmodelname is not None:
	#	logger.info("schw path: %s" % (os.path.join(modelpath, "schw", schwsetname)))
	#ini.read(filenames)
	#for filename in filenames:
	#	print filename, "exists?: ", os.path.exists(filename)
	ini, scope = configure.loadini(modelpath, scope=None, objectname=objectname, filenames=filenames, aliases=aliases, globals=globals)
	return ini, scope