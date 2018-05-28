from .. import configure
import os

import mab.gd.logging as logging
logger = logging.getLogger("gd.configure")

def loadini(modelpath, schwsetname=None, schwmodelname=None, scope=None, objectname="schwmodel", aliases=[], globals=None):
	if globals is None:
		globals = []
	if schwsetname:
		#scope["schwsetname"] = schwsetname
		globals.append(("schwsetname", "'"+schwsetname+"'"))
	if schwmodelname:
		#scope["schwmodelname"] = schwmodelname
		globals.append(("schwmodelname", "'"+schwmodelname+"'"))
	#ini, scope = configure.loadini(modelpath, scope, objectname=objectname)
	filenames = [os.path.join(modelpath, "schwmodel.ini")]
	if schwsetname is not None:
		filenames.append(os.path.join(modelpath, "schw", schwsetname, "schwmodel.ini"))
	if schwmodelname is not None:
		filenames.append(os.path.join(modelpath, "schw", schwsetname, schwmodelname, "galaxy.ini"))
	if schwmodelname is not None:
		logger.info("schw path: %s" % (os.path.join(modelpath, "schw", schwsetname)))
	#ini.read(filenames)
	for filename in filenames:
		pass #print filename, "exists?: ", os.path.exists(filename), `filename`
		#os.system("ls %s" % filename)
	ini, scope = configure.loadini(modelpath, scope=None, objectname=objectname, filenames=filenames,aliases=aliases, globals=globals)
	return ini, scope
