class __ParameterSweepCondor(object):
	def __init__(self, modelpath, runtemplate, setname, paramformats, rootdir, condor_requirements, condor_head_template, condor_item_template, execute=True, postfix="", **extra_params):
		self.modelpath = modelpath
		self.runtemplate = runtemplate
		self.setname = setname
		self.paramformats = paramformats
		self.execute = execute
		self.extra_params = extra_params
		self.condor_params = {"rootdir":rootdir, "condor_requirements":condor_requirements, "schwsetname":setname, "modelpath":modelpath}
		self.condor_head_template = condor_head_template
		self.condor_item_template = condor_item_template
		self.postfix = postfix
		
		self.paramlist = self.createparamlist(self.paramformats[0], self.paramformats[1:], "", [], [])
		self.paramnames = [param_name for param_name, param_format, param_values in self.paramformats]
		#print self.paramlist
		#for p in self.paramlist:
		#	print p
		#	
		#print len(self.paramlist)
		
	def createparamlist(self, paramformat_head, paramformat_tail, name, indices, values):
		paramlist = []
		param_name, param_format, param_values = paramformat_head
		for i, value in enumerate(param_values):
			indices_new = list(indices) + [i]
			values_new = list(values) + [value]
			if name:
				name_new = name + "/" + (param_format % value)
			else:
				name_new = (param_format % value)
			if len(paramformat_tail) >= 1:
				paramlist.extend(self.createparamlist(paramformat_tail[0], paramformat_tail[1:], name_new, indices_new, values_new))
			else:
				dirname = os.path.join(self.modelpath, "schw", self.setname, name_new)
				paramlist.append((name_new, dirname, indices_new, values_new))
		return paramlist
			
	
	def _prepare_dirs_and_files(self, paramformat_head, paramformat_tail, params, name, condor_arguments_list, execute=True):
		param_name, param_format, param_values = paramformat_head
		#print len(param_values)
		#count += len(param_values)
		for value in param_values:
			params_new = dict(params)
			params_new[param_name] = value
			if name:
				name_new = name + "/" + (param_format % value)
			else:
				name_new = (param_format % value)
			if len(paramformat_tail) >= 1:
				self._prepare_dirs_and_files(paramformat_tail[0], paramformat_tail[1:], params_new, name_new, condor_arguments_list, execute=execute)
			else:
				dirname = os.path.join(self.modelpath, "schw", self.setname, name_new)
				#print dirname, self.modelpath, self.setname, name_new
				if not os.path.exists(dirname):
					logger.debug("creating directory %s" % dirname)
					if execute:
						os.makedirs(dirname)
				else:
					logger.debug("directory %s already exists" % dirname)
				
				dirname = os.path.join(self.modelpath, "schw", self.setname, name_new, "log")
				if not os.path.exists(dirname):
					logger.debug("creating directory %s" % dirname)
					if execute:
						os.makedirs(dirname)
				else:
					logger.debug("directory %s already exists" % dirname)
				
				filename = os.path.join(self.modelpath, "schw", self.setname, name_new, "galaxy.ini")
				logger.debug("writing parameters to filename: %s" % filename)
				if execute:
					f = open(filename, "w")
					print >>f, "[parameters]"
					for key, value in params_new.items():
						print >>f, "%s=%r" % (key, value)
					
				script = file(self.runtemplate).read()
				vars = {}
				vars["modelpath"] = os.path.abspath(self.modelpath)
				vars["schwsetname"] = self.setname
				vars["schwmodelname"] = name_new
				script = script % vars
				filename = os.path.join(self.modelpath, "schw", self.setname, name_new, "run.sh")
				logger.debug("writing runscript to filename: %s" % filename)
				if execute:
					file(filename, "w").write(script)
					os.chmod(filename, 0700)
					
				condor_arguments = dict(self.condor_params)
				condor_arguments["schwmodelname"] = name_new
				condor_arguments_list.append(condor_arguments)
				
				filename = os.path.join(self.modelpath, "schw", self.setname, name_new, "condor-single.batch")
				condortxt = ""
				txt = file(self.condor_head_template).read()
				condortxt += string.Template(txt).substitute(**self.condor_params)
				condortxt += "\n"
				
				txt = file(self.condor_item_template).read()
				condortxt += string.Template(txt).substitute(**condor_arguments)
				condortxt += "\n"
				if execute:
					file(filename, "w").write(condortxt)
				

			
		
		
		
	def prepare_dirs_and_files(self):
		execute = self.execute
		
		dirname = os.path.join(self.modelpath, "schw", self.setname)
		condor_arguments_list = []
		if not os.path.exists(dirname):
			logger.debug("creating directory %s" % dirname)
			os.makedirs(dirname)
		else:
			logger.debug("directory %s already exists" % dirname)
			
		condortxt = ""
		txt = file(self.condor_head_template).read()
		condortxt += string.Template(txt).substitute(**self.condor_params)
		condortxt += "\n"
		
		for name, dirname, indices, values in self.paramlist:
			#dirname = os.path.join(self.modelpath, "schw", self.setname, name_new)
			#print dirname, self.modelpath, self.setname, name_new
			if not os.path.exists(dirname):
				logger.debug("creating directory %s" % dirname)
				if execute:
					os.makedirs(dirname)
			else:
				logger.debug("directory %s already exists" % dirname)
			
			for subdirname in ["log", "intermediate", "results"]:
				dirname = os.path.join(self.modelpath, "schw", self.setname, name, subdirname)
				if not os.path.exists(dirname):
					logger.debug("creating directory %s" % dirname)
					if execute:
						os.makedirs(dirname)
				else:
					logger.debug("directory %s already exists" % dirname)
			
			filename = os.path.join(self.modelpath, "schw", self.setname, name, "galaxy.ini")
			logger.debug("writing parameters to filename: %s" % filename)
			if execute:
				f = open(filename, "w")
				print >>f, "[parameters]"
				#for key, value in params_new.items():
				for key, value in zip(self.paramnames, values):
					print >>f, "%s=%r" % (key, value)
				
			script = file(self.runtemplate).read()
			vars = {}
			vars["modelpath"] = os.path.abspath(self.modelpath)
			vars["schwsetname"] = self.setname
			vars["schwmodelname"] = name
			script = script % vars
			filename = os.path.join(self.modelpath, "schw", self.setname, name, "run.sh")
			logger.debug("writing runscript to filename: %s" % filename)
			if execute:
				file(filename, "w").write(script)
				os.chmod(filename, 0700)
				
			condor_arguments = dict(self.condor_params)
			condor_arguments["schwmodelname"] = name
			#condor_arguments_list.append(condor_arguments)
			
			filename = os.path.join(self.modelpath, "schw", self.setname, name, "condor-single.batch")
			condor_singletxt = ""
			txt = file(self.condor_head_template).read()
			condor_singletxt += string.Template(txt).substitute(**self.condor_params)
			condor_singletxt += "\n"
			
			txt = file(self.condor_item_template).read()
			condor_singletxt_item = string.Template(txt).substitute(**condor_arguments) + "\n"
			condor_singletxt += condor_singletxt_item
			#condor_singletxt += "\n"
			if execute:
				file(filename, "w").write(condor_singletxt)
			
			condortxt += condor_singletxt_item
		
			
		filename = os.path.join(self.modelpath, "schw", self.setname, "condor.batch")
		logger.info("writing condor batch (%s)" % filename)
		file(filename, "w").write(condortxt)
		
		
		total = 1
		for paramformat in self.paramformats:
			name = paramformat[0]
			values = paramformat[-1]
			logger.info("parameter: %s, [%e, ..., %e] (#=%d)" % (name, min(values), max(values), len(values)))
			total *= len(values)
		logger.info("%d models" % total)
		
		 
	def collect(self):
		shape = []
		for paramformat in self.paramformats:
			name = paramformat[0]
			values = paramformat[-1]
			shape.append(len(values))
		pgrid = numpy.zeros(tuple(shape))
		logpgrid = numpy.zeros(tuple(shape)) - 1e500
		import scipy.interpolate as interpolate
		#betas =
			
		for name, dirname, indices, values in self.paramlist:
			#dirname = os.path.join(self.modelpath, "schw", self.setname, name_new)
			#print dirname, self.modelpath, self.setname, name_new
			
			filename = os.path.join(dirname, "results/probability" +self.postfix +".txt")
			if os.path.exists(filename):
				p = eval(file(filename).read().strip())
				pgrid[indices[0], indices[1]] = p
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
			#print indices, p
			
			filename = os.path.join(dirname, "results/logprobability" +self.postfix +".txt")
			if os.path.exists(filename):
				logp = eval(file(filename).read().strip())
				logpgrid[indices[0], indices[1]] = logp
			else:
				logger.error("file '%s' missing (model not finished?)" % filename)
			#print indices, p
		pgrid2 = exp(logpgrid)
		#print pgrid.max(), pgrid.min()
		#print pgrid2.max(), pgrid2.min()
		#print pgrid
		pgrid /= pgrid.max()
		pgrid2 = exp(logpgrid-logpgrid.max())
		pgrid = pgrid2
		if 1: 
			Nbeta = 500
			betamin, betamax = -3, 1
			#nr = 700
			#rmax = 1.0
			#xmin=(1.556303+0.45)
			#xmax=(5.181438-0.75)
			#oldrs = arange(0, nr)/(nr-1.)*rmax + 0.03
			#rs = arange(0, nr)/(nr-1.)*rmax + 0.03
			#linbetagrid = zeros((len(rs),Nbeta))
			#oldrs = self.storage3d
			nr = 50
			xmin=0
			xmax=1.5 # * 2610.9660574412533 
			rs = (arange(nr) + 0.5 )/ (nr) * (xmax - xmin) + xmin
			betagrid = zeros((len(rs), Nbeta))
			for name, dirname, indices, values in self.paramlist[1:]:
				filename = os.path.join(dirname, "results", "solution_moments3d" +self.postfix +".npy")
				if os.path.exists(filename):
					solution_moments3d = numpy.load(filename)
					varvr = solution_moments3d[4]
					varvphi = solution_moments3d[5]
					varvtheta = solution_moments3d[6]
					betas = 1 - (varvphi + varvtheta)/(2*varvr)
					if 0:
						print betas
						print varvr
						print varvphi
						print varvtheta
					betaindices = ((betas-betamin)/(betamax-betamin) * Nbeta).astype(int)
					#print ".", ~isnan(varvtheta)
					good_betaindices = (betaindices >= 0) & (betaindices < Nbeta) & (~isnan(varvtheta))
					betaindices = betaindices[good_betaindices]
					p = pgrid[indices[0], indices[1]]
					#print good_betaindices,betaindices
					betagrid[good_betaindices,betaindices] += p
					#print len(solution_moments3d)
					#f = interpolate.interp1d(oldrs, betagrid[:,i])
				else:
					#logger.error("file '%s' missing (model not finished?)" % filename)
					pass
		
		#print pgrid.max(), pgrid.min()
		#print pgrid2.max(), pgrid2.min()
		#print pgrid2
		from kaplot import *
		import scipy.misc.pilutil
		im = pgrid2
		import scipy.ndimage
		im = scipy.ndimage.gaussian_filter(im, 1.0)
		if 0:
			im = scipy.misc.pilutil.toimage(im)
			width, height = 500, 500
			import Image
			im = im.resize((width, height), Image.BILINEAR)
			im = scipy.misc.pilutil.fromimage(im)
			#im = pgrid2
		#im = scipy.ndimage.gaussian_filter(im, 0.5)
		#im = scipy.misc.pilutil.imresize(pgrid2, (200,200))
		#box()
		vsplit(box)
		#indexedimage(transpose(im))
		#print log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])
		#print (self.paramformats[0][-1]), (self.paramformats[1][-1])
		probimage2d(im, 0, 1, log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])) #, colormap="whiteblack")
		select(1)
		#probimage2d(betagrid, 0, 1, log10(self.paramformats[0][-1]), log10(self.paramformats[1][-1])) #, colormap="whiteblack")
		betagrid = scipy.ndimage.gaussian_filter(betagrid, [2, 3])
		if 0:
			for i in range(betagrid.shape[1]):
				if sum(betagrid[:,i]) > 0:
					betagrid[:,i] /= betagrid[:,i].sum()
		if 1:
			for i in range(betagrid.shape[0]):
				if sum(betagrid[i]) > 0:
					betagrid[i] /= betagrid[i].sum()
		#indexedimage(transpose(betagrid), colormap="whiterainbow")
		betas = (arange(Nbeta) + 0.5 )/ (Nbeta) * (betamax- betamin) + betamin
		probimage2d(betagrid, 0, 1, rs, betas) #, colormap="whiteblack")
		draw()
		
		
