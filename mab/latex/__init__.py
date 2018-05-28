# -*- coding: utf-8 -*-
class Table(object):
	def begin(self, header, *names):
		print r"\begin{tabular}{%s}" % header
		print " & ".join(names), r"\\"
		
	def end(self):
		print r"\end{tabular}"
		
	def hline(self):
		print r"\hline"

		
	def row(self, *values):
		print " & ".join(values), r"\\"
		
		
class Table2(Table):
	def __init__(self, scopes):
		self.scopes = scopes
		
		
	def run(self, args, opts, scope):
		for scope in self.scopes:
			scope.init()
			scope.readfiles()
		"""
\begin{tabular}{l|r|r|r|r}
Name & $v_\text{sys}$ & $\sigma_v$ & $R_\text{max}$ & foreground model \\
		
		
		"""
		self.begin("l|r|r|r|r|r", "Name", r"$v_\text{sys}$", r"$\sigma_v$", r"$R_{e,\text{max,gius}}$",r"$R_{e,\text{max,walker}}$", "foreground model")
		self.hline()
		for scope in self.scopes:
			rmax = scope["command_rmax"]
			rmax.load()
			Rmax1 = rmax.Rmaxs[0]
			Rmax2 = rmax.Rmaxs[1]
			model = scope["model_foreground"]
			fgmodel = r"$\mu = %.1f$ km/s, $\sigma = %.1f$ km/s " % (model.mu, model.sigma)
			self.row(scope["object_name_nice"], "%.1f" % scope["mean_vhelio"], "%.1f" % scope["mean_vsigma"], "%.1f" % Rmax1, "%.1f" % Rmax2, fgmodel) 
		
		self.end()
		

class TableCounts(Table):
	def __init__(self, scopes):
		self.scopes = scopes
		
		
	def run(self, args, opts, scope):
		for scope in self.scopes:
			scope.init()
			scope.readfiles()
		if 0:
			self.begin("l|r|r|r|r|r", "Name", r"$v_\text{sys}$", r"$\sigma_v$", r"$R_{e,\text{max,gius}}$",r"$R_{e,\text{max,walker}}$", "foreground model")
			self.hline()
		for scope in self.scopes:
			#rmax = scope["command_rmax"]
			#rmax.load()
			#Rmax1 = rmax.Rmaxs[0]
			#Rmax2 = rmax.Rmaxs[1]
			gius = scope["observation_helio_gius"]
			Ngius = len(gius.load())
			walker = scope["observation_helio_walker"]
			Nwalker = len(walker.load())
			Nmember = scope["Nmember"]
			#print Ngius, Nwalker, Nmember
			self.row(scope["object_name_nice"], "%d" % Ngius, "%d" % Nwalker, "%d" % Nmember)
			#fgmodel = r"$\mu = %.1f$ km/s, $\sigma = %.1f$ km/s " % (model.mu, model.sigma)
			#self.row(scope["object_name_nice"], "%.1f" % scope["mean_vhelio"], "%.1f" % scope["mean_vsigma"], "%.1f" % Rmax1, "%.1f" % Rmax2, fgmodel) 
		
		self.end()
		
