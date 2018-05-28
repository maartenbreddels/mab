import mab.gd.logging as logging

logger = logging.getLogger("gd.commands")



class SubMean(object):
	def __init__(self, input_observation, output_observation, vmean=None):
		self.input_observation = input_observation
		self.output_observation = output_observation
		self.vmean = vmean
	
	def run(self, args, opts, scope):
		stars = self.stars = self.input_observation.load()
		if self.vmean is None:
			vmean = stars.vlos_helio.mean()
		else:
			vmean = self.vmean
		logger.info("subtracting mean velocity: %f" % vmean)
		stars.vlos = stars.vlos_helio - vmean
		stars[0].attributes.append("vlos")
		self.output_observation.save(self.stars)

class Aperture(object):
	def __init__(self, aperture):
		self.aperture = aperture
		
	def run(self, args, opts, scope):
		self.aperture.init()
		self.aperture.process()
		self.aperture.save()
		
class BinData(object):
	def __init__(self, binned_data):
		self.binned_data = binned_data
		
	def run(self, args, opts, scope):
		self.binned_data.prepare()
		self.binned_data.store()		
		