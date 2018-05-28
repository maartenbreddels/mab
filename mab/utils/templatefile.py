import mab.gd.logging as logging

logger = logging.getLogger("gd.utils.templatefile")

class TemplateFile(object):
	def __init__(self, input, output, **kwargs):
		self.input = input
		self.output = output
		self.kwargs = kwargs
		for name, value in kwargs.items():
			setattr(self, name, value)
		
	def run(self, args, opts, scope):
		self.save()
		
	def save(self):
		template_data = file(self.input).read()
		data = template_data.format(**self.kwargs)
		logger.info("writing to file: %s" % self.output)
		f = open(self.output, "w")
		f.write(data)
		f.close() 
		