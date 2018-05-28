# -*- coding: utf-8 -*-
from types import MethodType

class Once(object):
	def __init__(self, func):
		self.func = func
		self.done = False
		self.result = None
		
	def __call__(self):
		if not self.done:
			#print "calling function"
			self.result = self.func()
			self.done = True
		#else:
		#	print "returning cached"
		return self.result
		
	def __get__(self, instance, owner):
		# make sure unbound method gets bound
		#self.f = MethodType(self.f, instance, owner)
		print self.func, instance, owner
		self.func = self.func.__get__(instance, owner)
		return self
		
		
class property_once(property):
	def __init__(self, *args, **kwargs):
		print "ctor", args
		super(property_once, self).__init__(*args, **kwargs)

	def __get__(self, instance, owner):
		#import pdb;
		#pdb.set_trace()
		try:
			return self.value
		except AttributeError:
			self.value = super(property_once, self).__get__(instance, owner)
		return self.value
		
	def ______set__(self, instance, value):
		print "super setter"
		super(property_once, self).__set__(instance, value)
		self.value = value
	
	def _setter(self, f):
		print "super setter", f
	#	def wrap(value):
			
			
		#self.value = value
		
    	
def once(func):
	# make sure a method gets called only once
	return Once(func)
	
	
if __name__ == "__main__":
	@once
	def load():
		print "loading..."
		return [1]
		
	a = load()
	print "a =", a
	b = load()
	print "b =", a
	assert a is b
	
	
	class Resource(object):
		@once
		def load(self):
			print "loading...", self
			return [2]
		
	
	res = Resource()
	a = res.load()
	print "a =", a
	b = res.load()
	print "b =", a
	assert a is b
	
	print res.load
	print res.load
	print res.load.func()
	
	
	if 1:
		print "#" * 70
		class Image(object):
			
			@property_once
			def data(self):
				print "loading data"
				return [3]
				
			#@data.setter
			#def data(self, value):
			#	print "saving", value
			#data = property_once(get_data)
				
		im = Image()
		a = im.data
		b = im.data
		print a, b
		assert a is b
		#im.data = [4]
		#print im.data
		
	
	