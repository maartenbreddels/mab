# -*- coding: utf-8 -*-
from numpy import *

class Record(object):
	def __init__(self):
		self.attributes = []
		
	def __getitem__(self, name):
		return self.__dict__[name]
	
	def __setitem__(self, name, value):
		self.__dict__[name] = value
		
	def clone(self):
		copy = Record()
		for attr in dir(self):
			if attr[0] != "_" and attr != "clone":
				setattr(copy, attr, getattr(self, attr))
				#print attr, getattr(copy, attr)
		copy.attributes = list(copy.attributes)
		return copy
			
	def _print(self):
		for attr in dir(self):
			if attr[0] != "_" and attr != "clone":
				print attr, getattr(self, attr)
				#print attr, getattr(copy, attr)
		
class CsvObject(object):
	def __init__(self, records):
		self.records = records
		self.cache = {}
		
	def __len__(self):
		return len(self.records)
		
	def __getitem__(self, i):
		if isinstance(i, slice):
			indices = i.indices(len(self))
			return CsvObject([self.records[i] for i in range(*indices)])
		elif isinstance(i, list) or hasattr(i, "__getitem__"):
			if hasattr(i, "dtype"):
				if i.dtype == bool:
					#return CsvObject([self.records[i] for i, b in enumerate(i) if b])
					return CsvObject([record for b, record in zip(i, self.records) if b])
			#indices = i.indices(len(self))
			return CsvObject([self.records[k] for k in i])
		else:
			return self.records[i]
		
	def __getattr__(self, name):
		if "records" not in self.__dict__:
			if name not in self.__dict__:
				raise AttributeError("bla")
			else:
				return self.__dict__[name]
		else:
			if name in dir(self.__dict__["records"][0]):
				if name in self.cache:
					return self.cache[name]
				else:
					ar = array([getattr(record, name) for record in self.records])
					self.cache[name] = ar
					return ar
			else:
				return self.__dict__[name]
		
	def __setattr__(self, name, values):
		if name == "records":
			self.__dict__[name] = values
		elif name == "cache":
			self.__dict__[name] = values
		else:
			for record, value in zip(self.__dict__["records"], values):
				record.__dict__[name] = value
				#pass
		#print "hooray" * 30
		#print dsadsada
			#record[name] = value
		
		
	def filter(self, pred):
		return CsvObject([k for k in self if pred(k)])
		
	def sort(self, cmp=None, key=None, reverse=False):
		#return CsvObject([k for k in self if pred(k)])
		#print self.records
		self.records.sort(key=key, reverse=reverse)
		self.cache = {} # empty cache
		#print self.records
		
	def clone(self):
		return CsvObject([k.clone() for k in self])
		
	def extend(self, object):
		if isinstance(object, CsvObject):
			self.records.extend(object.records)
		else:
			self.records.extend(object)
		self.cache = {} # empty cache
		
	
def readcsv(filename, columnMapping={}, maxlines=None, noeval=[], wrapper=False):
	file = open(filename)
	lines = file.readlines()
	lines = [line.strip() for line in lines]
	headers = [k.strip() for k in lines[0].split(",")]
	length = len(lines[1].split(","))
	if len(headers) != length:
		print "WARNNG, # headers and # columns in first row don't match"
	length = min(length, len(headers))
	records = []
	if maxlines:
		lines = lines[:maxlines+1]
	for line in lines[1:]:
		values = [k.strip() for k in line.split(",")]
		record = Record()
		for i, attrName in enumerate(headers):
			value = values[i]
			if not attrName in noeval:
				try:
					value = float(value)
				except:
					pass
			setattr(record, columnMapping.get(attrName, attrName), value)
			record.attributes.append(columnMapping.get(attrName, attrName))
		records.append(record)
	if wrapper:
		return CsvObject(records)
	else:
		return CsvObject(records)
		#return records
			
			
def writecsv(file, records, attributes=None, attributeMap={}):
	autoclose = False
	if isinstance(file, basestring):
		file = open(file, "w")
		autoclose = True  
	lines = []
	if attributes is None and len(records) > 0:
		attributes = records[0].attributes
		#attributes = [k for k in dir(records[0]) if k[0] != "_" and k not in ["clone", "attributes"]]
		#print attributes
	headers = [attributeMap.get(k,k) for k in attributes]
	lines.append(",".join(headers))
	for record in records:
		line = []
		for attr in attributes:
			line.append(str(getattr(record, attr)))
		lines.append(",".join(line))
	file.write("\n".join(lines))
	file.write("\n")
	if autoclose:
		file.close()

def writecolumns(file, columns, names):
	if names is not None:
		assert len(columns) == len(names) 
	autoclose = False
	if isinstance(file, basestring):
		file = open(file, "w")
		autoclose = True
	N = min([len(c) for c in columns])
	if names is not None:
		file.write(",".join(names))
		file.write("\n")
	for i in range(N):
		file.write(",".join([str(c[i]) for c in columns]))
		file.write("\n")
	if autoclose:
		file.close()
