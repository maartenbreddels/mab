# -*- coding: utf-8 -*-
import os
from mab.cvsfile import Record, CsvObject

class Column(object):
	def __init__(self, start, end, name, format, zerovalue=9999.9):
		self.start = start
		self.end = end
		self.name = name
		self.zerovalue = zerovalue
		if format[0] == "A":
			self.convert = lambda x: x.strip()
		elif format[0] == "I":
			self.convert = int
		elif format[0] == "F":
			self.convert = float
		else:
			raise Exception, "format does not start wiht A,F or I"
		
	def read(self, line):
		s = line[self.start-1:self.end]
		if len(s.strip()) == 0:
			return self.zerovalue
		else:
			return self.convert(s)

def read(file, headerstart, headerend, datastart, zerovalue=9999.9):
	if os.path.exists(file):
		lines = open(file).readlines()
	else:
		lines = file.split("\n")
	#lines = [line for line in lines]
	headerlines = lines[headerstart-1:headerend] # -1 for 0 offset counting
	datalines = lines[datastart-1:]
	#print headerlines
	
	columns = []
	for line in headerlines:
		if line[0:4].strip():
			start = int(line[0:4])
			end = int(line[5:8])
		else:
			start = int(line[5:8])
			end = int(line[5:8])
		name = line[22:30].strip()
		format = line[9:15]
		print start, end, `line[0:4]`, `line[5:8]`, `name`, name.strip()
		name = name.strip()
		columns.append(Column(start, end, name, format, zerovalue))
	
	records = []
	for line in datalines[:]:
		record = Record()
		for column in columns:
			#print `column.read(line)`,
			setattr(record, column.name, column.read(line))
		records.append(record)
		#print record
		
	return CsvObject(records)
		

class Simple(object):
	def __init__(self, filename, skipstart=0, skipend=None, printlines=False):
		self.filename = filename
		self.skipstart = skipstart
		self.skipend = skipend
		self.printlines = printlines
		
	def load(self):
		self.data = readsimple(self.filename, self.skipstart, self.skipend, self.printlines)
		
def readsimple(file, skipstart=0, skipend=None, printlines=False):
	if os.path.exists(file):
		lines = open(file).readlines()
	else:
		lines = file.split("\n")
	#lines = [line for line in lines]
	lines = lines[skipstart:]
	if skipend:
		lines = lines[:-(skipend)]
	datalines = lines[1:]
	#print headerlines
	
	#columns = []
	header = lines[0]
	header = header.strip()
	if header[0] == "#":
		header = header[1:].strip()
	columnnames = header.split()
	
	records = []
	for line in datalines[:]:
		if printlines:
			print line
		record = Record()
		for value_raw, name in zip(line.split(), columnnames):
			try:
				setattr(record, name, float(value_raw))
			except:
				setattr(record, name, value_raw)
		records.append(record)
	return CsvObject(records)		