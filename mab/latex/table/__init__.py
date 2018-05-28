# -*- coding: utf-8 -*-

class Table2(object):
	def __init__(self, scopes):
		self.scopes = scopes
		
		
	def run(self, args, opts, scope):
		for scope in self.scopes:
			print scope.name