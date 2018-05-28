# -*- coding: utf-8 -*-

class Four(object):
	def __init__(self, l1, l2, l3, l4):
		self.l1 = l1
		self.l2 = l2
		self.l3 = l3
		self.l4 = l4
		
		
	def R(self, p):
		return self.l1 + (p ** self.l3 - (1-p)**self.l4)  / self.l2
		
	def __call__(self, p):
		#p = self.R(x)
		a = self.l3 * p ** (self.l3 - 1.)  + self.l4 * (1-p) ** (self.l4 - 1)
		return self.l2 / a
		