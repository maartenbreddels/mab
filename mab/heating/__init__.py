from kaplot import *
from numpy import *
import numpy
# i index of room
# j index of 'other room'

class Section(object):
	def __init__(self, heat_capacity, T0, measured=None):
		self.heat_capacity = heat_capacity
		self.T0 = T0
		self.number = None # this will be provided by the 'System'
	
	def heating(self, t):
		return 0.
		
class Outside(Section):
	def __init__(self, T0, measured=None):
		Section.__init__(self, 1e100, T0, measured=measured)
		
Room = Section

class Wall(object):
	def __init__(self, room1, room2, alpha):
		self.room1 = room1
		self.room2 = room2
		self.alpha = alpha
		


class System(object):
	def __init__(self, rooms, connections):
		self.connections = connections
		self.room_count = len(rooms)
		self.rooms = rooms
		self.alpha = zeros((self.room_count, self.room_count))
		for i, room in enumerate(rooms):
			room.number = i
			
		for wall in self.connections:
			assert wall.room1.number != None
			assert wall.room2.number != None
			self.alpha[wall.room1.number, wall.room2.number] += wall.alpha
			self.alpha[wall.room2.number, wall.room1.number] += wall.alpha
		self.alpha = numpy.matrix(self.alpha).T
		#3	dsa
	def dTdt(self, T, t):
		"""
		dT_i/dt = \alpha_{i,j} * (T_j - T_i) + \sum_k \gamma_{i,k} P_{i, k}
		"""
		deltaT = zeros((self.room_count, self.room_count))
		for i in range(self.room_count):
			for j in range(self.room_count):
				deltaT[i,j] = -(T[i] - T[j])
		deltaT = numpy.matrix(deltaT)
		if 1:
			print "T=", T
			print "deltaT=", deltaT
			print "alpha=", self.alpha
			print deltaT*self.alpha
			print "..."
		
		dTdt = zeros((self.room_count))
		for i in range(self.room_count):
			for j in range(self.room_count):
				dTdt[i] += self.alpha[i,j] * deltaT[i,j]
		#dTdt = sum(deltaT*self.alpha, axis=1) # + array([room.heating(tindex) for room in self.rooms])
		#dTdt = array(dTdt.reshape(-1))[0]
		#print "    ", dTdt
		return array([dTdt[k]/self.rooms[k].heat_capacity for k in range(len(dTdt))])
		
	
	def run(self, args, opts, scope):
		T0 = array([room.T0 for room in self.rooms])
		t0 = 0.
		dt = 0.02/5
		N = 100
		temperatures = zeros((self.room_count, N+1))
		temperatures[:,0] = T0
		Ti = T0
		for i in range(N):
			t = t0 + i * dt
			#print "T", Ti
			Ti = Ti + self.dTdt(Ti, t) * dt
			print Ti
			temperatures[:,i+1] = Ti * 1
			
		for i, color in zip(range(self.room_count), nicecolors):
			graph(temperatures[i], color=color)
			
		draw()
			
			