import mab.cluster
import time
from numpy import *
import os

class MemWorker(mab.cluster.Worker):
	def setUp(self):
	
		self.count = 100000
		self.index = 2
		self.keep = []
		
	def returnJob(self, job):
		ar = job.result
		#self.keep.append(ar)
		
	def executeJob(self, job):
		return arange(0,1e6)
	
	def createJob(self):
		if self.index >= self.count:
			return None
		else:
			self.index += 1
			return mab.cluster.Job(nr=self.index-1)
			
	def finish(self):
		print "done"

if __name__ == "__main__":
	worker = MemWorker()
	node = mab.cluster.createNode(worker=worker, processors=1)
	if 1:
		for i in range(1,15):
			if i not in [16,15,11,5]:
				name = "virgo%02d" % i
				if name != os.environ["HOST"]:
					node.addRemote(name, processors=2)
	#node.addRemote("virgo03", processors=2)
	#node.addRemote("virgo04", processors=2)
	#node.addRemote("hercules03", processors=3)
	#node.addRemote("hypatia", processors=4)
	node.start()
