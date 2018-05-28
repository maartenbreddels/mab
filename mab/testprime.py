import mab.cluster
import time

class PrimeWorker(mab.cluster.Worker):
	def setUp(self):
	
		self.count = 100
		self.index = 2
		self.primes = []
		
	def returnJob(self, job):
		nr, prime = job.result
		if prime:
			print "%d is a prime" % nr
			self.primes.append(nr)
		
	def executeJob(self, job):
		#time.sleep(1)
		for i in range(2, job.nr-1):
			if job.nr % i == 0:
				return job.nr, False
		return job.nr, True
	
	def createJob(self):
		if self.index >= self.count:
			return None
		else:
			self.index += 1
			return mab.cluster.Job(nr=self.index-1)
			
	def finish(self):
		primes = ", ".join([str(k) for k in self.primes])
		print "all primes:", primes

if __name__ == "__main__":
	worker = PrimeWorker()
	node = mab.cluster.createNode(worker=worker, processors=1)
	for i in range(1,21):
		#if i not in [11]:
			node.addRemote("virgo%02d" % i, processors=2)
	#node.addRemote("virgo03", processors=2)
	#node.addRemote("virgo04", processors=2)
	#node.addRemote("hercules03", processors=1)
	node.addRemote("hercules03", processors=1)
	node.addRemote("herschel", processors=1)
	#node.addRemote("hypatia", processors=4)
	node.start()
