import threading
import Queue
import logging
import cPickle as pickle
import sys
from optparse import OptionParser
import os
#import subprocess
from cStringIO import StringIO
import time

realstdout = sys.stdout
readstdin = sys.stdin

pickle_protocol = 2
pickle_encoding = "hex"

logging.basicConfig()
#logging.basicConfig(level=logging.FATAL, stream=sys.stderr)
#logger = logging.getLogger("cluster")

#console = logging.StreamHandler()
##formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
##console.setFormatter(formatter)
#logging.getLogger('').addHandler(console)
#logging.getLogger('cluster.master').addHandler(console)
#logging.getLogger('').addHandler(console)

#logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)
#
level = logging.DEBUG
#level = logging.INFO
username = "breddels"

def getCpuCount(lines):
	cpucount = 0
	for line in lines:
		if line.startswith("processor"):
			cpucount += 1
	return cpucount

def getCpuInfo(lines):
	cpucount = 0
	models = []
	freqs = []
	for line in lines:
		if line.startswith("model name"):
			models.append(line[line.index(":")+2:].strip())
		if line.startswith("cpu MHz"):
			freqs.append(float(line[line.index(":")+2:].strip()))
	return models, freqs
	
def localMachineBusy():
	return False
	stdin, stdout = os.popen2("who -l -u")
	lines = stdout.readlines()
	loggedin, loggedinx = getUse(lines)
	if loggedinx:
		return True
	if os.path.exists("/tmp/nocluster"):
		return True
	return False
	
def getUse(lines):
	loggedin = False
	loggedinx = False
	for line in lines:
		line = line.strip()
		loggedin = True
		parts = line.split()
		#print `parts[1]`
		if parts[1][0] == ":":
			#print "!!!"
			loggedinx = True
	#print loggedin, loggedinx
	return loggedin, loggedinx

class Processor(threading.Thread):
	def __init__(self, getjob, callback, executeJob):
		super(Processor, self).__init__(name="processorThread")
		self.getjob = getjob
		self.callback = callback
		self.executeJob = executeJob
		self.logger = logging.getLogger("cluster.processor")
		self.logger.setLevel(level)
		self.finished = False
		self.error = False
		self.done = False

	def run(self):
		try:
			while not self.done:
				self.logger.debug("asking for job")
				job = self.getjob()
				self.logger.debug("got job")
				if job is None:
					self.done = True
					self.logger.debug("finished")
					self.finished = True
					self.callback(None)
					self.logger.debug("finished(ret)")
				else:
					job.result = self.executeJob(job)
					self.callback(job)
		except Exception, e:
			self.logger.error("error in processing job:", exc_info=True)
			self.finished = True
			self.error = True
			self.done = True
			self.callback(None)
			raise
		
class Job(object):
	def __init__(self, *args, **kwargs):
		self.args = args
		self.kwargs = kwargs
		for arg, value in kwargs.items():
			setattr(self, arg, value)
		self.id = hash(self)

	def __repr__(self):
		return "Job(id=%d, %r, %r)" % (self.id, self.args, self.kwargs)

class Remote(threading.Thread):
	def __init__(self, hostname, getjob, callback, executeJob, processors=1, forceUse=False):
		threading.Thread.__init__(self)
		self.hostname = hostname
		self.getjob = getjob
		self.callback = callback
		self.executeJob = executeJob
		self.processors = processors
		self.logger = logging.getLogger("cluster.%s" % self.hostname)
		self.logger.setLevel(level)
		self.finished = False
		self.stopped = 0
		self.error = False
		self.done = False
		self.forceUse = forceUse
		
	def run(self):
		try:
			self._run()
		except Exception, e:
			self.logger.error("error in remote:", exc_info=True)
			self.finished = True
			self.error = True
			self.done = True
			self.callback(None)
			raise
	
	def _run(self):
	
		self.jobMap = {}
		#self.stdin, self.stdout, stderr = os.popen3(self.getCommand())
		self.stdin, self.stdout = os.popen2(self.getCommand())
		#self.popen = subprocess.Popen(self.getCommand(),
		#	stdout=subprocess.PIPE, stdin=subprocess.PIPE)
		#self.stdout, self.stdin = self.popen.stdout, self.popen.stdin
		
		once = False
		#buffer = StringIO()
		while not self.done:
			self.logger.debug("waiting for command")
			line = self.stdout.readline()
			#buffer.write(line)
			#print "buffer: %r" % buffer.getvalue()
			self.logger.debug("command: %r", line)
			if line.strip() == "stopped":
				self.stopped += 1
				if self.stopped == self.processors:
					self.logger.debug("all processor finished")
					self.finished = True
					self.done = True
					self.callback(None)
					self.logger.debug("finished(ret)")
				else:
					self.logger.debug("one processor finished")
			elif line.strip() == "error":
				self.done = True
				self.finished = True
				self.error = True
				self.callback(None)
			elif line.strip() == "busy":
				job = self.getjob(peek=True)
				if job is None:
					self.stdin.write("stop\n")
					self.stdin.flush()
				else:
					self.stdin.write("ok\n")
					self.stdin.flush()
			elif line.strip() == "getjob":
				job = self.getjob()
				self.logger.debug("got asked for a job, sending job: %r" % job)
				if job is None:
					self.stdin.write("stop\n")
					self.stdin.flush()
				else:
					self.jobMap[job.id] = job
					data = pickle.dumps((job.args, job.kwargs), pickle_protocol).encode(pickle_encoding)
					self.stdin.write("job=%d length=%d\n" % (job.id, len(data)))
					self.stdin.write(data)
					self.stdin.flush()
					self.logger.debug("send: job=%d length=%d\n%s" % (job.id, len(data), data))
			elif line.startswith("jobresult"):
				parts = line.split()
				hash = int(parts[0].split("=")[1])
				length = int(parts[1].split("=")[1])
				rawdata = self.stdout.read(length)
				#buffer.write(rawdata)
				data = rawdata.decode(pickle_encoding)
				result = pickle.loads(data)
				#self.logger.debug("finished job %r, result=%r", job, result)
				self.logger.debug("finished job %r", job)
				job = self.jobMap[hash]
				del self.jobMap[hash]
				job.result = result
				self.callback(job)
			else:
				#if once:
				#	return
				#else:
				self.logger.error("unknown command: %r", line)
				return

	def getCommand(self):
		#return "c:\\python25\python e:\\studie\\lss\\twopoint\\calcprimes.py --slave --processors=%d" % (self.processors)
		#return "ssh %s cd studie/lss/twopoint/calcprimes.py --slave --processors=%d" % (self.hostname, self.processors)
		#return "ssh %s \"cd studie/lss/twopoint; /bin/nice -n 10 python twopointcluster2.py --slave --processors=%d\"" % (self.hostname, self.processors)
		command = " ".join(sys.argv)
		#dir = "studie/lss/twopoint"
		dir = os.getcwd()
		#"ko"
		return "ssh -q %s \"cd %s; /bin/nice -n 10 python %s --slave --processors=%d | tee /tmp/mabcluster_stdout_%s\"" % (self.hostname, dir, command, self.processors, username)
		#return "dir c:\\"
		
	def who(self):
		self.logger.info("host: %s", self.hostname)
		os.system("ssh -q %s who -l -u" % self.hostname)
	
	def getCpuCount(self):
		stdin, stdout = os.popen2("ssh -q %s cat /proc/cpuinfo" % self.hostname)
		lines = stdout.readlines()
		return getCpuCount(lines)
	
	def getCpuInfo(self):
		stdin, stdout = os.popen2("ssh -q %s cat /proc/cpuinfo" % self.hostname)
		lines = stdout.readlines()
		return getCpuInfo(lines)
		
	def getUse(self):
		#logger.info("host: %s", self.hostname)
		stdin, stdout = os.popen2("ssh -q %s who -l -u" % self.hostname)
		lines = stdout.readlines()
		return getUse(lines)
		
	def mayUse(self):
		return self.forceUse or self.getUse()[1]

class Worker(object):
	def createJob(self):
		pass
		
	def executeJob(self, job):
		pass
		
	def returnJob(self, job):
		pass
		
	def setUp(self):
		pass
		
	def finished(self):
		pass
		
class Machine(object):
	def __init__(self, worker=None, createJob=None, executeJob=None, returnJob=None, pre=None, post=None, processors=1, opts=None):
		if worker != None:
			self.createJob = worker.createJob
			self.executeJob = worker.executeJob
			self.returnJob = worker.returnJob
			self.pre = worker.setUp
			self.post = worker.finish
		else:
			self.createJob = createJob
			self.executeJob = executeJob
			self.returnJob = returnJob
			self.pre = pre
			self.post = post
		self.processors = processors
		self.jobQueue = Queue.Queue()
		#self.finishedJobQueue = Queue.Queue()
		self.commandQueue = Queue.Queue()
		self.logger = logging.getLogger("cluster.master")
		self.logger.setLevel(level)
		self.opts = opts
		self.remotes = []
		
	def addRemote(self, hostname, processors=1, forceUse=False):
		if self.opts and not self.opts.noremote:
			self.remotes.append(Remote(hostname, self.getJob, self.finishedJob, self.executeJob, processors=processors, forceUse=forceUse))

	def getJob(self, peek=False):
		# called from other threads to get a job
		self.commandQueue.put(self._addJob)
		job = self.jobQueue.get()
		#if peek and job is not None: # if only peeking, put existing job back in queue
		if peek: # if only peeking, put existing job back in queue
			self.jobQueue.put(job)
		return job
		
	def _addJob(self):
		# create a job and add it to the job queue
		self.logger.debug("added job to job queue")
		job = None
		try:
			job = self._createJob()
		except Exception, e:
			self.logger.error("error creating job:", exc_info=True)
			raise
		self.jobQueue.put(job)

	def _createJob(self):
		# can be overridden, directly creates a job
		return self.createJob()
		
		
	#def executeJob(self, job):
	#	# execute the job, and add it to the finished job queue
	#	result = self.execute(*job.args, **job.kwargs)
	#	job.result = result
	#	self.finishedJobQueue.put(job)
	#	def retJob(job=job):
	#		self._returnJob(job):
	#	self.commandQueue.put(retJob)
		
	def finishedJob(self, job):
		if job:
			def retJob(job=job):
				try:
					self.returnJob(job)
					del job
				except Exception, e:
					self.logger.error("error returning job:", exc_info=True)
			self.commandQueue.put(retJob)
		else:
			self.commandQueue.put(None) # indicate a processor/remote finished
		
	#def _returnJob(self, job):
	#	self.return
		
	def start(self):
		self.logger.debug("starting")
		self.doPre()
		if self.opts.test:
			job = self._createJob()
			while job != None:
				job.result = self.executeJob(job)
				self.returnJob(job)
				job = self._createJob()
		else:
			processors = [Processor(self.getJob, self.finishedJob, self.executeJob) for k in range(self.processors)]
			remotes = [k for k in self.remotes if not k.mayUse()]
			remotenames = [k.hostname for k in self.remotes]
			print >>sys.stderr,"remote machines: " +", ".join(remotenames)
			self.logger.debug("remote machines:" +", ".join(remotenames))
			for processor in processors:
				processor.start()
			self.logger.debug("processors created and started, waiting for commands")
			for remote in remotes:
				remote.start()
			self.logger.debug("remotes started, waiting for commands")
				
			done = False
			while not done:
				try:
					command = self.commandQueue.get()
					if command:
						command()
					else:
						errors = len([k for k in processors if k.error])
						errors += len([k for k in remotes if k.error])
						if errors > 0:
							self.onError(processors, remotes)
						unfinished = len([k for k in processors if not k.finished])
						unfinished += len([k for k in remotes if not k.finished])
						self.logger.debug("unfinished processors/remotes: %d", unfinished)
						if unfinished == 0:
							done = True
				except Exception, e:
					self.logger.error("error returning job:", exc_info=True)
					sys.exit(1)
	
			while not self.commandQueue.empty():
				command = self.commandQueue.get()
				if command:
					command()
		self.doPost()
		return True
			
	def doPost(self):
		if self.post:
			self.post()
			
	def doPre(self):
		if self.pre:
			self.pre()

class Slave(Machine):
	def __init__(self, *args, **kwargs):
		Machine.__init__(self, *args, **kwargs)
		self.logger = logging.getLogger("cluster.slave")
		
	def doPost(self):
		pass

#	def doPre(self):
#		pass
		
	def addRemote(self, hostame, processors=1, forceUse=False):
		self.logger.debug("slaves shouldn't have remote nodes")
		
	def finishedJob(self, job):
		if job:
			def retJob(job=job):
				# send job back
				s = `job.result`
				self.logger.debug("pickling: %s" % s)
				try:
					data = pickle.dumps(job.result, pickle_protocol).encode(pickle_encoding)
				except:
					s = `job.result`
					self.logger.error("error pickling: %s" % s, exc_info=True)
					sys.stdout.flush()
					sys.stderr.flush()
					realstdout.flush()
					realstderr.flush()
					raise
				realstdout.write("jobresult=%d length=%d\n" % (job.id, len(data)))
				realstdout.write(data)
				realstdout.flush()
			self.commandQueue.put(retJob)
		else:
			def retStopped(job=job):
				self.logger.debug("sending stop")
				realstdout.write("stopped\n")
				try:
					realstdout.flush()
				except:
					self.logger.error("flush failed", exc_info=True)
			self.commandQueue.put(retStopped)
			self.commandQueue.put(None) # indicate a processor/remote finished

	def _createJob(self):
		try:
			if 1:
				busy = localMachineBusy()
				while busy:
					self.logger.debug("busy, waiting")
					realstdout.write("busy\n")
					realstdout.flush()
					line = sys.stdin.readline().strip()
					if line == "ok":
						time.sleep(5)
					else: # error or stop
						return None
					busy = localMachineBusy()
					
			realstdout.write("getjob\n")
			realstdout.flush()
			line = sys.stdin.readline().strip()
			if line == "stop":
				return None
			if line == "error":
				return None
			else:
					parts = line.split()
					id = int(parts[0].split("=")[1])
					print >>sys.stderr, "job with id", id
					length = int(parts[1].split("=")[1])
					data = sys.stdin.read(length).decode(pickle_encoding)
					args, kwargs = pickle.loads(data)
					job = Job(*args, **kwargs)
					job.id = id
					return job
		except Exception, e:
			self.logger.error("error", exc_info=True)
			realstdout.write("error\n")
			realstdout.flush()

	def onError(self, processors, remotes):
		self.logger.error("error in processor")
		for k in processors + remotes:
			k.done = True
		self.logger.error("waiting...")
		for k in processors + remotes:
			k.join()
		self.logger.error("all processes stopped")
		realstdout.write("error\n")
	
					
class Master(Machine):
	def onError(self, processors, remotes):
		self.logger.error("error in processor or remote")
		for k in processors + remotes:
			k.done = True
		for remote in remotes:
			self.jobQueue.put(None) 
		self.logger.error("waiting...")
		for k in processors + remotes:
			k.join(3)
		#sys.exit(1)
	

def getParser():
	parser = OptionParser("")
	parser.add_option("--cmd", dest="cmd", type="string", default=None, help="execute 'cmd' on remote machines")
	parser.add_option("--test", dest="test", action="store_true", default=False, help="like --noremote, but without multiple threads")
	parser.add_option("--noremote", dest="noremote", action="store_true", default=False, help="doesn't use remote slaves")
	parser.add_option("--who", dest="who", action="store_true", default=False, help="shows who is logged in on the remote machine")
	parser.add_option("--listfree", dest="listfree", action="store_true", default=False, help="list the 'free' machines (no X logins)")
	parser.add_option("--cpucount", dest="cpucount", action="store_true", default=False, help="get cpu count for remote machines")
	parser.add_option("--cpuinfo", dest="cpuinfo", action="store_true", default=False, help="get cpu info for remote machines")
	parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="no debug output")
	#parser.add_option("--seed", dest="seed", type="int", help="seed for the machine")
	parser.add_option("--slave", dest="slave", action="store_true", default=False, help="puts machine in slave mode(used internally)")
	parser.add_option("--processors", dest="processors", type="int", help="# processors for machine(used internally)")
	return parser


class Who(Master):
	def start(self):
		for remote in self.remotes:
			remote.who()

class ListFree(Master):
	def start(self):
		for remote in self.remotes:
			loggedin, loggedinx = remote.getUse()
			if not loggedinx:
				self.logger.info("host: %s is free", remote.hostname)
			else:
				self.logger.info("host: %s NOT is free", remote.hostname)

class CpuCount(Master):
	def start(self):
		stdin, stdout = os.popen2("cat /proc/cpuinfo")
		lines = stdout.readlines()
		cpucountLocal = getCpuCount(lines)
		cpucountTotal = cpucountLocal
		for remote in self.remotes:
			cpucount = remote.getCpuCount()
			self.logger.info("%s: %s cpus", remote.hostname, cpucount)
			cpucountTotal += cpucount
		self.logger.info("localhost: %s cpus", cpucountLocal)
		self.logger.info("total: %s cpus", cpucountTotal)

class CpuInfo(Master):
	def start(self):
		freqTotal = 0
		for remote in self.remotes:
			models, freqs = remote.getCpuInfo()
			for i, (model, freq) in enumerate(zip(models, freqs)):
				self.logger.info("%s:%i:\t%s\t%.2f Ghz", remote.hostname, i, model, freq)
				freqTotal += freq
		stdin, stdout = os.popen2("cat /proc/cpuinfo")
		lines = stdout.readlines()
		models, freqs = getCpuInfo(lines)
		for i, (model, freq) in enumerate(zip(models, freqs)):
			freqTotal += freq
			self.logger.info("localhost:%i:\t%s %s Ghz", i, model, freq)
		self.logger.info("total freq: %s Ghz", freqTotal)

class CmdExecute(Master):
	def start(self):
		freqTotal = 0
		for remote in self.remotes:
			#loggedin, loggedinx = remote.getUse()
			#if remote.mayUse():
				stdin, stdout = os.popen2("ssh -q %s %s" % (remote.hostname, self.opts.cmd))
				#print ">>> output:"
				remote.logger.info("output:")
				print stdout.read()
			

def createNode(argv=None, parser=None, *args, **kwargs):
	global level
	if argv is None:
		argv = sys.argv
	
	if parser is None:
		parser = getParser()
	opts, arguments = parser.parse_args(argv)
	#sys.stdout = StringIO()
	#sys.stdin = StringIO()
	if opts.quiet:
		level = logging.INFO
	if opts.slave:
		sys.stderr = open("/tmp/mabcluster_stderr_" +username, "w")
		if opts.processors:
			kwargs["processors"] = opts.processors
		return Slave(opts=opts, *args, **kwargs)
	if opts.cmd:
		return CmdExecute(opts=opts, *args, **kwargs)
	if opts.who:
		return Who(opts=opts, *args, **kwargs)
	elif opts.cpucount:
		return CpuCount(opts=opts, *args, **kwargs)
	elif opts.cpuinfo:
		return CpuInfo(opts=opts, *args, **kwargs)
	elif opts.listfree:
		return ListFree(opts=opts, *args, **kwargs)
	else:
		return Master(opts=opts, *args, **kwargs)
		
