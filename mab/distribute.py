# -*- coding: utf-8 -*-
from parallelize import *
import time
computerscript = os.path.expanduser("~/.computers.py")

if os.path.exists(computerscript):
	execfile(computerscript)
else:
	print "warning:", computerscript, "not found"
	groups = {}
	
def addgroup(name, filename, extrahost=None):
	lines = file(filename).readlines()
	lines = [line.strip() for line in lines if line.strip()]
	lines = [line for line in lines if line[0] != "#"]
	groups[name] = []
	for line in lines:
		cpucount = 1
		if ":" in line:
			hostname, cpucount = line.split(":")
			hostname = hostname.strip()
			cpucount = int(cpucount)
		else:
			hostname = line.strip()
		groups[name].append((hostname, cpucount, extrahost))
			
	
#print groups
#sys.exit(0)
#hostlist = [("virgo04", 2), ("virgo05", 2)]
#hostlist = [("virgo04", 1)]
hostlist = [] 
for nr in range(2,20):
	if nr not in [10,16,18]:
		hostname = "virgo%02d" % nr
		hostlist.append((hostname, 1))
		hostlist.append((hostname, 1))
		#os.system("ssh -q %s uname -a" % hostname)
	
#sys.exit()

def log(*args):
	#print >>sys.stderr, "LOG: " + " ".join([str(k) for k in args])
	pass

#class RemoteExecution(threading.Thread):
username = os.environ["USER"]
magicstartline = "start the distributed system now!"

class SshExecutor(IOExecutor):
	def __init__(self, thread, function, hostname, extrahost):
		self.hostname = hostname
		self.extrahost = extrahost
		IOExecutor.__init__(self, thread, function)
		
	def init(self):
		started = False
		blanks = 0
		while not started:
			command = self.input.readline().strip()
			#print "Eating up bogus lines", `command`
			if command == magicstartline:
				#print "great.. started the real thing"
				started = True
			elif len(command) == 0:
				blanks += 1
			if blanks == 4:
				print "failure?"
				return False
		return True
		#output.write("\n")
		#output.write(magicstartline)
		#output.write("\n")
		#IOExecutor.run(self)
	
	def createIO(self):
		command = " ".join(sys.argv)
		dir = os.getcwd()
		homedir = os.path.expanduser("~")
		if dir.startswith(homedir):
			dir = dir[len(homedir)+1:]
		#cmd = "ssh -q %s \"cd %s; /bin/nice -n 10 python %s --slave | tee /tmp/mabcluster_stdout_%s\"" % (self.hostname, dir, command, username)
		if self.extrahost:
			extra = "ssh " + self.extrahost + " "
		else:
		  	extra = " "
		nicecmd = "%s ssh -x -q %s /bin/nice -n 10 sleep 0" % (extra, self.hostname)
		niceok = True
		if os.system(nicecmd) != 0:
			print >>sys.stderr, "nice not supported: %s" % self.hostname
			niceok = False
		#else:
		
		if niceok:
			cmd = "ssh -x -q %s \"cd %s; /bin/nice -n 10 python %s --slaveside | tee /tmp/mabcluster_stdout_%s\"" % (self.hostname, dir, command, username)
		else:
			cmd = "ssh -x -q %s \"cd %s; python %s --slaveside | tee /tmp/mabcluster_stdout_%s\"" % (self.hostname, dir, command, username)
		if self.extrahost:
			cmd = extra + "\"" + cmd.replace("\"", "\\\"") + "\""
		print >>sys.stderr, cmd
		
		#os.system(cmd)
		self.output, self.input = os.popen2(cmd)

	def join(self):
		pass
		
class SshExecution(Execution):
	def __init__(self, hostname, taskQueue, extrahost, args, kwargs):
	#def __init__(self, hostname, cores, args, kwargs):
		Execution.__init__(self, taskQueue, False, None, args, kwargs)
		self.hostname = hostname
		self.extrahost = extrahost
		self.executor = SshExecutor(self, None, hostname, extrahost)

def distribute(groupname, cores=2, fork=True, flatten=False, info=False, *args, **kwargs):
	def server_wrapper(f):
		def execute(*multiargs):
			results = []
			len(zip(*multiargs))
			noArgs = len(zip(*multiargs))
			taskQueue = Queue.Queue(noArgs)
			#for timenr in range(times):
			#	taskQueue.put(timenr)
			for tasknr, _args in enumerate(zip(*multiargs)):
				taskQueue.put((tasknr, list(_args)))
			#for timenr in range(times):
			#	result = f(*args, **kwargs)
			#	results.append(result)
			hostlist = groups[groupname]
			executions = []
			if noArgs > 0:
				for hostname, _cores, extrahost in hostlist:
				  	if extrahost:
						extra = "ssh " + extrahost + " "
					else:
						extra = ""
						
					cmd = "%sssh -x -q %s sleep 0" % (extra, hostname)
					#print cmd, cores
					if os.system(cmd) != 0:
						print "Error connecting to hostname: %s" % hostname
					else:
						for _ in range(_cores):
							executions.append(SshExecution(hostname, taskQueue, extrahost, args, kwargs))
				executions.extend([Execution(taskQueue, fork, f, args, kwargs) for corenr in range(cores)])
			#executions.extend(remote_executions)
				#lala
				if info:
					InfoThread(len(multiargs[0]), taskQueue).start()
			for i, execution in enumerate(executions):
				#execution.setName("T-%d" % i)
				execution.start()
			for execution in executions:
				log("joining:",execution.getName())
				execution.join()
				results.extend(execution.results)
			log("sorting..")
			results.sort(cmp=lambda a, b: cmp(a[0], b[0]))
			results = [k[1] for k in results]
			if flatten:
				log("flattening..")
				flatresults = []
				for result in results:
					flatresults.extend(result)
				results = flatresults
			return results
		return execute
	def client_wrapper(f):
		def execute(*multiargs):
			done = False
			input = sys.stdin
			output = sys.stdout
			sys.stdout = sys.stderr
			started = False
			output.write("\n")
			output.write(magicstartline)
			output.write("\n")
			output.flush()
			#while not started:
			#	command = input.readline().strip()
			#	if command == magicstartline:
			#		started = True
			while not done:
				log("waiting for command")
				command = input.readline().strip()
				if command == "stop":
					log("stopping")
					done = True
					sys.exit(0)
				elif command == "execute":
					response = input.readline().strip()
					args, kwargs = deserialize(response)
					
					info = "ok"
					exc_info = None
					result = None
					try:
						result = f(*args, **kwargs)
					except Exception, e:
						info = "exception"
						done = True
						exc_info = traceback.format_exc()
					# encode result
					output.write(serialize((info, exc_info, result)))
					output.write("\n")
					output.flush()
				else:
					log("unkown command:", command)
					done = True
			sys.exit(0)
			#results.sort(cmp=lambda a, b: cmp(a[0], b[0]))
			#results = [k[1] for k in results]
			#if flatten:
			#	flatresults = []
			#	for result in results:
			#		flatresults.extend(result)
			#	results = flatresults
			#return results
		return execute
	slaveside = "--slaveside" in sys.argv
	if slaveside:
		return client_wrapper
	else:
		return server_wrapper

if __name__ == "__main__":
	#print "lalalalalaaa"
	@timed
	@distribute("go", cores=0, flatten=True)
 	#@parallelize(cores=8, fork=True, flatten=True, text="common argument")
	def testprime(from_nr, to_nr, text=""):
		#print text, from_nr, to_nr
		primes = []
		from_nr = max(from_nr, 2)
		for p in range(from_nr, to_nr):
			isprime = True
			for i in range(2, p):
				if p % i == 0:
					isprime = False
					break
			if isprime:
				primes.append(p)
		return primes
			
	
	testnumbers = range(0, 50001, 1000)
	testnumbers += range(5000, 70001, 2)
	from_nrs = testnumbers[:-1]
	to_nrs = testnumbers[1:]
	results = testprime(from_nrs, to_nrs)
		
	print len(results)
		