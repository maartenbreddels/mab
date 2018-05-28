# -*- coding: utf-8 -*-
import os
import subprocess

import mab.gd.logging as logging

logger = logging.getLogger("gd.utils.run")


class Run(object):
	def __init__(self, command, working_directory=None):
		self.command = command
		self.working_directory = working_directory
		
	def run(self, args, opts, scope):
		oldcwd = os.getcwd()
		try:
			if self.working_directory:
				logger.info("changing working directory to: %s" % self.working_directory)
				os.chdir(self.working_directory)
			logger.info("running %s" % self.command)
			os.system(self.command)
		finally:
			os.chdir(oldcwd)
		
class Batch(object):
	def __init__(self, commands):
		self.commands = commands
		
	def run(self, args, opts, scope):
		returncode = 0
		for command in self.commands:
			returncode = command.run(args, opts, scope)
			if returncode not in [0, None]:
				return -1
		return returncode	
		
class Print(object):
	def run(self, args, opts, scope):
		for arg in args[1:]:
			print "%s=%r" % (arg, scope[arg])
		return len(args[1:])
class Print2(object):
	def __init__(self, names):
		self.names = names
	def run(self, args, opts, scope):
		for name in self.names:
			print "%s=%r" % (name, scope[name])
		return 0	
		
class CommandMakeDirs(object):
	def __init__(self, directories):
		self.directories = directories
		
	def run(self, args, opts, scope):
		for dirname in self.directories:
			if not os.path.exists(dirname):
				logger.info("creating directory %s" % dirname)
				os.makedirs(dirname)
			else:
				logger.info("directory %s already exists" % dirname)
		
		
		
class MpiCommand(object):
	def __init__(self, binary, working_directory, cores, arguments, stdout=None, stderr=None, name="mpirun", nice=None):
		self.name = name
		self.binary = binary
		self.working_directory = working_directory
		self.cores=cores
		self.arguments = arguments
		self.stdout = stdout
		self.stderr = stderr
		self.nice = nice
		
		
	def run(self, args, opts, scope):
		oldcwd = os.getcwd()
		try:
			logger.info("changing working directory to: %s" % self.working_directory)
			os.chdir(self.working_directory)
			
			filename_check = "finished.%s" % self.name
			if self.nice:
				cmd = "nice -n%d mpirun -np %d %s %s" % (self.nice, self.cores, self.binary, " ".join(self.arguments))
			else:
				cmd = "mpirun -np %d %s %s" % (self.cores, self.binary, " ".join(self.arguments))
			#if os.path.exists(filename_check):
			#	logger.info("skipping cmd: %s" % cmd)
			#	return
			#else:
			#	logger.info("running command: %s" % cmd)
			logger.info("running command: %s" % cmd)
			
			
			if self.stdout:
				stdout = open(self.stdout, "w")
				filename = os.path.join(os.getcwd(), self.stdout)
				logger.info("redirecting stdout to: %s" % filename)
			else:
				stdout = None
			if self.stderr:
				if self.stderr == self.stdout:
					stderr = stdout
					logger.info("redirecting stderr to stdout")
				else:
					stderr = open(self.stderr, "w")
					filename = os.path.join(os.getcwd(), self.stderr)
					logger.info("redirecting stderr to: %s" % filename)
			else:
				stderr = None
			popen = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
			returncode = popen.wait()
			if returncode == 0:
				os.system("touch %s" % filename_check)
			else:
				logger.error("errorcode is: %d" % returncode)
			return returncode
			#os.system(cmd)
		finally:
			os.chdir(oldcwd)