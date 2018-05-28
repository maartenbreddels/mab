from kaplot import *

class QHullPoints(object):
	def __init__(self, filename):
		self.filename = filename
		
	def load(self):
		lines = file(self.filename).readlines()
		self.dimension = int(lines[0].split()[0])
		points = int(lines[1].split()[0])
		self.points = []
		for line in lines[2:]:
			p = [float(k) for k in line.split()]
			self.points.append(p)
		self.array = array(self.points)
			
"""2
15 10 1
-10.101 -10.101
-0.1294381801544404 -0.07247409101984714
0.08267689532419747 -0.2397644955865706
0.1295260566906465 1.716033573116837
0.1740355150742391 0.5317519038435655
0.1851415205797575 0.3882545794457364
-0.9065939866848107 -0.2962957610652135
-0.1954805620516266 -0.07111892482963184
-0.1407581310832468 0.7233857048236082
-0.1676297826663962 0.2080621273999375
0.05868313821742954 0.06632066014880154
0.08806341399736994 0.1054080604689985
0.4761588899009253 -0.03168366595227294
3.094213357897477 -0.064721945677682
0.5410515627308725 0.2115615434955919
5 7 1 2 0 6
4 8 3 0 6
3 14 5 11
4 9 4 3 8
6 14 5 4 3 0 13
4 9 7 6 8
5 14 11 10 12 13
7 11 5 4 9 7 1 10
4 13 0 2 12
4 12 2 1 10
"""
import subprocess
import StringIO

class Test(object):
	def __init__(self, voronoi, parameterset): 
		self.voronoi = voronoi 
		self.parameterset = parameterset
		
	def run(self, *args):
		points = [k.values_org for k in self.parameterset.parameter_values]
		self.voronoi.process(points)
	

	 
class QHullVoronoi(object):
	def __init__(self, filename):
		self.filename = filename
		
	def process(self, points):
		output = "%d rbox 10 d2\n" % (len(points[0]))
		output += "%d\n" % (len(points))
		for p in points:
			output += " ".join([str(k) for k in p])
			output += "\n"
		print output
		p = subprocess.Popen(["qvoronoi", "o" "FA"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		out, err = p.communicate(output)
		lines = out.split("\n")
		print lines
		self.process_lines(lines)
		
	def load(self):
		lines = file(self.filename).readlines()
		self.process_lines(lines)
		
	def process_lines(self, lines):
		self.dimension = int(lines[0].split()[0])
		vertices = int(lines[1].split()[0])
		regions = int(lines[1].split()[1])
		self.points = []
		for line in lines[2:2+vertices]:
			print "vertex", line
			p = [float(k) for k in line.split()]
			self.points.append(p)
		self.array = array(self.points)
		self.regions = []
		for line in lines[2+vertices:2+vertices+regions]:
			print "region", line
			indices = [int(k) for k in line.split()[1:]]
			self.regions.append(indices)
			#p = [float(k) for k in line.split()]
			#self.points.append(p)
		
class QHullVoronoiParameterset(QHullVoronoi):
	def __init__(self, parameterset):
		QHullVoronoi.__init__(self, None)
		self.parameterset = parameterset 
		
	def load(self):
		points = [k.values_org for k in self.parameterset.parameter_values]
		self.original_points = array(points)
		self.process(points)
	
class Plot(object):
	def __init__(self, points, voronoi):
		self.points = points
		self.voronoi = voronoi
			
	def run(self, *args):
		self.points.load()
		self.voronoi.load()
		#dsa
		box()
		#end = 3
		#scatter(self.points.array[:,0][:],self.points.array[:,1][:], symbolsize="20pt")
		scatter(self.voronoi.original_points[:,0][:],self.voronoi.original_points[:,1][:], symbolsize="20pt")
		#scatter(self.voronoi.array[1:,0],self.voronoi.array[1:,1], symbolsize="20pt", color="red", symbolName="circle")
		for region in self.voronoi.regions[:]:#[2:4]:
			xs = []
			ys = []
			print region
			for index in region[:]:
				if index != 0:
				#if 1:
					xs.append(self.voronoi.array[index,0])
					ys.append(self.voronoi.array[index,1])
						#for i1, i2 in zip(region[:-1], region[1:]):
			#	if (i1 != 0) and (i2 != 0):
			#		#line(self.voronoi.array[i1,0], self.voronoi.array[i1,1], self.voronoi.array[i2,0], self.voronoi.array[i2,1])
				#xs.append(self.voronoi.array[index,0])
				#ys.append(self.voronoi.array[index,1])
			print xs, ys
			polyline(xs, ys, close=True)
		grow(1.1)
		draw()  