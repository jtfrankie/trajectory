import pandas as pd
import numpy as np
from math import sin, cos, radians, degrees, tan, acos, atan, atan2, sqrt, fabs, pi
import xlrd
import xlwings as xw
import openpyxl

_tantra = None

def createTrajectory(strategy, md, inc, azi):
	if strategy == 'tangential':
		#if _tantra == None:
			#_tantra = TangentialTrajectory(md, inc, azi)
		#return _tantra
		return TangentialTrajectory(md, inc, azi)

	if strategy == 'balanced_tangential':
		return BalancedTangentialTrajectory(md, inc, azi)

	if strategy == 'average_angles':
		return AverageAnglesTrajectory(md, inc, azi)

	if strategy == 'vector_average':
		return VectorAverageTrajectory(md, inc, azi)

	if strategy == 'radii_of_curvature':
		return RadiiOfCurvatureTrajectory(md, inc, azi)

	if strategy == 'minimum_curvature':
		return MinimumCurvatureTrajectory(md, inc, azi)

def getAvailStrategies():
	return ['tangential', 'balanced_tangential', 'average_angles', 'vector_average', 'radii_of_curvature', 'minimum_curvature']

def getBestStrategy():
	return 'tangential'

class Trajectory():
	def __init__(self, md, inc, azi):
		self.md, self.inc, self.azi = md, inc, azi
		self.del_m, self.del_v, self.del_h, self.del_n, self.del_e = [], [], [], [], []
		self.vertical, self. northing, self.easting = [], [], []
		self.dls, self.cl_dep, self.cl_azi, self.cl_vs_dep = None, [], [], []
		self.clo_dist, self.clo_azi, self.vert_sec = [], [], []
		self.calculate_all()
		

	def store_data_to_sheet(self, sheet):
		sheet.range('A3:A167').options(transpose = True).value = self.md
		sheet.range('B3:B167').options(transpose = True).value = self.inc
		sheet.range('C3:C167').options(transpose = True).value = self.azi

		sheet.range('D3:D167').options(transpose = True).value = self.vertical
		sheet.range('E3:E167').options(transpose = True).value = self.northing
		sheet.range('F3:F167').options(transpose = True).value = self.easting

		if self.dls!= None:
			sheet.range('G3:G167').options(transpose = True).value = self.dls
			sheet.range('H3:H167').options(transpose = True).value = self.clo_dist
			sheet.range('I3:I167').options(transpose = True).value = self.clo_azi
			sheet.range('J3:J167').options(transpose = True).value = self.vert_sec

		else:
			sheet.range('G3:G167').options(transpose = True).value = self.clo_dist
			sheet.range('H3:H167').options(transpose = True).value = self.clo_azi
			sheet.range('I3:I167').options(transpose = True).value = self.vert_sec


	def calculate(self):
		pass

	def calculate_all(self):
		self.calculate()
		self.merge_NEV()
		self.closure_properties()


	def closure_properties(self):
		for i in range(len(self.northing)):
			target_azi = radians(-158.00)
			closure_distant = sqrt(self.northing[i] **2 + self.easting[i]**2)
			self.clo_dist.append(round(closure_distant, 2))

			if self.northing[i] == 0:
				closure_azi = 0
				self.clo_azi.append(round(closure_azi,2))
				vertical_sec = self.clo_dist[i] * cos(target_azi - self.clo_azi[i])
				
			else:
				closure_azi = atan2(self.northing[i], self.easting[i])
				
				if closure_azi <= radians(180) and closure_azi > radians(90):
					closure_azi = radians(450) - closure_azi
				elif closure_azi <= radians(90) and closure_azi >= radians(0):
					closure_azi = radians(90) - closure_azi
				else:
					closure_azi = radians(90) + fabs(closure_azi)

				self.clo_azi.append(closure_azi)
				vertical_sec = self.clo_dist[i] * cos(target_azi - self.clo_azi[i])

			
			self.clo_dist[i] = round(self.clo_dist[i], 2)
			self.clo_azi[i] = round(degrees(self.clo_azi[i]), 2)
			self.vert_sec.append(round(vertical_sec, 2))
		

	def calculate_delta_md(self):
		for i in range(len(self.md) - 1):
			delta_m = self.md[i + 1] - self.md[i]
			self.del_m.append(delta_m)

	def merge_NEV(self):
		vertical = 900
		northing = 0
		easting = 0

		for i in range(len(self.md)):
			if i == 0:
				self.vertical.append(vertical)
				self.northing.append(northing)
				self.easting.append(easting)

			else:
				vertical = vertical + self.del_v[i - 1]
				northing = northing + self.del_n[i - 1]
				easting = easting + self.del_e[i - 1]

				self.vertical.append(round(vertical, 2))
				self.northing.append(round(northing, 2))
				self.easting.append(round(easting, 2))

class TangentialTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()
		
		for i in range(len(self.md) - 1):
				delta_v = self.del_m[i] * cos(self.inc[i])
				delta_h = self.del_m[i] * sin(self.inc[i])

				self.del_v.append(float(delta_v))
				self.del_h.append(float(delta_h))

				delta_n = self.del_h[i] * cos(self.azi[i])
				delta_e = self.del_h[i] * sin(self.azi[i])
							
				self.del_n.append(float(delta_n))
				self.del_e.append(float(delta_e))


class BalancedTangentialTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()

		for i in range(len(self.md) - 1):
				delta_v = self.del_m[i]/2 * (cos(self.inc[i + 1]) + cos(self.inc[i]))
				self.del_v.append(float(delta_v))
				
				delta_n = self.del_m[i]/2 * sin(self.inc[i]) * cos(self.azi[i]) + self.del_m[i]/2 * sin(self.inc[i + 1]) * cos(self.azi[i + 1])
				delta_e = self.del_m[i]/2 * sin(self.inc[i]) * sin(self.azi[i]) + self.del_m[i]/2 * sin(self.inc[i + 1]) * sin(self.azi[i + 1])
				
				self.del_n.append(float(delta_n))
				self.del_e.append(float(delta_e))

		

class AverageAnglesTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()

		for i in range(len(self.md) - 1):
				delta_v = self.del_m[i] * cos((self.inc[i] + self.inc[i + 1])/2)
				delta_h = self.del_m[i] * sin((self.inc[i] + self.inc[i + 1])/2)
				
				self.del_v.append(float(delta_v))
				self.del_h.append(float(delta_h))

				average_azi = atan2(sin(self.azi[i]) + sin(self.azi[i + 1]), cos(self.azi[i]) + cos(self.azi[i + 1]))
				

				delta_n = self.del_h[i] * cos(average_azi)
				delta_e = self.del_h[i] * sin(average_azi)

				self.del_n.append(float(delta_n))
				self.del_e.append(float(delta_e))

		

class VectorAverageTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()

		first_survey_point, second_survey_point = np.array([0, 0, 0]), np.array([0, 0, 0])
		del_p = []
		for i in range(len(self.md) - 1):
			first_survey_point = np.array([cos(self.inc[i]), sin(self.inc[i]) * cos(self.azi[i]), sin(self.inc[i]) * sin(self.azi[i])])
			second_survey_point = np.array([cos(self.inc[i + 1]), sin(self.inc[i + 1]) * cos(self.azi[i + 1]), sin(self.inc[i + 1]) * sin(self.azi[i + 1])])

			delta_p = self.del_m[i] * np.divide(np.add(first_survey_point, second_survey_point), np.linalg.norm(np.add(first_survey_point, second_survey_point)))
			del_p.append(delta_p)

		for i in range(len(del_p)):
			self.del_v.append(del_p[i][0])
			self.del_n.append(del_p[i][1])
			self.del_e.append(del_p[i][2])

		


class RadiiOfCurvatureTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()

		radii_v, radii_h = [], []
		for i in range(len(self.md) - 1):
			if self.inc[i + 1] != self.inc[i]:
				rv = self.del_m[i]/(self.inc[i + 1] - self.inc[i])
				delta_v = rv * (sin(self.inc[i + 1]) - sin(self.inc[i])) 
				delta_h = rv * (cos(self.inc[i]) - cos(self.inc[i + 1])) 

			else:
				rv = 1
				delta_v = cos(self.inc[i + 1]) * self.del_m[i]
				delta_h = sin(self.inc[i + 1]) * self.del_m[i]

			radii_v.append(rv)

			self.del_v.append(delta_v)
			self.del_h.append(delta_h)

			if self.azi[i + 1] != self.azi[i]:
				rh = self.del_h[i]/(self.azi[i + 1] - self.azi[i])
				delta_n = rh * (sin(self.azi[i + 1]) - sin(self.azi[i]))
				delta_e = rh * (cos(self.azi[i]) - cos(self.azi[i + 1]))

			else:
				rh = 1
				delta_n = cos(self.azi[i + 1]) * self.del_h[i]
				delta_e = sin(self.azi[i + 1]) * self.del_h[i]

			radii_h.append(rh)

			self.del_n.append(delta_n)
			self.del_e.append(delta_e)

		
		self.dogleg(radii_v, radii_h)

	def dogleg(self, radii_v, radii_h):	
		self.dls = []
		self.dls.append(0)
		for i in range(len(self.md) - 1):
			dls = 5400/pi * sqrt((1/radii_v[i])**2 + (sin((self.inc[i] + self.inc[i +1])/2)**2/radii_h[i])**2)
			self.dls.append(round(dls, 2))

class MinimumCurvatureTrajectory(Trajectory):
	def __init__(self, md, inc, azi):
		Trajectory.__init__(self, md, inc, azi)

	def calculate(self):
		self.calculate_delta_md()
		del_p, radii_m = [], []

		for i in range(len(self.md) - 1):

			first_survey_point = np.array([cos(self.inc[i]), sin(self.inc[i]) * cos(self.azi[i]), sin(self.inc[i]) * sin(self.azi[i])])
			second_survey_point = np.array([cos(self.inc[i + 1]), sin(self.inc[i + 1]) * cos(self.azi[i + 1]), sin(self.inc[i + 1]) * sin(self.azi[i + 1])])

			alfa = 0.5 * np.arccos(np.dot(first_survey_point, second_survey_point))
			rm = self.del_m[i]/(2*alfa)
			radii_m.append(rm)

			delta_p = self.del_m[i]/2 * np.tan(alfa)/alfa * np.add(first_survey_point, second_survey_point)
			del_p.append(delta_p)

		for i in range(len(del_p)):
			self.del_v.append(del_p[i][0])
			self.del_n.append(del_p[i][1])
			self.del_e.append(del_p[i][2])
			

		
		self.dogleg(radii_m)
	
	def dogleg(self, radii_m):
		self.dls = []
		self.dls.append(0)
		for i in range(len(self.md) - 1):
			dls = 5400/pi * 1/radii_m[i]
			self.dls.append(round(dls, 2))
	

class Trajectories_Repository():
	def __init__(self):
		self.trajectories = {}
		self.md, self.inc, self.azi = [], [], []

	def read_file(self, fname):
		wb = openpyxl.load_workbook(fname)
		sheet = wb.get_sheet_by_name('Sheet1')
		directions = []
		for rowOfCellObjects in sheet['B6':'B170']:
			for cellObj in rowOfCellObjects:
				self.md.append(float(cellObj.value))

		for rowOfCellObjects in sheet['C6':'C170']:
			for cellObj in rowOfCellObjects:
				self.inc.append(float(radians(cellObj.value)))
		
		for rowOfCellObjects in sheet['D6':'D170']:
			for cellObj in rowOfCellObjects:
				directions.append(cellObj.value)

		for i in range(len(directions)):
			
			if directions[i][0] == 'N' and directions[i][-1] == 'E':
				azi = float(directions[i][1:-1])
				self.azi.append(radians(azi))

			elif directions[i][0] == 'S' and directions[i][-1] == 'E':
				azi = 180 - float(directions[i][1:-1])
				self.azi.append(radians(azi))

			elif directions[i][0] == 'S' and directions[i][-1] == 'W':
				azi = 180 + float(directions[i][1:-1])
				self.azi.append(radians(azi))

			elif directions[i][0] == 'N' and directions[i][-1] == 'W':
				azi = 360 - float(directions[i][1:-1])
				self.azi.append(radians(azi))
	
	def create_spreadsheet(self, newfile):
		#self.wb = xw.Book('Drilling_Project_1.xlsx')
		print('newfile', newfile, 'string', 'C:\\Users\\jtfra\\Desktop\\Python\\Drilling\\' + newfile )
		self.output_file = xw.Book('C:\\Users\\jtfra\\Desktop\\Python\\Drilling\\' + newfile)

	
		
	def calculate_trajectories(self):

		for strategy in getAvailStrategies():
			print('strategy', strategy)
			self.trajectories[strategy] = createTrajectory(strategy, self.md, self.inc, self.azi)
			sheet = self.output_file.sheets[strategy]
			self.trajectories[strategy].store_data_to_sheet(sheet)
		'''
		self.trajectories['balanced_tangential'] = BalancedTangentialTrajectory(self.md, self.inc, self.azi)
		sheet = self.output_file.sheets['balanced_tangential']
		self.trajectories['balanced_tangential'].store_data_to_sheet(sheet)

		self.trajectories['average_angles'] = AverageAnglesTrajectory(self.md, self.inc, self.azi)
		sheet = self.output_file.sheets['average_angles']
		self.trajectories['average_angles'].store_data_to_sheet(sheet)

		self.trajectories['vector_average'] = VectorAverageTrajectory(self.md, self.inc, self.azi)
		sheet = self.output_file.sheets['vector_average']
		self.trajectories['vector_average'].store_data_to_sheet(sheet)

		self.trajectories['radii_of_curvature'] = RadiiOfCurvatureTrajectory(self.md, self.inc, self.azi)
		sheet = self.output_file.sheets['radii_of_curvature']
		self.trajectories['radii_of_curvature'].store_data_to_sheet(sheet)

		self.trajectories['minimum_curvature'] = MinimumCurvatureTrajectory(self.md, self.inc, self.azi)
		sheet = self.output_file.sheets['minimum_curvature']
		self.trajectories['minimum_curvature'].store_data_to_sheet(sheet)
		'''


tra = Trajectories_Repository()
tra.read_file('Project_1.xlsx')
tra.create_spreadsheet('Drilling_Project_1.xlsx')


tra.calculate_trajectories()
