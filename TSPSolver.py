#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

# Written by Jesse total time - 1.5 hours
	# Time complexity: O(n^2) - nearest neighbor * num cities
	# Space complexity: O(n) - route of cities
	def greedy( self,time_allowance=60.0 ):
		# Get the related information out of the class (make conceptually easier to understand)
		cities: List[City] = self._scenario.getCities()
		bssf = None

		# Spec assumes path exists for discussion of complexity, so this loop is not included in time complexity analysis
		startTime = time.time()
		solutionFound = False
		while not solutionFound and startTime - time.time() < time_allowance:

			# Pick a random city - constant time
			randCityIndex = random.randrange(0, len(cities) - 1)
			
			# Initialize route and start city
			# Time complexity: O(1)
			# Space complexity: O(n) Eventually the route and set of cities visited will include all cities 2n
			startCity = cities[randCityIndex]
			cityVisitedSet = set()
			route: List[City] = []
			cityVisitedSet.add(startCity)
			route.append(startCity)

			# Get the first nearest neighbor
			# Time complexity: O(n) - must traverse all edges in matrix (exist or not)
			# Space complexity: O(1)
			nearestNeighbor = self.getGreedyNeighbor(startCity, cityVisitedSet)

			# While there is a neighboring city that hasn't been visited and all cities haven't been visited
			# Time complexity: O(n) - runs once for each city (except first city)
			while nearestNeighbor != None and len(cityVisitedSet) != len(cities):

				# Mark the current city as visited
				cityVisitedSet.add(nearestNeighbor)
				route.append(nearestNeighbor)

				# Go to that city
				nearestNeighbor = self.getGreedyNeighbor(nearestNeighbor, cityVisitedSet)
			
			# If all cities haven't been visited, restart (failed because greedy path didn't work)
			if len(cityVisitedSet) == len(cities) and route[len(route) - 1].costTo(startCity) != math.inf:

				# Get the solution in the correct return format
				# Time complexity: O(n) - traverses city list to calculate cost
				# Space complexity: O(n) - holds a list with all cities
				bssf = TSPSolution(route)
				solutionFound = True

		endTime = time.time()

		results = {}
		results['cost'] = bssf.cost if solutionFound else math.inf 	# Cost of best solution
		results['time'] = endTime - startTime 						# Time spent to find the best solution
		results['count'] = 1 if solutionFound else 0 				# Total number of solutions found
		results['soln'] = bssf 										# The best solution found
		results['max'] = None 										# Null
		results['total'] = None 									# Null
		results['pruned'] = None 									# Null
		return results

	# Return the city's neighbor that is the closest (and not visited)
	# Return None if all reachable cities have been visited
	# Note that edges are represented in a matrix, not a list
	# Time complexity: O(n) - loop through all of the possible edges of a city
	# Space complexity: O(1) - assuming that the cities are stored elsewhere, no space added just accessing initial city and edge lists
	def getGreedyNeighbor(self, city: City, cityVisitedSet: Set) -> City:
		cities: List[City] = self._scenario.getCities()
		edges: List[List[Boolean]] = self._scenario._edge_exists

		# Loop through all possible city connections
		lowestCost = math.inf
		lowestCity = None

		# Go through all city connections
		# Time complexity: O(n)
		for i in range(0, len(edges[city._index])):

			# If there exists an edge between two cities
			# Time complexity: O(1)
			if edges[city._index][i]:
				cityCost = city.costTo(cities[i])
				if cityCost < lowestCost and not cities[i] in cityVisitedSet :
					lowestCost = cityCost
					lowestCity = cities[i]

		return lowestCity



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		pass



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		pass
