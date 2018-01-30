import numpy as np
import math
from scipy.spatial import distance

class ProblemSpec:
    def __init__(self, fileName):
        self.max_vehicles_per_depot = 0
        self.num_depots = 0
        self.num_customers = 0
        self.depots = []
        self.customers = []
        self.solution_cost = 0
        self.readProblemFile(fileName=fileName)

        self.num_stops = self.num_customers + self.num_depots
        self.cost_matrix = np.zeros((self.num_stops, self.num_stops))
        self.constructCostMatrix()

    def readProblemFile(self, fileName):
        with open('../data/Data Files/' + fileName, 'r') as f:
            main_data = list(map(int,(f.readline().strip('\n').split())))
            self.max_vehicles_per_depot, self.num_customers, self.num_depots = main_data
            for _ in range(0, self.num_depots):
                depot_data = list(map(int,(f.readline().strip('\n').split())))
                D, Q = depot_data
                self.depots.append(Depot(D, Q))
            for _ in range(0, self.num_customers):
                customer_data = list(map(int, (f.readline().strip('\n').split())))
                i, x, y, d, q = customer_data[:5]
                self.customers.append(Customer(i, x, y, d, q))
            for i in range(0, self.num_depots):
                depot_coords = list(map(int, (f.readline().strip('\n').split())))
                x, y = depot_coords[1:3]
                self.depots[i].x = x
                self.depots[i].y = y
        with open('../data/Solution Files/' + fileName + '.res', 'r') as f:
            solution_cost = float(f.readline().strip('\n'))
            self.solution_cost = solution_cost


    def constructCostMatrix(self):
        stops = self.customers + self.depots
        for j in range(0, len(stops)):
            for i in range(0, len(stops)):
                x_dist = abs(stops[i].x - stops[j].x)
                y_dist = abs(stops[i].y - stops[j].y)
                self.cost_matrix[i, j] = math.sqrt(math.pow(y_dist, 2) + math.pow(x_dist, 2))







class Depot:
    def __init__(self, D, Q):
        self.x = 0         # x coordinate
        self.y = 0         # y coordinate
        self.D = D         # max route duration
        self.Q = Q         # max vehicle load


class Customer:
    def __init__(self, i, x, y, d, q):
        self.i = i          # customer number
        self.x = x          # x coordinate
        self.y = y          # y coordinate
        self.d = d          # necessary service duration
        self.q = q          # demand for customer

problem_spec = ProblemSpec('p01')
yo = 5