import numpy as np
from itertools import chain

class ProblemSpec:
    def __init__(self):
        self.max_vehicles_per_depot = 0
        self.nr_depots = 0
        self.nr_customers = 0
        self.depots = []
        self.customers = []
        self.solution_cost = 0

    def readProblemFile(self, fileName):
        with open('../data/Data Files/' + fileName, 'r') as f:
            main_data = list(map(int,(f.readline().strip('\n').split())))
            self.max_vehicles_per_depot, self.nr_customers, self.nr_depots = main_data
            for _ in range(0, self.nr_depots):
                depot_data = list(map(int,(f.readline().strip('\n').split())))
                D, Q = depot_data
                self.depots.append(Depot(D, Q))
            for _ in range(0, self.nr_customers):
                customer_data = list(map(int, (f.readline().strip('\n').split())))
                i, x, y, d, q = customer_data[:5]
                self.customers.append(Customer(i, x, y, d, q))
            for i in range(0, self.nr_depots):
                depot_coords = list(map(int, (f.readline().strip('\n').split())))
                x, y = depot_coords[1:3]
                self.depots[i].x = x
                self.depots[i].y = y
        with open('../data/Solution Files/' + fileName + '.res', 'r') as f:
            solution_cost = float(f.readline().strip('\n'))
            self.solution_cost = solution_cost



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

problem_spec = ProblemSpec()
problem_spec.readProblemFile('p01')