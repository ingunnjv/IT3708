import numpy as np
import math
import itertools

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
        self.max_cost = 0
        self.max_load = 0
        self.constructCostMatrix()
        self.findMaxCostAndLoad()

        self.swappable_customers = []
        self.constructCandidateListsAndSwappableCustomersList()


    def readProblemFile(self, fileName):
        with open('../data/Data Files/' + fileName, 'r') as f:
            main_data = list(map(int,(f.readline().strip('\n').split())))
            self.max_vehicles_per_depot, self.num_customers, self.num_depots = main_data
            for _ in range(0, self.num_depots):
                depot_data = list(map(int,(f.readline().strip('\n').split())))
                D, Q = depot_data
                if D == 0: D = float("Inf")
                self.depots.append(Depot(D, Q))
            for _ in range(0, self.num_customers):
                customer_data = list(map(int, (f.readline().strip('\n').split())))
                i, x, y, d, q = customer_data[:5] # only the first 5 datas of ustomers are used in this project
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

    def findMaxCostAndLoad(self):
        self.max_cost = self.cost_matrix.max()
        self.max_load = max([c.q for c in self.customers])

    def constructCandidateListsAndSwappableCustomersList(self):
        for c in self.customers:
            depots = list(range(0, self.num_depots))
            distances_to_depots = []

            # Find distance to each depot
            for depot_num in depots:
                customer_to_depot_cost = self.cost_matrix[c.i - 1, depot_num + self.num_customers]
                distances_to_depots.append(customer_to_depot_cost)
            sorted_depots = [x for _, x in sorted(zip(distances_to_depots, depots))]
            sorted_distances = sorted(distances_to_depots)

            # Add closest depot(s) to candidate list
            for i, depot in enumerate(sorted_depots):
                if sorted_distances[i] / sorted_distances[0] <= 3:
                    c.candidate_list.append(depot)
                else:
                    break

            # Add to swappable customers if swappable
            if len(c.candidate_list) > 1:
                self.swappable_customers.append(c)



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
        self.candidate_list = []    # list containing the nearest depot and any other depot close enough

    def __eq__(self, other):
        return self.i == other.i

    def __ne__(self, other):
        return self.i != other.i



if __name__ == "__main__":
    pr = problemSpec('p23')