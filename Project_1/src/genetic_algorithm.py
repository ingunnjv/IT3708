import problem_reader as pr
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import time


colors = ['crimson', 'green', 'blue',
        'gold', 'deeppink', 'aquamarine',  'blueviolet', 'brown',
        'chartreuse', 'coral',  'darkblue', 'darkcyan',
        'greenyellow', 'grey', 'aqua', 'olive',  'teal']

class GA:
    def __init__(self, fileName, population_size):
        self.problem_spec = pr.ProblemSpec(fileName)
        self.population_size = population_size

        self.population = []


    def initializePopulation(self):
        for _ in range(0, self.population_size):
            genotype = Genotype(self.problem_spec)
            self.population.append(genotype)

    def fitness(self):
        pass

    def parentSelection(self):
        pass

    def mutation(self, parent):
        pass

    def crossover(self, parent1, parent2):
        pass

    def survivorSelection(self):
        pass



class Genotype:
    def __init__(self, problem_spec):
        # Visualization data
        self.ax = None
        self.background = None
        self.fig = None

        # Genes
        max_vehicles = problem_spec.max_vehicles_per_depot * problem_spec.num_depots
        self.vehicle_routes = [ [] for i in range(0, max_vehicles)]
        random.seed()
        for customer in problem_spec.customers:
            inserted = False
            while(not inserted):
                vehicle_nr = random.randrange(max_vehicles)  # Choose random vehicle
                # Only append if it doesnt cause vehicle overload :
                if not self.vehicleOverloaded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer.q, problem_spec):
                    self.vehicle_routes[vehicle_nr].append(customer)
                    inserted = True
        self.problem_spec = problem_spec

    def vehicleOverloaded(self, vehicle_route, vehicle_nr, customer_demand, problem_spec):
        demand_sum = customer_demand
        for customer in vehicle_route:
            demand_sum += customer.q
        depot_number = self.getDepotNumber(vehicle_nr)
        if demand_sum > problem_spec.depots[depot_number - 1].Q:
            return True
        else:
            return False

    def vizualizeGenes(self):
        fig, ax = plt.subplots(1, 1)
        ax.set_aspect('equal')

        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        background = fig.canvas.copy_from_bbox(ax.bbox)

        for vehicle_nr, route in enumerate(self.vehicle_routes):
            route_x_coords = np.zeros(len(route) + 2)
            route_y_coords = np.zeros(len(route) + 2)
            depot_nr = self.getDepotNumber(vehicle_nr)
            route_x_coords[0] = self.problem_spec.depots[depot_nr].x
            route_y_coords[0] = self.problem_spec.depots[depot_nr].y
            for customer_num in range(0, len(route)):
                route_x_coords[customer_num + 1] = route[customer_num].x
                route_y_coords[customer_num + 1] = route[customer_num].y
            route_x_coords[-1] = route_x_coords[0]
            route_y_coords[-1] = route_y_coords[0]
            ax.plot(route_x_coords,route_y_coords, 'x-', color=colors[vehicle_nr % self.problem_spec.max_vehicles_per_depot], linewidth=0.8)

        plt.pause(1)
        self.background = background
        self.fig = fig
        self.ax = ax

    ###########################################
    # Get depot number of the specified vehicle
    def getDepotNumber(self, vehicle_nr):
        return math.floor(vehicle_nr / self.problem_spec.num_depots)

    ######################################################
    # Finds the duration (total cost) of a single solution
    def duration(self):
        duration = 0
        # Go through all routes
        for vehicle_nr, route in enumerate(self.vehicle_routes):
            depot_nr = self.getDepotNumber(vehicle_nr)
            # The depot is considered the very first "customer"
            previous_customer = depot_nr + self.problem_spec.num_customers
            for customer in route:
                # Find distance between the current customer and the previous
                duration += self.problem_spec.cost_matrix[previous_customer, customer]
                previous_customer = customer
            # Remember to add the distance from last customer to depot
            duration += self.problem_spec.cost_matrix[previous_customer, depot_nr + self.problem_spec.num_customers]
        return duration





ga = GA('p01', 10)
ga.initializePopulation()
#ga.population[0].vizualizeGenes()
#ga.population[0].vizualizeGenes()
#ga.population[1].vizualizeGenes()
#plt.close("all")

yo = 5