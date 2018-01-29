import problem_reader as pr
import random
import math

class GA:
    def __init__(self, fileName, population_size):
        self.problem_object = pr.ProblemSpec()
        self.problem_object.readProblemFile(fileName)
        self.population_size = population_size

        self.population = []


    def initializePopulation(self):
        for _ in range(0, self.population_size):
            genotype = Genotype(self.problem_object)
            self.population.append(genotype)
        pass

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
    def __init__(self, problem_object):
        max_vehicles = problem_object.max_vehicles_per_depot * problem_object.nr_depots
        self.vehicle_routes = [ [] for i in range(0, max_vehicles)]
        random.seed()
        for customer in problem_object.customers:
            inserted = False
            while(not inserted):
                vehicle = random.randrange(max_vehicles)  # Choose random vehicle
                if not self.vehicleOverloaded(self.vehicle_routes[vehicle], vehicle, customer.q, problem_object): # Only append if it doesnt cause vehicle overload
                    self.vehicle_routes[vehicle].append(customer)
                    inserted = True

    def vehicleOverloaded(self, vehicle_route, vehicle, customer_demand, problem_object):
        demand_sum = customer_demand
        for customer in vehicle_route:
            demand_sum += customer.q
        depot_number = math.floor(vehicle / problem_object.nr_depots)
        if demand_sum > problem_object.depots[depot_number - 1].Q:
            return True
        else:
            return False






ga = GA('p01', 10)
ga.initializePopulation()

yo = 5