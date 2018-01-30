import problem_reader as pr
import random
import math

class GA:
    def __init__(self, fileName, population_size):
        self.problem_spec = pr.ProblemSpec()
        self.problem_spec.readProblemFile(fileName)
        self.population_size = population_size

        self.population = []


    def initializePopulation(self):
        for _ in range(0, self.population_size):
            genotype = Genotype(self.problem_spec)
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
    def __init__(self, problem_spec):
        max_vehicles = problem_spec.max_vehicles_per_depot * problem_spec.nr_depots
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

    def vehicleOverloaded(self, vehicle_route, vehicle_nr, customer_demand, problem_spec):
        demand_sum = customer_demand
        for customer in vehicle_route:
            demand_sum += customer.q
        depot_number = math.floor(vehicle_nr / problem_spec.nr_depots)
        if demand_sum > problem_spec.depots[depot_number - 1].Q:
            return True
        else:
            return False






ga = GA('p01', 10)
ga.initializePopulation()

yo = 5