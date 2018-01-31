import problem_reader as pr
import random
import math
import matplotlib.pyplot as plt
import numpy as np
from timeit import default_timer as timer
from operator import itemgetter


colors = ['crimson', 'green', 'blue',
        'gold', 'deeppink', 'aquamarine',  'blueviolet', 'brown',
        'chartreuse', 'coral',  'darkblue', 'darkcyan',
        'greenyellow', 'grey', 'aqua', 'olive',  'teal']

class GA:
    def __init__(self, fileName, population_size, generations, elite_ratio, tournament_ratio, crossover_prob, intra_mutation_prob, inter_mutation_prob):
        self.problem_spec = pr.ProblemSpec(fileName)

        # Parameters
        self.population_size = population_size
        self.generations = generations
        self.elite_ratio = elite_ratio
        self.tournament_ratio = tournament_ratio
        self.crossover_prob = crossover_prob
        self.intra_mutation_prob = intra_mutation_prob
        self.inter_mutation_prob = inter_mutation_prob

        # Data storage
        self.population = []
        self.average_fitness_history = []
        self.best_fitness_history = []

    ######################################################
    # Create an initial population by creating Genotype()'s and storing them in self.population
    def initializePopulation(self):
        start = timer()
        for _ in range(0, self.population_size):
            genotype = Genotype(self.problem_spec)
            self.population.append(genotype)
        end = timer()
        print("initalizePopulation timer: " + str(end-start))


    def intraMutation(self, offspring):
        pass

    def interMutation(self, offspring):
        pass

    ######################################################
    # Perform Best Cost Route Crossover on two parents. Produces two offspring
    def crossover(self, parent1, parent2):
        start = timer()
        # Randomly select a route from each parent, r1 in p1, r2 in p2
        total_num_vehicles = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
        route1 = parent1.vehicle_routes[random.randint(0, total_num_vehicles)]
        route2 = parent2.vehicle_routes[random.randint(0, total_num_vehicles)]

        # Remove all customers in r1 from p2 (pop(index))
        for c1 in route1:
            for vehicle_nr, route in enumerate(parent2.vehicle_routes):
                parent2.vehicle_routes[vehicle_nr] = [c2 for c2 in route if c1 != c2]

        # Remove all customers in r2 from p1
        for c2 in route2:
            for vehicle_nr, route in enumerate(parent1.vehicle_routes):
                parent1.vehicle_routes[vehicle_nr] = [c1 for c1 in route if c1 != c2]

        # List of position in p2 (index), insertion cost and feasibility[0,1] of the insertion
        insertion_cost = []

        # For each customer in r1:
        for route1_customer in route1:

            # Find insertion cost of the customer at each location in p2.
            for vehicle_nr, r in enumerate(parent2.vehicle_routes):
                depot_index = parent2.getDepotNumber(vehicle_nr) + parent2.problem_spec.num_customers
                prev_customer_index = depot_index

                for c_pos in range(0, len(r) + 1):
                    if c_pos < len(r): next_customer_index = r[c_pos].i - 1
                    else: next_customer_index = depot_index

                    cost = parent2.problem_spec.cost_matrix[prev_customer_index, route1_customer.i - 1] \
                           + parent2.problem_spec.cost_matrix[route1_customer.i - 1, next_customer_index]
                    prev_customer_index = next_customer_index

                    # Feasible?
                    feasible = (not (parent2.vehicleOverloaded(r, vehicle_nr, route1_customer.q, parent2.problem_spec)
                                or parent2.routeDurationLimitExceeded(r, vehicle_nr, route1_customer, c_pos)))
                    # Save the position in p2
                    pos = (vehicle_nr, c_pos)
                    insertion_cost.append((pos, cost, feasible))

            # Store in ordered list (increasing cost)
            insertion_cost = sorted(insertion_cost, key=itemgetter(1))
        # Random number (0,1)
        # If k <= 0.8:
            # choose first feasible insertion location
            # no such exists:
                # new route with this as only customer
        # else:
            # choose first entry in list

        # Repeat for r2.
        end = timer()
        print("Time: ", 300*(end - start))

    def survivorSelection(self):
        pass
    ######################################################
    # Returns the parent with the highest fitness from the parents parameter
    def getBestParent(self, parents):
        fitnesses = []
        for parent in parents:
            fitnesses.append(parent.fitness)
        return parents[fitnesses.index(max(fitnesses))]


    ######################################################
    # Return the top self.elites percent of the poplation
    def getPopulationElite(self):
        population = list(self.population)
        fitnesses = []
        for individual in population:
            fitnesses.append(individual.fitness)
        sorted_population = [x for _, x in sorted(zip(fitnesses, population))] #sorts population in descending order based on fitnesses
        elites_num = int(math.ceil(self.population_size * self.elite_ratio))
        elites = sorted_population[0:elites_num]
        return elites

    ######################################################
    # Return two genotypes from the population based on a tournament selection
    def parentSelection(self):
        num_individuals = int(math.ceil(self.tournament_ratio * self.population_size))
        if num_individuals < 2:
            num_individuals = 2

        selected = random.sample(range(0, self.population_size), num_individuals)
        tournament_population = []
        for i in selected:
            tournament_population.append(self.population[i])
        fitnesses = []
        for i in tournament_population:
            fitnesses.append(i.fitness)
        sorted_tournament_population = [x for _, x in sorted(zip(fitnesses, tournament_population))]  # sorts population in descending order based on fitnesses

        return sorted_tournament_population[0:2]

    ######################################################
    # Evolves an initial population by means of mutation and recombinations over a given number of generations
    def evolutionCycle(self):
        random.seed(None)

        # Generate the initial population
        self.initializePopulation()

        # Start the evolution process
        for generation in range(0, self.generations):
            # Evaluate the fitness of all individuals in the population and compute the average
            population_fitness = []
            for individual in self.population:
                population_fitness.append(individual.fitness)

            # Store fitness data
            average_fitness = sum(population_fitness) / self.population_size
            best_fitness = min(population_fitness)
            self.average_fitness_history.append(average_fitness)
            self.best_fitness_history.append(best_fitness)

            # Print generation data to terminal
            print("Generation: %d" % generation)
            print("Best fitness score: %f\n" % best_fitness)

            # Create a new population and directly copy the elites of the previous generation to the new one
            new_population = []
            elites = self.getPopulationElite()
            new_population += elites

            # Populate the new population until it has 200 individuals by doing recombination, mutation and survivors
            while len(new_population) < (self.population_size):
                # Do tournament selection to choose parents
                parents = self.parentSelection()

                # Roll dice on what event shall happen (recombination or survivors)
                event = random.random()
                if event < self.crossover_prob:
                    # Do crossover and pass offspring to next generation
                    offspring = self.crossover(parents[0], parents[1])

                else:
                    offspring = self.getBestParent(parents)

                # Roll dice on what event shall happen (inter-mutation, intra-mutation or no mutation)
                event = random.random()
                if event < self.inter_mutation_prob:
                    offspring = self.interMutation(offspring)
                elif (self.inter_mutation_prob < event and event < self.inter_mutation_prob + self.intra_mutation_prob):
                    offspring = self.intraMutation(offspring)
                new_population += offspring














class Genotype:
    def __init__(self, problem_spec):
        # Visualization data
        self.ax = None
        self.background = None
        self.fig = None

        # Randomly initialize the chromsome
        max_vehicles = problem_spec.max_vehicles_per_depot * problem_spec.num_depots
        self.vehicle_routes = [ [] for i in range(0, max_vehicles)]
        self.problem_spec = problem_spec
        random.seed()
        for customer in problem_spec.customers:
            # Choose depot based on tournament selection
            random_depots = random.sample(range(0, problem_spec.num_depots), int(math.ceil(problem_spec.num_depots * 0.60)))
            closest_depot_distance = float('Inf')
            closest_depot = random_depots[0]
            for depot_num in random_depots:
                if problem_spec.cost_matrix[customer.i - 1, depot_num] < closest_depot_distance:
                    closest_depot_distance = problem_spec.cost_matrix[customer.i - 1, depot_num]
                    closest_depot = depot_num

            inserted = False
            # Attempt to put the customer in the next depot if the selected one fails
            for depot in range(closest_depot, closest_depot + problem_spec.num_depots):
                depot = depot % problem_spec.num_depots
                vehicle_start_index = depot * problem_spec.max_vehicles_per_depot
            while(not inserted):
                vehicle_nr = random.randrange(max_vehicles)  # Choose random vehicle
                # Only append if it doesnt cause vehicle overload :
                if ((not self.vehicleOverloaded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer.q, problem_spec))
                    and (not self.routeDurationLimitExceeded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer, len(self.vehicle_routes[vehicle_nr])))):
                    self.vehicle_routes[vehicle_nr].append(customer)
                    inserted = True

                # Proceed by putting the customer in the first available vehicle
                for vehicle_nr in range(vehicle_start_index, vehicle_start_index + problem_spec.max_vehicles_per_depot - 1):
                    # Only append if it doesnt cause vehicle overload :
                    if not self.vehicleOverloaded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer.q, problem_spec):
                        self.vehicle_routes[vehicle_nr].append(customer)
                        inserted = True
                        break
                if inserted:
                    break


        self.fitness = self.duration()

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __gt__(self, other):
        return self.fitness > other.fitness

    ######################################################
    # Determine if a vehicle becomes overloaded if its assigned an additional customer, returns True/False
    def vehicleOverloaded(self, vehicle_route, vehicle_nr, customer_demand, problem_spec):
        demand_sum = customer_demand
        for customer in vehicle_route:
            demand_sum += customer.q
        depot_number = self.getDepotNumber(vehicle_nr)
        if demand_sum > problem_spec.depots[depot_number].Q:
            return True
        else:
            return False

    ######################################################
    # Determine if a route becomes too long if an additional customer is assigned in a certain position
    def routeDurationLimitExceeded(self, route, vehicle_nr, new_customer, new_customer_position):
        new_route = list(route)
        new_route.insert(new_customer_position, new_customer)

        duration = self.routeDuration(new_route, vehicle_nr)
        max_route_duration = self.problem_spec.depots[self.getDepotNumber(vehicle_nr)].D
        if duration > max_route_duration and max_route_duration != 0:
            return True
        else:
            return False

    ######################################################
    # Plots all the vehicle routes of the chromosome in a single plot
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
            duration += self.routeDuration(route, vehicle_nr)
        return duration

    ######################################################
    # Finds the duration of a single route
    def routeDuration(self, route, vehicle_nr):
        duration = 0
        # The depot is considered the very first "customer"
        depot_nr = self.getDepotNumber(vehicle_nr)
        prev_customer_index = depot_nr + self.problem_spec.num_customers
        for customer in route:
            customer_index = customer.i - 1
            # Find distance between the current customer and the previous
            duration += self.problem_spec.cost_matrix[prev_customer_index, customer_index]
            prev_customer_index = customer_index
        # Remember to add the distance from last customer to depot
        duration += self.problem_spec.cost_matrix[prev_customer_index, depot_nr + self.problem_spec.num_customers]

        return duration



ga = GA(fileName = 'p01', population_size = 500, generations = 1000,
        elite_ratio = 0.02, tournament_ratio = 0.15,
        crossover_prob = 0.6, intra_mutation_prob = 0.2, inter_mutation_prob = 0.25)


ga.evolutionCycle()
#ga.initializePopulation()
#dur = ga.population[0].duration()

#ga.population[0].vizualizeGenes()
#ga.population[1].vizualizeGenes()
#plt.close("all")

yo = 5