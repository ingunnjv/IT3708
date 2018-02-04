import problem_reader as pr
import random
import math
import matplotlib.pyplot as plt
import numpy as np
import copy
from timeit import default_timer as timer
from operator import itemgetter
import copy
import itertools

import _pickle as cPickle

# colors = ['crimson', 'green', 'blue',
#         'gold', 'deeppink', 'aquamarine',  'blueviolet', 'brown',
#         'chartreuse', 'coral',  'darkblue', 'darkcyan',
#         'greenyellow', 'grey', 'aqua', 'olive',  'teal']

depot1_colors = ['crimson', 'coral', 'red', 'tomato']
depot2_colors = ['grey', 'blue', 'darkblue', 'teal']
depot3_colors = ['green', 'darkgreen', 'seagreen', 'lime']
depot4_colors = ['chartreuse', 'gold', 'yellow', 'olive']

colors = []
colors.append(depot1_colors)
colors.append(depot2_colors)
colors.append(depot3_colors)
colors.append(depot4_colors)


class GA:
    def __init__(self, fileName, population_size, generations, elite_ratio, tournament_ratio, crossover_prob, intra_mutation_prob, inter_mutation_prob,
                 inter_mutation_attempt_rate):
        self.problem_spec = pr.ProblemSpec(fileName)

        # Parameters
        self.population_size = population_size
        self.generations = generations
        self.elite_ratio = elite_ratio
        self.tournament_ratio = tournament_ratio
        self.crossover_prob = crossover_prob
        self.intra_mutation_prob = intra_mutation_prob
        self.inter_mutation_prob = inter_mutation_prob
        self.inter_muation_attempt_rate = inter_mutation_attempt_rate

        # Data storage
        self.population = []
        self.average_fitness_history = []
        self.best_fitness_history = []

        # Workload timers
        self.crossover_time = 0
        self.intra_mutation_time = 0
        self.inter_mutation_time = 0
        self.update_fitness_time = 0
        self.initial_loop_time = 0
        self.deep_copy_time = 0

    ######################################################
    # Create an initial population by creating Genotype()'s and storing them in self.population
    def initializePopulation(self):
        start = timer()
        for _ in range(0, self.population_size):
            genotype = Genotype()
            genotype.initGenes(self.problem_spec)
            genotype.tooManyCustomers(self.problem_spec)
            self.population.append(genotype)
        end = timer()
        print("initalizePopulation timer: " + str(end-start))


    ######################################################
    # Applies a change to one random vehicle route (inversion)
    def intraMutation(self, offspring, type):
        start = timer()
        if type == "inversion":
            num_vehicles = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
            range_num_vehicles = list(range(0, num_vehicles))
            vehicle_route = None
            attempts = 1
            while range_num_vehicles:
                vehicle_nr = range_num_vehicles.pop(random.randint(0, num_vehicles - attempts))
                vehicle_route = offspring.vehicle_routes[vehicle_nr]
                if len(vehicle_route) >= 2:
                    vehicle_route_length = len(vehicle_route)
                    cutpoints = random.sample(range(0, vehicle_route_length), 2)
                    material = vehicle_route[min(cutpoints): max(cutpoints) + 1]
                    material.reverse()
                    vehicle_route[min(cutpoints): max(cutpoints) + 1] = material

                    new_route_duration = offspring.routeDuration(vehicle_route, vehicle_nr, self.problem_spec)
                    depot_nr = offspring.getDepotNumber(vehicle_nr, self.problem_spec)
                    if new_route_duration > self.problem_spec.depots[depot_nr].D:
                        material.reverse()
                        vehicle_route[min(cutpoints): max(cutpoints) + 1] = material
                        attempts += 1
                    else:
                        break
                else:
                    attempts += 1

        offspring.tooManyCustomers(self.problem_spec)
        end = timer()
        self.intra_mutation_time += (end - start)


    ######################################################
    # Applies a swapping of one customer from one vehicle route to another
    def interMutation(self, offspring):
        start = timer()
        num_vehicles = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
        random_vehicle_route = None
        # Acquire a non-empty vehicle route
        while not random_vehicle_route:
            random_vehicle_route_nr = random.randint(0, num_vehicles - 1)
            random_vehicle_route = offspring.vehicle_routes[random_vehicle_route_nr]
        random_customer_nr = random.randint(0, len(random_vehicle_route) - 1)
        customer = random_vehicle_route[random_customer_nr]

        range_num_vehicles = list(range(0, num_vehicles))
        attempts = 0
        # Attempt to insert the customer into a random vehicle route
        while range_num_vehicles:
            vehicle_nr = range_num_vehicles.pop(random.randint(0, num_vehicles - attempts - 1))

            if ((not offspring.vehicleOverloaded(offspring.vehicle_routes[vehicle_nr], vehicle_nr, customer.q, self.problem_spec))
                and (not offspring.routeDurationLimitExceeded(offspring.vehicle_routes[vehicle_nr], vehicle_nr, customer,
                                                         len(offspring.vehicle_routes[vehicle_nr]), self.problem_spec))):
                random_vehicle_route.pop(random_customer_nr)
                random_placement = random.randint(0, len(offspring.vehicle_routes[vehicle_nr]))
                offspring.vehicle_routes[vehicle_nr].insert(random_placement, customer)
                break
            attempts += 1
        offspring.tooManyCustomers(self.problem_spec)
        end = timer()
        self.inter_mutation_time += (end - start)


    ######################################################
    # Perform Best Cost Route Crossover on two parents. Produces two offspring
    def crossover(self, parent1, parent2):
        random.seed()
        start = timer()

        # Initialize offspring
        start = timer()
        offspring1 = copy.deepcopy(parent1)
        offspring2 = copy.deepcopy(parent2)
        end = timer()
        self.deep_copy_time += (end-start)

        # Randomly select a route from each parent, r1 in p1, r2 in p2
        total_num_vehicles = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
        route1 = offspring1.vehicle_routes[random.randint(0, total_num_vehicles - 1)]
        route2 = offspring2.vehicle_routes[random.randint(0, total_num_vehicles - 1)]

        ## What if a route is empty? No change will happen in the other parent... ##

        # Remove all customers in r1 from p2
        for c1 in route1:
            for vehicle_nr, route in enumerate(offspring2.vehicle_routes):
                offspring2.vehicle_routes[vehicle_nr] = [c2 for c2 in route if c1 != c2]

        # Remove all customers in r2 from p1
        for c2 in route2:
            for vehicle_nr, route in enumerate(offspring1.vehicle_routes):
                offspring1.vehicle_routes[vehicle_nr] = [c1 for c1 in route if c1 != c2]

        # Find location and insert customers from route 1 into offspring 2
        for route1_customer in route1:
            cost = self.insertionCost(route1_customer, offspring2)
            offspring2 = self.insertCustomerInRoute(route1_customer, offspring2, cost)

        # Find location and insert customers from route 2 into offspring 1
        for route2_customer in route2:
            cost = self.insertionCost(route2_customer, offspring1)
            offspring1 = self.insertCustomerInRoute(route2_customer,offspring1, cost)


        offspring1.tooManyCustomers(self.problem_spec)
        offspring2.tooManyCustomers(self.problem_spec)

        end = timer()
        self.crossover_time += (end-start)
       # print("Time: ",(end - start))

        return offspring1, offspring2

    ######################################################
    # Calculate the cost of adding a new customer to each feasible location in an individual's routes
    def insertionCost(self, new_customer, individual):
        # List of position (index) and insertion cost of route_customer into all locations in the individual's routes
        insertion_cost = []

        # Find insertion cost of the customer at each location in the individual.
        for vehicle_nr, r in enumerate(individual.vehicle_routes):
            depot_index = individual.getDepotNumber(vehicle_nr, self.problem_spec) + self.problem_spec.num_customers
            prev_customer_index = depot_index

            for c_i in range(0, len(r) + 1):
                if c_i < len(r): next_customer_index = r[c_i].i - 1
                else: next_customer_index = depot_index

                cost = self.problem_spec.cost_matrix[prev_customer_index, new_customer.i - 1] \
                       + self.problem_spec.cost_matrix[new_customer.i - 1, next_customer_index]
                prev_customer_index = next_customer_index

                # Feasible?
                feasible = (not (individual.vehicleOverloaded(r, vehicle_nr, new_customer.q, self.problem_spec)
                                 or individual.routeDurationLimitExceeded(r, vehicle_nr, new_customer, c_i, self.problem_spec)))
                # Save the position
                #if feasible:
                location = (vehicle_nr, c_i)
                insertion_cost.append((location, cost, feasible))

        # Store in ordered list (increasing cost)
        insertion_cost = sorted(insertion_cost, key=itemgetter(1))
        return insertion_cost


    ######################################################
    # Add a new customer to a location in another individual's routes
    def insertCustomerInRoute(self, new_customer, individual, insertion_cost):
        k = random.random()

        # Choose best insertion location
        #TODO: Handle when there is no feasible insertions
        if k <= 0.8:
            inserted = False
            for i in range(len(insertion_cost)):
                if insertion_cost[i][2]:
                    insertion_location = insertion_cost[i][0]
                    inserted = True
                    break
            if not inserted:
                insertion_location = insertion_cost[random.randint(0, len(insertion_cost) - 1)][0]
        # Choose random feasible entry in list
        else: insertion_location = insertion_cost[0][0]

        # Insert the new customer in the chosen route and position in the individual
        insertion_vehicle, insertion_customer_pos = insertion_location
        individual.vehicle_routes[insertion_vehicle].insert(insertion_customer_pos, new_customer)

        return individual



    def survivorSelection(self):
        pass

    ######################################################
    # Returns the parent with the highest fitness from the parents parameter
    def getBestParent(self, parents):
        fitnesses = []
        for parent in parents:
            fitnesses.append(parent.fitness)
        start = timer()
        best_parent = copy.deepcopy(parents[fitnesses.index(max(fitnesses))])
        end = timer()
        self.deep_copy_time += (end-start)
        return best_parent


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
        start = timer()
        elites = copy.deepcopy(elites)
        end = timer()
        self.deep_copy_time += (end - start)
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
        for generation in range(1, self.generations):
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
            print("Best fitness score: %f" % best_fitness)
            print("Average fitness score: %f" % average_fitness)

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
                    offsprings = self.crossover(parents[0], parents[1])

                else:
                    start = timer()
                    offsprings = copy.deepcopy(parents)
                    end = timer()
                    self.deep_copy_time += (end - start)

                for offspring in offsprings:
                    # Roll dice on what event shall happen (inter-mutation, intra-mutation or no mutation)
                    event = random.random()
                    if event < self.inter_mutation_prob and generation % self.inter_muation_attempt_rate == 0:
                        self.interMutation(offspring)
                    elif self.inter_mutation_prob < event and event < self.inter_mutation_prob + self.intra_mutation_prob:
                        self.intraMutation(offspring, 'inversion')

                    start = timer()
                    offspring.updateFitnessVariables(self.problem_spec)
                    offspring.updateFitness(self.problem_spec)
                    end = timer()
                    self.update_fitness_time += (end-start)
                    new_population += [offspring]

            # Stop population growing in case everything is not even numbers
            if len(new_population) > (self.population_size):
                new_population.pop()

            self.population = new_population



            print("Crossover work: %f [s]" % self.crossover_time)
            print("Intra-mutation work: %f [s]" % self.intra_mutation_time)
            print("Inter-mutation work: %f [s]" % self.inter_mutation_time)
            print("Update fitness work: %f [s]" % self.update_fitness_time)
            print("Deep copy work: %f [s]\n" % self.deep_copy_time)
            self.crossover_time = 0
            self.intra_mutation_time = 0
            self.inter_mutation_time = 0
            self.update_fitness_time = 0
            self.deep_copy_time = 0
        self.population[0].visualizeGenes(self.problem_spec)
        self.population[0].printSolutionData(self.problem_spec)

class Genotype:
    def __init__(self):
        # Visualization data
        self.ax = None
        self.background = None
        self.fig = None

        self.fitness = float("Inf")
        self.vehicle_routes = None

        self.demand_ol = 0           # demand overload
        self.duration_ol = 0   #duration_overload
        self.duration = 0

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __gt__(self, other):
        return self.fitness > other.fitness

    def __deepcopy__(self, memodict={}):
        return cPickle.loads(cPickle.dumps(self, -1))

    def initGenes(self,  problem_spec):
        # Initialization of placement of customers to vehicle routes
        max_vehicles = problem_spec.max_vehicles_per_depot * problem_spec.num_depots
        self.vehicle_routes = [[] for i in range(0, max_vehicles)]

        self.closestDepotInit(problem_spec)
        self.updateFitnessVariables(problem_spec)
        self.updateFitness(problem_spec)

    def randomInit(self):
        pass

    def closestDepotInit(self, problem_spec):
        random.seed()
        range_num_customers = list(range(0, len(problem_spec.customers)))
        customers_placed = 0
        while range_num_customers:
            customer_nr = range_num_customers.pop(random.randint(0, len(problem_spec.customers) - customers_placed - 1))
            customer = problem_spec.customers[customer_nr]
            # Choose depot based on tournament selection
            random_depots = random.sample(range(0, problem_spec.num_depots),
                                          int(math.ceil(problem_spec.num_depots * 0.60)))
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
                # Proceed by putting the customer in the first available vehicle
                for vehicle_nr in range(vehicle_start_index,
                                        vehicle_start_index + problem_spec.max_vehicles_per_depot):
                    # Only append if it doesnt cause vehicle overload :
                    if ((not self.vehicleOverloaded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer.q,
                                                    problem_spec))
                        and (not self.routeDurationLimitExceeded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer,
                                                                 len(self.vehicle_routes[vehicle_nr]), problem_spec))):
                        self.vehicle_routes[vehicle_nr].append(customer)
                        inserted = True
                        break
                if inserted:
                    break
            # If no feasible position found, insert in random position
            if not inserted:
                random_vehicle_number = random.randint(vehicle_start_index, vehicle_start_index + problem_spec.max_vehicles_per_depot - 1)
                self.vehicle_routes[random_vehicle_number].append(customer)

            customers_placed += 1




    ######################################################
    # Calculate the amount of overload of a solution.
    # - Duration longer than max allowed duration of a route
    # - Load above max allowed load in a route
    def updateFitnessVariables(self, problem_spec):
        self.duration_ol = 0
        self.demand_ol = 0
        self.duration = 0
        for vehicle_nr, route in enumerate(self.vehicle_routes):
            depot_num = self.getDepotNumber(vehicle_nr, problem_spec)
            route_duration = self.routeDuration(route, vehicle_nr, problem_spec)
            self.duration += route_duration

            duration_ol = route_duration - problem_spec.depots[depot_num].D
            demand_ol = self.routeDemand(route) - problem_spec.depots[depot_num].Q

            if duration_ol > 0: self.duration_ol += duration_ol
            if demand_ol > 0: self.demand_ol += demand_ol



    ######################################################
    # Determine if a vehicle becomes overloaded if its assigned an additional customer, returns True/False
    def vehicleOverloaded(self, vehicle_route, vehicle_nr, customer_demand, problem_spec):
        demand_sum = self.routeDemand(vehicle_route)
        demand_sum += customer_demand
        depot_number = self.getDepotNumber(vehicle_nr, problem_spec)
        if demand_sum > problem_spec.depots[depot_number].Q:
            return True
        else:
            return False


    ######################################################
    # Determine if a route becomes too long if an additional customer is assigned in a certain position
    def routeDurationLimitExceeded(self, route, vehicle_nr, new_customer, new_customer_position, problem_spec):
        new_route = list(route)
        new_route.insert(new_customer_position, new_customer)

        duration = self.routeDuration(new_route, vehicle_nr, problem_spec)
        max_route_duration = problem_spec.depots[self.getDepotNumber(vehicle_nr, problem_spec)].D
        if duration > max_route_duration and max_route_duration != 0:
            return True
        else:
            return False

    ######################################################
    # Plots all the vehicle routes of the chromosome in a single plot
    def visualizeGenes(self, problem_spec):
        fig, ax = plt.subplots(1, 1)
        ax.set_aspect('equal')

        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        background = fig.canvas.copy_from_bbox(ax.bbox)

        for vehicle_nr, route in enumerate(self.vehicle_routes):
            route_x_coords = np.zeros(len(route) + 2)
            route_y_coords = np.zeros(len(route) + 2)
            depot_nr = self.getDepotNumber(vehicle_nr, problem_spec)
            route_x_coords[0] = problem_spec.depots[depot_nr].x
            route_y_coords[0] = problem_spec.depots[depot_nr].y
            for customer_num in range(0, len(route)):
                route_x_coords[customer_num + 1] = route[customer_num].x
                route_y_coords[customer_num + 1] = route[customer_num].y
            route_x_coords[-1] = route_x_coords[0]
            route_y_coords[-1] = route_y_coords[0]
            ax.plot(route_x_coords,route_y_coords, 'x-', color=colors[int(depot_nr)][ vehicle_nr % problem_spec.max_vehicles_per_depot], linewidth=0.8)

        plt.show(block=True)
        self.background = background
        self.fig = fig
        self.ax = ax

    ######################################################
    #
    def printSolutionData(self, problem_spec):
        vehicle_num_from_depot = 1
        prev_depot_nr = 1
        print(self.duration)
        for vehicle_nr, route  in enumerate(self.vehicle_routes):
            depot_nr = self.getDepotNumber(vehicle_nr, problem_spec) + 1
            if depot_nr != prev_depot_nr:
                vehicle_num_from_depot = 1

            if route:
                route_duration = self.routeDuration(route, vehicle_nr, problem_spec)
                route_demand = self.routeDemand(route)
                #customers =
                print(str(depot_nr) + '\t' + str(vehicle_num_from_depot) + '\t' + str(route_duration) + '\t' + str (route_demand) +'\t' + str([0]+ [c.i for c in route] +[0]))
                vehicle_num_from_depot += 1

            prev_depot_nr = self.getDepotNumber(vehicle_nr, problem_spec)



    ###########################################
    # Get depot number of the specified vehicle
    def getDepotNumber(self, vehicle_nr, problem_spec):
        return math.floor(vehicle_nr / problem_spec.max_vehicles_per_depot)



    def updateFitness(self, problem_spec):
        self.fitness = self.duration + self.duration_ol + self.demand_ol

    ######################################################
    # Finds the duration of a single route
    def routeDuration(self, route, vehicle_nr, problem_spec):
        duration = 0
        # The depot is considered the very first "customer"
        depot_nr = self.getDepotNumber(vehicle_nr, problem_spec)
        prev_customer_index = depot_nr + problem_spec.num_customers
        for customer in route:
            customer_index = customer.i - 1
            # Find distance between the current customer and the previous
            duration += problem_spec.cost_matrix[prev_customer_index, customer_index]
            prev_customer_index = customer_index
        # Remember to add the distance from last customer to depot
        duration += problem_spec.cost_matrix[prev_customer_index, depot_nr + problem_spec.num_customers]

        return duration


    ######################################################
    # Finds the demand/load of a single route
    def routeDemand(self, vehicle_route):
        demand = 0
        for customer in vehicle_route:
            demand += customer.q
        return demand

    def tooManyCustomers(self, problem_spec):
        flattened = list(itertools.chain.from_iterable(self.vehicle_routes))

        for i, c in enumerate(flattened):
            flattened[i] = c.i
        customers = len(flattened)
        unique_customers = np.unique(flattened).size

        if customers > problem_spec.num_customers:
            print("\nERROR:\tToo many customers in routes: \t%d customers\n" % customers)
        elif customers < problem_spec.num_customers:
            print("\nERROR:\tToo few customers in routes: \t%d customers\n" % customers)
        if customers > unique_customers:
            print("\nERROR:\tThere are repeated customers in the routes\n")




ga = GA(fileName = 'p01', population_size = 300, generations = 100,
        elite_ratio = 0.01, tournament_ratio = 0.08,
        crossover_prob = 0.4, intra_mutation_prob = 0.10, inter_mutation_prob = 0.10,
        inter_mutation_attempt_rate = 10)
#ga.initializePopulation()
#ga.population[3].tooManyCustomers(ga.problem_spec)
ga.evolutionCycle()
#ga.initializePopulation()
#A = ga.population[0]
#B = copy.deepcopy(ga.population[0])
#B.problem_spec = None

#ga.population[0].vizualizeGenes()
#ga.population[1].vizualizeGenes()
#plt.close("all")

yo = 5