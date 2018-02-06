import problem_reader as pr
from genotype import *
import random
import math
from timeit import default_timer as timer
from operator import itemgetter, attrgetter
import copy
import matplotlib.pyplot as plt
import numpy as np




class GA:
    def __init__(self, fileName, population_size, generation_size, generations, elite_ratio, tournament_ratio, div_bound, crossover_prob, intra_mutation_prob, inter_mutation_prob,
                 inter_mutation_attempt_rate):
        self.problem_spec = pr.ProblemSpec(fileName)
        self.early_stopping_percents = [30, 20, 10, 5]

        # Parameters
        self.min_population_size = population_size
        self.generation_size = generation_size      # > population size
        self.generations = generations
        self.elite_ratio = elite_ratio
        self.tournament_ratio = tournament_ratio
        self.div_bound = div_bound
        self.crossover_prob = crossover_prob
        self.intra_mutation_prob = intra_mutation_prob
        self.inter_mutation_prob = inter_mutation_prob
        self.inter_muation_attempt_rate = inter_mutation_attempt_rate

        # Data storage
        self.population = []
        self.average_fitness_history = []
        self.best_fitness_history = []


    ######################################################
    # Create an initial population by creating Genotype()'s and storing them in self.population
    def initializePopulation(self):
        for _ in range(0, self.min_population_size):
            genotype = Genotype()
            genotype.initGenes(self.problem_spec)
            genotype.tooManyCustomers(self.problem_spec)
            self.population.append(genotype)


    ######################################################
    # Applies a change to one random vehicle route (inversion, single customer re-routing, swapping)
    def intraMutation(self, offspring):
        random.seed()
        type = random.randint(1,3)
        if type == 1: # inversion with one vehicle route
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
        if type == 2: # single customer re-route from one route to another within the entire chromosome
            num_vehicles = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
            random_vehicle_route = None
            # Acquire a non-empty vehicle route
            while not random_vehicle_route:
                random_vehicle_route_nr = random.randint(0, num_vehicles - 1)
                random_vehicle_route = offspring.vehicle_routes[random_vehicle_route_nr]
            random_customer_nr = random.randint(0, len(random_vehicle_route) - 1)
            customer = random_vehicle_route[random_customer_nr]
            random_vehicle_route.pop(random_customer_nr)

            # Find the best feasible location and insert it into the route
            cost = self.insertionCost(customer, offspring)
            self.insertCustomerInRoute(customer, offspring, cost, 1)

        if type == 3: # swapping of two random customers within one particular depot
            depot_nr = random.randint(0, self.problem_spec.num_depots - 1)
            vehicle_start_index = depot_nr * self.problem_spec.max_vehicles_per_depot
            max_num_vehicles_per_depot = self.problem_spec.max_vehicles_per_depot

            customers_nrs_to_swap = []
            routes = []
            for _ in range (0,2):
                route = None
                # Acquire a random, non-empty vehicle route
                while not route:
                    route_nr = random.randint(vehicle_start_index, vehicle_start_index + max_num_vehicles_per_depot - 1)
                    route = offspring.vehicle_routes[route_nr]
                customer_nr = random.randint(0, len(route) - 1)
                # Acquire a random customer from said route
                customers_nrs_to_swap.append(customer_nr)
                routes.append(route)
            temp_customer = routes[0][customers_nrs_to_swap[0]]
            routes[0][customers_nrs_to_swap[0]] = routes[1][customers_nrs_to_swap[1]]
            routes[1][customers_nrs_to_swap[1]] = temp_customer
        offspring.tooManyCustomers(self.problem_spec)


    ######################################################
    # Applies a swapping of two random customers within the entire chromosome
    def interMutation(self, offspring):
        total_num_vehicle_routes = self.problem_spec.max_vehicles_per_depot * self.problem_spec.num_depots
        customers_nrs_to_swap = []
        routes = []
        for _ in range(0, 2):
            route = None
            # Acquire a random, non-empty vehicle route
            while not route:
                route_nr = random.randint(0, total_num_vehicle_routes - 1)
                route = offspring.vehicle_routes[route_nr]
            customer_nr = random.randint(0, len(route) - 1)
            # Acquire a random customer from said route
            customers_nrs_to_swap.append(customer_nr)
            routes.append(route)
        temp_customer = routes[0][customers_nrs_to_swap[0]]
        routes[0][customers_nrs_to_swap[0]] = routes[1][customers_nrs_to_swap[1]]
        routes[1][customers_nrs_to_swap[1]] = temp_customer


    ######################################################
    # Perform Best Cost Route Crossover on two parents. Produces two offspring
    def crossover(self, parent1, parent2):
        random.seed()

        # Initialize offspring
        offspring1 = copy.deepcopy(parent1)
        offspring2 = copy.deepcopy(parent2)

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
            self.insertCustomerInRoute(route1_customer, offspring2, cost, 0.8)

        # Find location and insert customers from route 2 into offspring 1
        for route2_customer in route2:
            cost = self.insertionCost(route2_customer, offspring1)
            self.insertCustomerInRoute(route2_customer, offspring1, cost, 0.8)


        offspring1.tooManyCustomers(self.problem_spec)
        offspring2.tooManyCustomers(self.problem_spec)


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
    def insertCustomerInRoute(self, new_customer, individual, insertion_cost, p_best_loc):
        k = random.random()

        # Choose best insertion location
        #TODO: Handle when there is no feasible insertions
        if k <= p_best_loc:
            inserted = False
            for i in range(len(insertion_cost)):
                if insertion_cost[i][2]:
                    insertion_location = insertion_cost[i][0]
                    inserted = True
                    break
            # If no locations were feasible, choose a random one
            if not inserted:
                insertion_location = insertion_cost[random.randint(0, len(insertion_cost) - 1)][0]
        # Choose first entry in list regardless of feasibility
        else:
            randomy_list_entry = random.randint(0, len(insertion_cost) - 1)
            insertion_location = insertion_cost[randomy_list_entry][0]

        # Insert the new customer in the chosen route and position in the individual
        insertion_vehicle, insertion_customer_pos = insertion_location
        individual.vehicle_routes[insertion_vehicle].insert(insertion_customer_pos, new_customer)


    ######################################################
    # Returns the parent with the highest fitness from the parents parameter
    def getBestParent(self, parents):
        fitnesses = []
        for parent in parents:
            fitnesses.append(parent.fitness)
        best_parent = copy.deepcopy(parents[fitnesses.index(max(fitnesses))])
        return best_parent


    ######################################################
    # Return the top self.elites percent of the population
    def getPopulationElite(self, fitnesses):
        population = self.population
        sorted_population = [x for _, x in sorted(zip(fitnesses, population))] #sorts population in descending order based on fitnesses
        elites_num = int(math.ceil(self.min_population_size * self.elite_ratio))
        elites = sorted_population[0:elites_num]
        elites = copy.deepcopy(elites)
        return elites

    ######################################################
    # Return two genotypes from the population based on a tournament selection
    def parentSelection(self):
        num_individuals = int(math.ceil(self.tournament_ratio * len(self.population)))
        if num_individuals < 2:
            num_individuals = 2

        selected = random.sample(range(0, self.min_population_size), num_individuals)
        tournament_population = []
        for i in selected:
            tournament_population.append(self.population[i])
        fitnesses = []
        for i in tournament_population:
            fitnesses.append(i.fitness)
        sorted_tournament_population = [x for _, x in sorted(zip(fitnesses, tournament_population))]  # sorts population in descending order based on fitnesses

        return sorted_tournament_population[0:2]

    ######################################################
    # Successively eliminates clones, then bad individuals, until the population size is at minimum
    def survivorSelection(self):
        for i in range(0, len(self.population) - self.min_population_size):
            clones = self.findAllClonesInPopulation()
            if clones:
                # Remove individual in clones with worst fitness
                fitnesses = []
                for i, _ in clones:
                    fitnesses.append(i.fitness)
                sorted_clones = [x for _, x in sorted(zip(fitnesses, clones))]
                self.population.pop(sorted_clones[-1][1])
            else:
                # Remove individual in whole population with worst fitness
                sorted_population = sorted(self.population, key=attrgetter('fitness'))
                sorted_population.pop()
                self.population = sorted_population

    def findAllClonesInPopulation(self):
        clones = []
        clone_indices = []
        for i, individual1 in enumerate(self.population):
            if i in clone_indices:
                continue
            for j, individual2 in enumerate(self.population):
                if j in clone_indices:
                    continue
                elif individual1.isClone(individual2, self.problem_spec) and i != j:
                    clones.append((individual1, i))
                    clones.append((individual2, j))
                    clone_indices += [i, j]
        return clones


    ######################################################
    # Eliminates all but the 1/3 best individuals, and creates new ones as in the initialization phase
    def diversify(self, fitnesses):
        population = self.population
        sorted_population = [x for _, x in sorted(zip(fitnesses, population))]  # sorts population in ascending order based on fitnesses

        # Save the best individuals
        diversification_num = int(math.ceil(len(population) / 3))
        best_individuals = sorted_population[:diversification_num]
        new_generation = copy.deepcopy(best_individuals)

        # Generate new individuals: 2/3
        generation_num = len(population) - diversification_num
        for _ in range(0, generation_num):
            genotype = Genotype()
            genotype.initGenes(self.problem_spec)
            genotype.tooManyCustomers(self.problem_spec)
            new_generation.append(genotype)
        return new_generation


    def plotHistoryGraph(self, yvals, xvals=None, xtitle='X', ytitle='Y', title='Y = F(X)', label=''):
        xvals = xvals if xvals is not None else list(range(len(yvals)))
        plt.plot(xvals, yvals, label=label)
        plt.pause(0.01)
        plt.xlabel(xtitle); plt.ylabel(ytitle); plt.title(title)
        plt.legend()
        plt.draw()
        plt.pause(.001)

    ######################################################
    # Evolves an initial population by means of mutation and recombinations over a given number of generations
    def evolutionCycle(self):
        start_evolution_time = timer()
        random.seed(None)

        # Generate the initial population
        self.initializePopulation()
        prev_best = None
        it_div = 0

        # Start the evolution process
        for generation in range(1, self.generations + 1):
            # Evaluate the fitness of all individuals in the population and compute the average
            population_fitness = []
            for individual in self.population:
                population_fitness.append(individual.fitness)

            # Store fitness data
            average_fitness = sum(population_fitness) / self.min_population_size
            best_fitness = min(population_fitness)
            best = np.argmin(population_fitness)
            self.average_fitness_history.append(average_fitness)
            self.best_fitness_history.append(best_fitness)

            # Create a new population and directly copy the elites of the previous generation to the new one
            new_population = []
            elites = self.getPopulationElite(population_fitness)
            new_population += elites

            # Print generation data to terminal
            print("Generation: %d" % generation)
            print("Best duration: %f " % self.population[int(best)].duration)
            print("Best fitness score: %f (infeasibility_count: %d)" % (best_fitness, elites[0].infeasibility_count))
            print("Average fitness score: %f\n" % average_fitness)

            # Check for early stop
            for index, percent in enumerate(self.early_stopping_percents):
                if elites[0].checkForSatisfyingSolution(self.problem_spec, percent):
                    self.early_stopping_percents.pop(index)
                    ask_user_to_continue =  input("Feasible solution found within "+str(percent)+"% from optimal solution. "
                                                  "Continue? [y/n]: ")
                    break
                else:
                    ask_user_to_continue = 'y'
                    break
            if ask_user_to_continue == 'n':
                break

            # Diversification process
            current_best = elites[0]
            if prev_best:
                if current_best.fitness == prev_best.fitness: # no improvement in best solution
                    it_div += 1
                else:
                    it_div = 0
            if it_div >= self.div_bound:   # no improvement in best solution for it_div_bound iterations
                # Populate the new population through diversification
                new_population = self.diversify(fitnesses=population_fitness)
                it_div = 0

            else:
                # Populate the new population by doing recombination, mutation and survivors
                while len(new_population) < self.generation_size:
                    # Do tournament selection to choose parents
                    parents = self.parentSelection()

                    # Roll dice on what event shall happen (recombination or survivors)
                    event = random.random()
                    if event < self.crossover_prob:
                        offsprings = self.crossover(parents[0], parents[1])

                    else:
                        offsprings = copy.deepcopy(parents)

                    for offspring in offsprings:
                        # Roll dice on what event shall happen (inter-mutation, intra-mutation or no mutation)
                        event = random.random()
                        if event < self.inter_mutation_prob and generation % self.inter_muation_attempt_rate == 0:
                            self.interMutation(offspring)
                        elif self.inter_mutation_prob < event and event < self.inter_mutation_prob + self.intra_mutation_prob \
                                and not generation % self.inter_muation_attempt_rate == 0:
                            self.intraMutation(offspring)

                        offspring.updateFitnessVariables(self.problem_spec)
                        offspring.updateFitness(self.problem_spec)
                        new_population += [offspring]
            prev_best = current_best

            # Stop population growing in case everything is not even numbers
            if len(new_population) > (self.generation_size):
                new_population.pop()

            self.population = new_population
            self.survivorSelection()

        end_evolution_time = timer()
        print("Algorithm terminated - elapsed time: %f" % (end_evolution_time-start_evolution_time))
        bestGenotype = elites[0]
        self.plotHistoryGraph(yvals = self.best_fitness_history, xtitle = "Generations", ytitle = "Fitness",
                              title = "Best individual fitness history")
        bestGenotype.printGenotypeData(self.problem_spec)
        bestGenotype.visualizeGenes(self.problem_spec)




if __name__ == '__main__':
    ga = GA(fileName='p01', population_size=25, generation_size=25, generations=10000,
            elite_ratio=0.4, tournament_ratio=0.25, div_bound=200,
            crossover_prob=0.6, intra_mutation_prob=0.2, inter_mutation_prob=0.25,
            inter_mutation_attempt_rate=10)

    ga.evolutionCycle()
    y = 0