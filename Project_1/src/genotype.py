import random
import itertools
import numpy as np
import matplotlib.pyplot as plt

depot1_colors = ['crimson', 'brown', 'coral', 'red', 'deeppink', 'tomato', 'darkred']
depot2_colors = ['grey', 'blue', 'darkblue', 'teal', 'black', 'royalblue', 'indigo']
depot3_colors = ['green', 'darkgreen', 'seagreen', 'lime', 'yellowgreen', 'olivedrab', 'darkslategrey']
depot4_colors = ['chartreuse', 'gold', 'yellow', 'olive', 'darkorange', 'orange', 'yellowgreen']

colors = []
colors.append(depot1_colors)
colors.append(depot2_colors)
colors.append(depot3_colors)
colors.append(depot4_colors)

class Genotype:
    def __init__(self):
        self.fitness = float("Inf")
        self.vehicle_routes = []

        self.demand_ol = 0          # demand overload
        self.duration_ol = 0        # duration_overload
        self.duration = 0

        self.infeasibility_count = 0    # number of infeasible insertions (occur in crossover and init)

    def __lt__(self, other):
        return self.fitness < other.fitness

    def __gt__(self, other):
        return self.fitness > other.fitness

    def __eq__(self, other):
        return (self.fitness == other.fitness)

    def __deepcopy__(self, memodict={}):
        copy = Genotype()
        copy.vehicle_routes = [list(r) for r in self.vehicle_routes]
        copy.fitness = self.fitness
        copy.infeasibility_count = self.infeasibility_count
        copy.demand_ol = self.demand_ol
        copy.duration_ol = self.duration_ol
        copy.duration = self.duration
        return copy

    ######################################################
    # A clone of an individual has the same cost and the same depot assignments
    def isClone(self, other, problem_spec):
        if self.fitness == other.fitness and self.infeasibility_count == other.infeasibility_count:
            return True
        else:
            return False

    ######################################################
    #
    def initGenes(self,  problem_spec):
        # Initialization of placement of customers to vehicle routes
        max_vehicles = problem_spec.max_vehicles_per_depot * problem_spec.num_depots
        self.vehicle_routes = [[] for i in range(0, max_vehicles)]

        self.assignCustomers(problem_spec)
        self.updateFitnessVariables(problem_spec)
        self.updateFitness(problem_spec)

    ######################################################
    #
    def assignCustomers(self, problem_spec):
        random.seed()
        range_num_customers = list(range(0, len(problem_spec.customers)))
        customers_placed = 0

        while range_num_customers:
            customer_nr = random.choice(range_num_customers)#.pop(random.randint(0, customers_left))
            range_num_customers.remove(customer_nr)
            customer = problem_spec.customers[customer_nr]

            depot = customer.candidate_list[0]
            vehicle_start_index = depot * problem_spec.max_vehicles_per_depot
            inserted = False
            # Proceed by putting the customer in the first available vehicle
            for vehicle_nr in range(vehicle_start_index, vehicle_start_index + problem_spec.max_vehicles_per_depot):
                # Only append if it doesnt cause vehicle overload :
                if ((not self.vehicleOverloaded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer.q, problem_spec))
                    and (not self.routeDurationLimitExceeded(self.vehicle_routes[vehicle_nr], vehicle_nr, customer,
                                                             len(self.vehicle_routes[vehicle_nr]), problem_spec))):
                    self.vehicle_routes[vehicle_nr].append(customer)
                    inserted = True
                    break
            # If no feasible position found, insert in random route in closest depot
            if not inserted:
                vehicle_start_index = depot * problem_spec.max_vehicles_per_depot
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
        self.infeasibility_count = 0
        for vehicle_nr, route in enumerate(self.vehicle_routes):
            depot_num = problem_spec.vehicle_to_depot_dict[vehicle_nr]

            max_duration = problem_spec.depots[depot_num].D
            max_load = problem_spec.depots[depot_num].Q
            route_duration = self.routeDuration(route, vehicle_nr, problem_spec)
            route_demand = self.routeDemand(route)

            if route_duration > max_duration:
                self.duration_ol += route_duration - max_duration
                self.infeasibility_count += (route_duration) / max_duration
            if route_demand > max_load:
                self.demand_ol += route_demand - max_load
                self.infeasibility_count += (route_demand) / max_load
            self.duration += route_duration


    ######################################################
    # Determine if a vehicle becomes overloaded if its assigned an additional customer, returns True/False
    def vehicleOverloaded(self, vehicle_route, vehicle_nr, customer_demand, problem_spec):
        demand_sum = self.routeDemand(vehicle_route)
        demand_sum += customer_demand
        depot_number = problem_spec.vehicle_to_depot_dict[vehicle_nr]
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
        max_route_duration = problem_spec.depots[problem_spec.vehicle_to_depot_dict[vehicle_nr]].D
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

        #background = fig.canvas.copy_from_bbox(ax.bbox)
        prev_depot_nr = -1
        for vehicle_nr, route in enumerate(self.vehicle_routes):
            route_x_coords = np.zeros(len(route) + 2)
            route_y_coords = np.zeros(len(route) + 2)
            depot_nr = problem_spec.vehicle_to_depot_dict[vehicle_nr]
            route_x_coords[0] = problem_spec.depots[depot_nr].x
            route_y_coords[0] = problem_spec.depots[depot_nr].y
            for customer_num in range(0, len(route)):
                route_x_coords[customer_num + 1] = route[customer_num].x
                route_y_coords[customer_num + 1] = route[customer_num].y
            route_x_coords[-1] = route_x_coords[0]
            route_y_coords[-1] = route_y_coords[0]
            if depot_nr != prev_depot_nr:
                ax.plot(route_x_coords[0], route_y_coords[0], 's', color='black', fillstyle='none', markersize=7, linewidth=1)
            ax.plot(route_x_coords[1:-1], route_y_coords[1:-1], 'x', color=colors[int(depot_nr) % 4][ vehicle_nr % problem_spec.max_vehicles_per_depot % 4], linewidth=1)
            ax.plot(route_x_coords,route_y_coords, '-', color=colors[int(depot_nr) % 4][ vehicle_nr % problem_spec.max_vehicles_per_depot % 4], linewidth=2)
            prev_depot_nr = depot_nr
        plt.show(block=True)

    ######################################################
    # Print the solution (routes), demand over routes and total durations of all routes
    def printGenotypeData(self, problem_spec):
        vehicle_num_from_depot = 1
        prev_depot_nr = 1
        percent_from_optimal = 100*(self.duration / problem_spec.solution_cost)
        print('This solution cost: %.2f' % self.duration)
        print('Optimal solution cost: %.2f' % problem_spec.solution_cost)
        print('Percent within the optimal solution: %.2f%%\n' % (percent_from_optimal - 100))
        print(self.duration)
        for vehicle_nr, route  in enumerate(self.vehicle_routes):
            depot_nr = problem_spec.vehicle_to_depot_dict[vehicle_nr] + 1
            if depot_nr != prev_depot_nr:
                vehicle_num_from_depot = 1
            if route:
                route_duration = self.routeDuration(route, vehicle_nr, problem_spec)
                route_demand = self.routeDemand(route)
                full_route = (['0'] + [str(c.i) for c in route] + ['0'])
                print('%d\t%d\t%.2f\t%d\t' % (depot_nr, vehicle_num_from_depot, route_duration, route_demand), end='')
                print(' '.join((full_route)))
                vehicle_num_from_depot += 1
            prev_depot_nr = depot_nr

    ######################################################
    #
    def updateFitness(self, problem_spec):
        self.fitness = self.duration + self.duration_ol + self.demand_ol + 2*problem_spec.max_cost*(self.infeasibility_count)

    ######################################################
    # Finds the duration of a single route
    def routeDuration(self, route, vehicle_nr, problem_spec):
        duration = 0
        # The depot is considered the very first "customer"
        depot_nr = problem_spec.vehicle_to_depot_dict[vehicle_nr]
        prev_customer_index = depot_nr + problem_spec.num_customers
        for customer in route:
            customer_index = customer.i - 1
            # Find distance between the current customer and the previous
            duration += problem_spec.cost_matrix[prev_customer_index, customer_index]
            duration += customer.d
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

    ######################################################
    #
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
        elif customers < unique_customers:
            print("\nERROR:\tThere are unique customers missing in the routes\n")

    ######################################################
    #
    def checkForSatisfyingSolution(self, problem_spec, percent):
        percent_from_optimal = (100 * (self.duration / problem_spec.solution_cost)) - 100

        if percent_from_optimal <= percent and not self.infeasibility_count:
            return True
        else:
            return False
