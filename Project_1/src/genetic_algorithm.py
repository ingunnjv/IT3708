import problem_reader as pr

class GA:
    def __init__(self, fileName):
        self.mdvrp = pr.ProblemSpec()
        self.mdvrp.readProblemFile(fileName)
        self.population_size = 0

    def initialise_population(self):
        pass

    def fitness(self):
        pass

    def parent_selection(self):
        pass

    def mutation(self, parent):
        pass

    def crossover(self, parent1, parent2):
        pass

    def survivor_selection(self):
        pass

ga = GA('p01')