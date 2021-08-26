from importing import *
import Auxiliary

class Individual:
    def __init__(self, chrom = [], front = -1, crowd = 0, domcount = infinite, ObjVal = [], valid = True):
        self.chrom = chrom
        self.front = front
        self.crowd = crowd
        self.domcount = domcount
        self.ObjVal = ObjVal
        self.valid = valid

    def set_valid(self, boolean):
        self.valid = boolean
    def get_valid(self):
        return self.valid

    def set_chrom(self, chrom):
        self.chrom = chrom
    def get_chrom(self):
        return self.chrom

    def set_Front(self, front):
        self.front = front
    def get_Front(self):
        return self.front

    def set_Crowding(self, crowd):
        self.crowd = crowd
    def get_Crowding(self):
        return self.crowd

    def set_Dominated_counter(self, counter):
        self.domcount = counter
    def get_Dominated_counter(self):
        return self.domcount

    def set_ObjVal(self, value, index = -1):
        if index == -1:
            self.ObjVal = value
        else:
            self.ObjVal[index] = value
    def get_ObjVal(self, index = -1):
        if index == -1:
            return self.ObjVal
        return self.ObjVal[index]




#main functions
def run_GA(n_gen, n_ind, gene_limits, evaluate_func, valid_func = Auxiliary.return_true, mut_rate = 0.05, t_size = 2):
    print("Geração 1")
    pop_ini = create_first_gen(n_ind, evaluate_func, valid_func, gene_limits)
    current_pop = pop_ini[:]
    for generation in range(2, n_gen + 1): 
        print(f"Geração {generation}")
        for ind in current_pop:
            #print(ind.get_chrom())
            funcs = evaluate_func(ind.get_chrom())
            plt.plot(funcs[0], funcs[1], 'ro')
        offspring = create_offspring(current_pop, gene_limits, evaluate_func, valid_func, mut_rate = mut_rate, t_size = t_size)
        current_pop = reinsert(current_pop, offspring)
        plt.show()
    #plt.title("Todos os indivíduos")
    return current_pop

def create_first_gen(n_ind, evaluate_func, valid_func, limits):
    individuals = 0 
    chrom_pop = []
    first_gen = []
    while individuals < n_ind:
        chromossome = generate_chrom(limits)
        if valid_func(chromossome) and chromossome not in chrom_pop:
            chrom_pop.append(chromossome)
            individuals += 1
    for chrom in chrom_pop:
        first_gen.append(Individual(chrom = chrom))
    evaluate_population(first_gen, evaluate_func, valid_func)
    first_gen = assign_fronts(first_gen)
    first_gen = assign_crowding(first_gen)
    return first_gen
    
def create_offspring(population, limits, evaluate_func, valid_func, mut_rate = 0.05, t_size = 2):
    init_chrom_pop = []
    chrom_pop = []
    length = len(population)
    for individual in population:
        init_chrom_pop.append(individual.get_chrom())
    population = strip_rejected(population)
    population = strip_equal(population)
    while len(chrom_pop) < length:
        mother = select_tournament(population, tournament_size = t_size)
        father = select_tournament(population, tournament_size = t_size)
        son_chrom = crossoverFull(mother.get_chrom(), father.get_chrom())
        son_chrom = mutate(son_chrom, mut_rate, limits)
        if valid_func(son_chrom) and (son_chrom not in init_chrom_pop): #and (son_chrom not in chrom_pop):
            chrom_pop.append(son_chrom)
    offspring = []
    for chrom in chrom_pop:
        offspring.append(Individual(chrom = chrom))
    evaluate_population(offspring, evaluate_func, valid_func)
    return offspring

def evaluate_population(population, evaluate_func, valid_func):
    for individual in population:
        valid = valid_func(individual.get_chrom())
        if not valid:
            individual.set_valid(False)
            continue
        obj_vals = evaluate_func(individual.get_chrom())
        individual.set_ObjVal(obj_vals)




#Reinsertion
def reinsert(pop_ini, pop_new):
    length = len(pop_ini)
    pop_final = []
    pop_ini.extend(pop_new)
    pop_ini = strip_rejected(pop_ini)
    pop_ini = strip_equal(pop_ini)
    pop_ini = assign_fronts(pop_ini)
    pop_ini = assign_crowding(pop_ini)
    while len(pop_final) < length:
        best = pop_ini[0]
        i = 1
        hold = 0
        while i < len(pop_ini):
            if pop_ini[i].get_Front() < best.get_Front():
                best = pop_ini[i]
                hold = i
            elif pop_ini[i].get_Front() == best.get_Front():
                if pop_ini[i].get_Crowding() > best.get_Crowding():
                    best = pop_ini[i]
                    hold = i
            i += 1
        pop_final.append(best)
        pop_ini.pop(hold)
    return pop_final



#Crossover
def crossover1Point(mother, father):
    crossover_point = random.randint(0, len(mother) - 1)
    son1_genes = mother[:crossover_point]
    son1_genes.extend(father[crossover_point:])
    son2_genes = father[:crossover_point]
    son2_genes.extend(mother[crossover_point:])
    return son1_genes, son2_genes

def crossover2Point(mother, father):
    crossover_point1 = random.randint(0, len(mother) - 1)
    crossover_point2 = random.randint(0, len(mother) - 1)
    if (crossover_point2 == crossover_point1) or (abs(crossover_point1 - crossover_point2) == len(mother) - 1):
        return mother, father
    if crossover_point2 < crossover_point1:
        temp = crossover_point1
        crossover_point1 = crossover_point2
        crossover_point2 = temp
    son1_genes = mother[:crossover_point1]
    son1_genes.extend(father[crossover_point1:crossover_point2])
    son1_genes.extend(mother[crossover_point2:])
    son2_genes = father[:crossover_point1]
    son2_genes.extend(mother[crossover_point1:crossover_point2])
    son2_genes.extend(father[crossover_point2:])
    return son1_genes, son2_genes

def crossoverFull(mother, father):
    son_genes = []
    length = len(mother)
    for i in range(length):
        chosen = random.random()
        if chosen > 0.5:
            son_genes.append(mother[i])
        else:
            son_genes.append(father[i])
    return son_genes



#Mutation
def mutate(chrom, mut_rate, limits):
    length = len(chrom)
    for i in range(length):
        mut = random.random()
        if mut < mut_rate:
            chrom[i] = Auxiliary.generate_random_parameter(limits[i])
    return chrom



#Selection
def select_tournament(population, tournament_size = 2):
    index_list = []
    for _ in range(tournament_size):
        rand_index = random.randint(0, len(population) - 1)
        while rand_index in index_list:
            rand_index = random.randint(0, len(population) - 1)
        index_list.append(rand_index)
    best = population[index_list[0]]
    for i in range(1, tournament_size): 
        new = population[index_list[i]]
        best = battle(best, new)
    return best



#Auxiliary
def assign_domination(population, n_objvals):
    if len(population) == 1:
        population[0].set_Dominated_counter(0)
        return
    elif len(population) == 0:
        return
    length = len(population)
    for i in range(length):
        dominated_counter = 0
        for j in range(length):
            if j == i:
                continue
            if calculate_domination(population[j], population[i], n_objvals):
                dominated_counter += 1
        population[i].set_Dominated_counter(dominated_counter)

def assign_fronts(population):
    n_objvals = len(population[0].get_ObjVal())
    copy = population[:]
    new_pop = []
    front = 1
    while len(copy) > 0:
        assign_domination(copy, n_objvals)
        i = 0
        while i < len(copy):
            if copy[i].get_Dominated_counter() == 0:
                copy[i].set_Front(front)
                new_pop.append(copy[i])
                copy.pop(i)
            else:
                i += 1
        front += 1
    return new_pop

def assign_crowding(population):
    n_objvals = len(population[0].get_ObjVal())
    length = len(population)
    reset_crowding(population)
    for i in range(n_objvals):
        copy = sort_ObjVal(population, i)
        copy[0].set_Crowding(infinite)
        copy[-1].set_Crowding(infinite)
        minval = copy[0].get_ObjVal(i)
        maxval = copy[-1].get_ObjVal(i)
        for j in range(1, length - 1):
            crowding = abs((copy[j + 1].get_ObjVal(i) - copy[j - 1].get_ObjVal(i))/(maxval - minval)) #abs(Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j + 1].get_ObjVal(i)) - Auxiliary.linear_interpolate(minval, maxval, 0, 1, copy[j - 1].get_ObjVal(i)))
            copy[i].set_Crowding(copy[i].get_Crowding() + crowding)
    return copy

def reset_crowding(population):
    for ind in population:
        ind.set_Crowding(0)

def battle(indA, indB):
    if indA.get_Front() < indB.get_Front():
        return indA
    elif indA.get_Front() == indB.get_Front():
        if indA.get_Crowding() > indB.get_Crowding():
            return indA
    return indB

def sort_ObjVal(population, objval_index):
    new_pop = []
    copy = population[:]
    while len(copy) > 0:
        min_obj_val = infinite
        i = 0
        hold = 0
        while i < len(copy):
            obj_val = copy[i].get_ObjVal(objval_index)
            if obj_val < min_obj_val:
                min_obj_val = obj_val
                hold = i
            i += 1
        new_pop.append(copy[hold])
        copy.pop(hold)
    return new_pop

def calculate_domination(indA, indB, n_objvals):
    boolA = indA.get_ObjVal(0) <= indB.get_ObjVal(0)
    boolB = indA.get_ObjVal(0) < indB.get_ObjVal(0)
    for i in range(1, n_objvals):
        boolA = boolA and (indA.get_ObjVal(i) <= indB.get_ObjVal(i))
        boolB = boolB or (indA.get_ObjVal(i) < indB.get_ObjVal(i))
    return (boolA and boolB) #if True, A dominates B, if False, A does not dominate B

def strip_rejected(population):
    i = 0
    while i < len(population): 
        if not population[i].get_valid():
            population.pop(i)
            i -= 1
        i += 1
    return population

def strip_equal(population):
    res = []
    for i in population:
        if i not in res:
            res.append(i)
    return res

def generate_chrom(limits):
    chrom = []
    for limit in limits:
        param = Auxiliary.generate_random_parameter(limit)
        chrom.append(param)
    return chrom
