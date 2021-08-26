from importing import *

#functions import
import PARSEC_functions
import optimization_NSGA2 as optim
import optimization_GA as optimGA

from pymoo.factory import get_problem
from pymoo.util.plotting import plot


absolute_limits = [[0, 5], [0, 3]]

n_ind = 50
n_gen = 25
mut_rate = 0.1
t_size = 2

def evaluate(chrom):
    x1, x2 = chrom[0], chrom[1]
    f1 = 4*(x1**2) + 4*(x2**2)
    f2 = (x1 - 5)**2 + (x2 - 5)**2
    return [f1, f2]

def is_valid(chrom):
    #return True
    x1, x2 = chrom[0], chrom[1]
    C1 = (x1 - 5)**2 + x2**2 #≤ 25
    C2 = (x1 - 8)**2 + (x2 + 3)**2 #≥ 7.7
    if (C1 > 25) or C2 < 7.7:
        return False
    return True

optimized_pop = optim.run_GA(n_gen, n_ind, absolute_limits, evaluate, mut_rate = mut_rate, t_size = t_size, valid_func = is_valid) 

for individual in optimized_pop:
    funcs = evaluate(individual.get_chrom())
    print(f'chrom: {individual.get_chrom()} Front: {individual.get_Front()} Crowding: {individual.get_Crowding()} Dominated_counter: {individual.get_Dominated_counter()}')
    plt.plot(funcs[0], funcs[1], 'bo')
problem = get_problem("bnh")
plt.title("População final")
plot(problem.pareto_front(), no_fill = True)
plt.show()