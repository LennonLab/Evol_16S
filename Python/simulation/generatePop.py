from __future__ import division
import numpy as np
import itertools
import math
import argparse
from collections import Counter
import random
#from collections import defaultdict


def generate_base_haplotype(seq_length):
    base_haplotype = ''.join(["A" for i in range(seq_length)])
    return base_haplotype

def random_dna_sequence(seq_length):
    return ''.join(np.random.choice(list('ACTG')) for _ in range(seq_length))

def generate_pop(pop_size, seq_length):
    pop = {}
    base_haplotype = generate_base_haplotype(seq_length)
    pop[base_haplotype] = pop_size
    return pop

def choose_by_weight(weights):
    weights = map(float, weights)
    rndm = random.random() * sum(weights)
    for i, j in enumerate(weights):
        rndm -= j
        if rndm < 0:
            return i

# mutation
def get_mutation_count(mutation_rate, pop_size, seq_length):
    mean = mutation_rate * pop_size * seq_length
    return np.random.poisson(mean)

def get_random_haplotype(pop):
    # Need random haplotype for active only
    haplotypes = pop.keys()
    active_size_step = sum(pop.values())
    frequencies = [x/float(active_size_step) for x in pop.values()]
    total = sum(frequencies)
    frequencies = [x / total for x in frequencies]
    return np.random.choice(haplotypes, p=frequencies)

def get_mutant(haplotype, seq_length, rates):
    alphabet = ['A', 'T', 'G', 'C']
    if rates == 0:
        site = np.random.randint(seq_length)
    else:
        site = np.random.choice(seq_length, p = rates)
    possible_mutations = list(alphabet)
    possible_mutations.remove(haplotype[site])
    mutation = np.random.choice(possible_mutations)
    new_haplotype = haplotype[:site] + mutation + haplotype[site+1:]
    return new_haplotype

def mutation_event(pop, seq_length, rates):
    haplotype = get_random_haplotype(pop)
    if pop[haplotype] > 1:
        pop[haplotype] -= 1
        new_haplotype = get_mutant(haplotype, seq_length, rates)
        if new_haplotype in pop:
            pop[new_haplotype] += 1
        else:
            pop[new_haplotype] = 1
    else:
        new_haplotype = get_mutant(haplotype, seq_length, rates)
        if new_haplotype in pop:
            pop[new_haplotype] += 1
        else:
            pop[new_haplotype] = 1
        del pop[haplotype]


def mutation_step(pop, mutation_rate, seq_length, pop_size, rates):
    mutation_count = get_mutation_count(mutation_rate, pop_size, seq_length)
    for i in range(mutation_count):
        mutation_event(pop, seq_length, rates)

# reproduce active pop, drift acts here
def get_offspring_counts(pop, pop_size):
    haplotypes = pop.keys()
    frequencies = [x/float(pop_size) for x in pop.values()]
    total = sum(frequencies)
    frequencies = [x / total for x in frequencies]
    return list(np.random.multinomial(pop_size, frequencies))

def offspring_step(pop, pop_size):
    counts = get_offspring_counts(pop, pop_size)
    for (haplotype, count) in zip(pop.keys(), counts):
        if (count > 0):
            pop[haplotype] = count
        else:
            del pop[haplotype]

def time_step(pop, mutation_rate, seq_length, pop_size, rates = 0):
    if mutation_rate != 0:
        mutation_step(pop, mutation_rate, seq_length, pop_size, rates)
    #offspring_step(pop, pop_size)


def simulate(pop, mutation_rate, generations, seq_length, pop_size, rates):
    rates = [(x / seq_length) for x in rates]
    for i in range(generations):
        time_step(pop, mutation_rate, seq_length, pop_size, rates)
    return pop

def simulate_save_history(OUT, pop, mutation_rate, generations, seq_length, pop_size):
    for i in range(generations):
        time_step(pop, mutation_rate, seq_length, pop_size)
        return pop

def simulate_history(pop, mutation_rate, generations, seq_length, pop_size):
    history = []
    history.append(pop)
    for i in range(generations):
        print active_size, dormant_size, c, i
        time_step(pop, mutation_rate, seq_length, active_size, dormant_size, c)
        history.append(pop)
    return history
