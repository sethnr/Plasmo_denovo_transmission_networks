#!/bin/python

import sys
#sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h


# specify parameters
L = 10000  #no_loci
#C=int(10**8)    #capacity
C=int(10**4)    #capacity
print str(C)
N=1       #init pop

# set up population
pop = h.haploid_highd(L)                        # produce an instance of haploid_highd with L loci
pop.carrying_capacity = C                  # set the average population size to 50000
pop.outcrossing_rate = 0                        # remove outcrossing
pop.crossover_rate = 0                          # set the crossover rate to zero (core genome)
pop.mutation_rate = 0.01 / pop.carrying_capacity # per locus mutation rate equal to 0.1/N

# set fitness landscape
selection_coefficients = 0.0*np.ones(pop.L)     # most loci are neutral

#m = 10
#selection_coefficients[::m] = -0.1             # every m-th locus is strongly deleterious
pop.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness

# initialize the population in linkage equilibrium with the specified allele frequencies
initial_allele_frequencies = np.zeros(pop.L)    # set all alleles to zero (clonal founder)

#pop.set_allele_frequencies(initial_allele_frequencies, init_size)
pop.set_wildtype(N)


# evolve for 2000 generations and track the allele frequencies
maxgen = 10000
allele_frequencies = [pop.get_allele_frequencies()]
tp = [pop.generation]
while pop.generation < maxgen:
#    print "generation:", pop.generation,pop.get_allele_frequencies()
    
    pop.evolve(1)

    # save allele frequencies and time
    allele_frequencies.append(pop.get_allele_frequencies()) 
    tp.append(pop.generation)
    # every 200 generations, make one of the deleterious mutations beneficial
    if (pop.generation % 50 == 0):
        print "generation:", pop.generation, ':',
        af = pop.get_allele_frequencies()
#        afInt = (af*1000).astype(int)
        (hist,edges) = np.histogram(af,bins=10,range=(0,1))
        print pop.N," individuals ", sum(af > 0), " mutations",hist
        # update fitness function
#        selection_coefficients[m*np.random.randint(0,25)] = 0.01
#        pop.set_trait_additive(selection_coefficients)

# convert to an array to enable slicing
allele_frequencies = np.array(allele_frequencies)

print >>sys.stderr, allele_frequencies
