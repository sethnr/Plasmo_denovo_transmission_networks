#!/bin/python

import sys
#sys.path.insert(0, '../pkg/python')

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import FFPopSim as h


# specify parameters
L = 5e5  #no_loci
#C=int(10**8)    #capacity
C=1e6    #capacity
N=1       #init pop

SNPmutRate=4.07e-10
INDELmutRate=2.4e-9

#locus multiplier
m = 23e6/(L*2)


print str(C),"inds",L,"bp","("+str(m)+")"
print "SNPrate:   "+str(SNPmutRate)+" ("+str(SNPmutRate*m)+")"
print "INDELrate: "+str(INDELmutRate)+" ("+str(INDELmutRate*m)+")"

SNPmutRate=SNPmutRate*m
INDELmutRate=INDELmutRate*m

# set up population
pop = h.haploid_highd(L)                        # produce an instance of haploid_highd with L loci
pop.carrying_capacity = C                  # set the average population size to 50000
pop.outcrossing_rate = 0                        # remove outcrossing
pop.crossover_rate = 0                          # set the crossover rate to zero (core genome)
pop.mutation_rate = SNPmutRate                  # per locus mutation rate equal to 0.1/N
selection_coefficients = 0.0*np.ones(pop.L)     # most loci are neutral
pop.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness
initial_allele_frequencies = np.zeros(pop.L)    # set all alleles to zero (clonal founder)
pop.set_wildtype(N)


popI = h.haploid_highd(L)                        # produce an instance of haploid_highd with L loci
popI.carrying_capacity = C                  # set the average population size to 50000
popI.outcrossing_rate = 0                        # remove outcrossing
popI.crossover_rate = 0                          # set the crossover rate to zero (core genome)
popI.mutation_rate = INDELmutRate                # per locus mutation rate ~10X SNP rate
selection_coefficients = 0.0*np.ones(popI.L)     # most loci are neutral
popI.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness
initial_allele_frequencies = np.zeros(popI.L)    # set all alleles to zero (clonal founder)
popI.set_wildtype(N)


def _printGen():
    af = pop.get_allele_frequencies()
    afI = popI.get_allele_frequencies()
#        afInt = (af*1000).astype(int)
    posAf = af[af > 0]
    posAfI = afI[afI > 0]
    
    posAc = sum(af > 0)
    posAcI = sum(afI > 0)
    detAc = sum(af > detLimit)
    detAcI = sum(afI > detLimit)
    
    SIratio="0:0"
    detSIratio="0:0"
    if posAcI > 0: SIratio = str(posAc/ posAcI)+":1"
    if detAcI > 0: detSIratio = str(detAc/ detAcI)+":1"
    
    #        posAfD = af[af > detlimit]
    #        posAfID = afI[afI > detlimit]
    
    (hist,edges) = np.histogram(posAf,bins=10,range=(0,1))
    (histI,edges) = np.histogram(posAfI,bins=10,range=(0,1))
    noMutations = sum(af > 0) + sum(afI>0)
    noDMutations = sum(af > detLimit) + sum(afI>detLimit)
    
    print "generation:", pop.generation, "days:",pop.generation*2,"\t",
    print "N:",pop.N, 
#    print "clones:",pop.number_of_clones,
    print "mutations:",noMutations, SIratio, 
    print "detectable:",noDMutations,detSIratio,
    print hist,histI


detLimit=0.3

# evolve for 2000 generations and track the allele frequencies
maxgen = 150
allele_frequencies = [pop.get_allele_frequencies()]
tp = [pop.generation]
popStep=1
while pop.generation < maxgen:
#    print "generation:", pop.generation,pop.get_allele_frequencies()
    
    pop.evolve(popStep)
    popI.evolve(popStep)

    # save allele frequencies and time
    # every 200 generations, make one of the deleterious mutations beneficial
    if pop.generation<10:
        allele_frequencies.append(pop.get_allele_frequencies()) 
        tp.append(pop.generation)
        _printGen()

    if (pop.generation % 10 == 0):
        popStep=10
        allele_frequencies.append(pop.get_allele_frequencies()) 
        tp.append(pop.generation)
        _printGen()
        # update fitness function
#        selection_coefficients[m*np.random.randint(0,25)] = 0.01
#        pop.set_trait_additive(selection_coefficients)

# convert to an array to enable slicing
allele_frequencies = np.array(allele_frequencies)

print >>sys.stderr, allele_frequencies
