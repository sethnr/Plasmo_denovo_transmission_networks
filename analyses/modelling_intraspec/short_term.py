#!/bin/python

import sys
#sys.path.insert(0, '../pkg/python')
import argparse
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import FFPopSim as h

##defaults:
L = int(1e4)  #no_loci
lowC=int(1.5e6)    #capacity
highC=1e8    #capacity
C=int(highC)
N=1       #init pop


parser = argparse.ArgumentParser(description='run drift model of within-host allele frequencies')

parser.add_argument('-L','--loci', action="store", dest='L', type=int, help='number of (bp) loci', nargs='?', default=L)
parser.add_argument('-C','--capacity', action="store", dest='C', type=int, help='max indivs in host', nargs='?', default=C)
parser.add_argument('-g','--genome', action="store", dest='genome', type=int, help='size of full genome (bp)', nargs='?', default=23e6)
parser.add_argument('-G','--growth', action="store", dest='G', type=int, help='growth rate (merozoites per schizont)', nargs='?', default=16)
parser.add_argument('--nomultiply', action="store_false", dest='multip', help='multiple by genome eize')
parser.add_argument('--bottleneck', action="store_true", dest='bottle', help='bottleneck 10pc before each generation')


args = parser.parse_args()

SNPmutRate=4.07e-10
INDELmutRate=2.4e-9

#locus multiplier
#m = 23e6/(L*2)
if args.multip:
    m = args.genome/args.L
else:
    m=1

print str(args.C),"inds",args.L,"bp","("+str(m)+")"
print "SNPrate:   "+str(SNPmutRate)+" ("+str(SNPmutRate*m)+")"
print "INDELrate: "+str(INDELmutRate)+" ("+str(INDELmutRate*m)+")"

F=""
if args.bottle:
    F = F+".btl"
if args.multip:
    F = F+".mtp"

outfile = open("qstages.AF.C"+str(args.C/1e6)+"m.L"+str(args.L/1e6)+"mb.G"+str(args.G)+F+".txt","w")
afreqs = open("qstages.AF.C"+str(args.C/1e6)+"m.L"+str(args.L/1e6)+"mb.G"+str(args.G)+F+".AFs.txt","w")


SNPmutRate=SNPmutRate*m
INDELmutRate=INDELmutRate*m

# set up population
pop = h.haploid_highd(args.L)                        # produce an instance of haploid_highd with L loci
pop.growth_rate=36
pop.carrying_capacity = args.C                       # set the average population size to 50000

pop.outcrossing_rate = 0                        # remove outcrossing
pop.crossover_rate = 0                          # set the crossover rate to zero (core genome)
pop.mutation_rate = SNPmutRate + INDELmutRate   # per locus mutation rate equal to 0.1/N
#pop.mutation_rate = INDELmutRate   # per locus mutation rate equal to 0.1/N
selection_coefficients = 0.0*np.ones(pop.L)     # most loci are neutral
pop.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness
initial_allele_frequencies = np.zeros(pop.L)    # set all alleles to zero (clonal founder)
pop.set_wildtype(N)


# popI = h.haploid_highd(L)                        # produce an instance of haploid_highd with L loci
# popI.carrying_capacity = C                  # set the average population size to 50000
# popI.outcrossing_rate = 0                        # remove outcrossing
# popI.crossover_rate = 0                          # set the crossover rate to zero (core genome)
# popI.mutation_rate = INDELmutRate                # per locus mutation rate ~10X SNP rate
# selection_coefficients = 0.0*np.ones(popI.L)     # most loci are neutral
# popI.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness
# initial_allele_frequencies = np.zeros(popI.L)    # set all alleles to zero (clonal founder)
# popI.set_wildtype(N)


def _printGen():
    af = pop.get_allele_frequencies()
#    afI = popI.get_allele_frequencies()
#        afInt = (af*1000).astype(int)
    posAf = af[af > 0]
#    posAfI = afI[afI > 0]
    
    posAc = sum(af > 0)
#    posAcI = sum(afI > 0)
    detAc = sum(af > detLimit)
#    detAcI = sum(afI > detLimit)
    
#    SIratio="0:0"
#    detSIratio="0:0"
#    if posAcI > 0: SIratio = str(posAc/ posAcI)+":1"
#    if detAcI > 0: detSIratio = str(detAc/ detAcI)+":1"
    
    #        posAfD = af[af > detlimit]
    #        posAfID = afI[afI > detlimit]
    
    (hist,edges) = np.histogram(posAf,bins=20,range=(0,maxAF))
#    (histI,edges) = np.histogram(posAfI,bins=10,range=(0,1))
#    noMutations = sum(af > 0) + sum(afI>0)
#    noDMutations = sum(af > detLimit) + sum(afI>detLimit)
    noMutations = sum(af > 0)
    noDMutations = sum(af > detLimit)
    noMaxMutations = sum(af > maxAF)
    if pop.generation==1:
        print >>outfile, "\t".join(map(str,["gen","size"] + [x/float(100) for x in range(0,int(maxAF*100),int(maxAF*5))] + [">"+str(maxAF)]))

    print "generation:", pop.generation, "days:",pop.generation*2,"\t",
    print "N:",pop.N, 
    print "%:",round(pop.N/float(pop.carrying_capacity),2), 
    print "part_ratio:",round(pop.participation_ratio,2),
##    print "clones:",pop.number_of_clones,
    print "mutations:",noMutations, # SIratio, 
    print "detectable:",noDMutations, # detSIratio,
    print hist #,histI

    posAf.sort()
    print >>outfile, "\t".join(map(str,[pop.generation,pop.N] + hist.tolist() + [noMaxMutations]))
    print >>afreqs, "\t".join(map(str,[pop.generation] + posAf.tolist()))

detLimit=0.05

# evolve for 2000 generations and track the allele frequencies
maxgen = 150
allele_frequencies = [pop.get_allele_frequencies()]
tp = [pop.generation]
popStep=1
#bottlenext (attrition due to lack of RBC)
#nb: no of generation to KILL
BN=0.8

maxAF =0.25


while pop.generation < maxgen:
#    print "generation:", pop.generation,pop.get_allele_frequencies()
    pop.evolve(popStep)
#    popI.evolve(popStep)

    # save allele frequencies and time
    # every 200 generations, make one of the deleterious mutations beneficial
#    if pop.generation<10:
    if pop.generation in [14,90,180]:
        allele_frequencies.append(pop.get_allele_frequencies()) 
        tp.append(pop.generation)
        _printGen()

#    if (pop.generation % 10 == 0):
#        allele_frequencies.append(pop.get_allele_frequencies()) 
#        tp.append(pop.generation)
#        _printGen()
    pcMax = float(pop.N)/pop.carrying_capacity
    if pcMax > BN and args.bottle:        
#        print "bottlenecking: "+" ".join(map(str,[args.C,pop.N,round(pcMax,2),int(pop.N*(1-(BN*pcMax)))]))
        print "bottlenecking: "+" ".join(map(str,[args.C,pop.N,round(pcMax,2),int(pop.carrying_capacity*BN)])),
        print str(pop.N)+" -> ",
        pop.bottleneck(pop.N)
        print pop.N
# convert to an array to enable slicing
allele_frequencies = np.array(allele_frequencies)

print >>sys.stderr, allele_frequencies
