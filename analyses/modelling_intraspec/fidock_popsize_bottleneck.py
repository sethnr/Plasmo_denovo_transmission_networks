#!/bin/python

import sys
import os 
#sys.path.insert(0, '../pkg/python')
import argparse
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import FFPopSim as h
import random

##defaults:
L = int(1e4)  #no_loci
lowC=int(1.5e6)    #capacity
highC=1e10    #capacity
C=int(highC)
N=36       #init pop


parser = argparse.ArgumentParser(description='run drift model of within-host allele frequencies')

parser.add_argument('-L','--loci', action="store", dest='L', type=int, help='number of (bp) loci', nargs='?', default=L)
parser.add_argument('-C','--capacity', action="store", dest='C', type=int, help='max indivs in host', nargs='?', default=C)
parser.add_argument('-g','--genome', action="store", dest='genome', type=int, help='size of full genome (bp)', nargs='?', default=23e6)
parser.add_argument('-G','--growth', action="store", dest='G', type=int, help='growth rate (merozoites per schizont)', nargs='?', default=16)
parser.add_argument('--nomultiply', action="store_false", dest='multip', help='multiple by genome eize')
#parser.add_argument('--bottleneck', action="store_true", dest='bottle', help='bottleneck 10pc before each generation')
parser.add_argument('--bottleneck','-b', action="store", dest='bottle', type=float, help='bottleneck 10pc before each generation', default=-1, nargs='?')
parser.add_argument('--popsizes', action="store", dest='popsizes', help='popsizes derived from fidock model')
parser.add_argument('--array', action="store_true", dest='array', help='part of job array')


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
if args.bottle > 0:
    F = F+".btl"+str(args.bottle)
if args.multip:
    F = F+".mtp"

if args.array:
    taskid = os.getenv('SGE_TASK_ID')
    if taskid is None:
        taskid = os.getenv('LSB_JOBINDEX')        
    F=F+".tsk"+str(taskid)
outfile = open("fidock.C"+str(args.C/1e6)+"m.L"+str(args.L/1e6)+"mb.G"+str(args.G)+F+".txt","w")
afreqs = open("fidock.C"+str(args.C/1e6)+"m.L"+str(args.L/1e6)+"mb.G"+str(args.G)+F+".AFs.txt","w")

SNPmutRate=SNPmutRate*m
INDELmutRate=INDELmutRate*m

# set up population
pop = h.haploid_highd(args.L)                        # produce an instance of haploid_highd with L loci
pop.growth_rate=args.G
pop.carrying_capacity = args.C                       # set the average population size to 50000

pop.outcrossing_rate = 0                        # remove outcrossing
pop.crossover_rate = 0                          # set the crossover rate to zero (core genome)
pop.mutation_rate = SNPmutRate + INDELmutRate   # per locus mutation rate equal to 0.1/N
#pop.mutation_rate = INDELmutRate   # per locus mutation rate equal to 0.1/N
selection_coefficients = 0.0*np.ones(pop.L)     # most loci are neutral
pop.set_trait_additive(selection_coefficients)  # trait 0 is by default fitness
initial_allele_frequencies = np.zeros(pop.L)    # set all alleles to zero (clonal founder)
pop.set_wildtype(N)

popSizes = open(args.popsizes,'r')
maxC=dict()
for L in popSizes:
    (T,C) = L.rstrip().split("\t")
#    print T+" -> "+C
    if T != "T" and C != "NA":
        maxC[int(T)] = long(C)

def _printGen():
    af = pop.get_allele_frequencies()
    posAf = af[af > 0]
    
    posAc = sum(af > 0)
    detAc = sum(af > detLimit)
    
    (hist,edges) = np.histogram(posAf,bins=20,range=(0,maxAF))

    noMutations = sum(af > 0)
    noDMutations = sum(af > detLimit)
    noMaxMutations = sum(af > maxAF)
    if pop.generation==1:
        print >>outfile, "\t".join(map(str,["gen","size"] + [x/float(100) for x in range(0,int(maxAF*100),int(maxAF*5))] + [">"+str(maxAF)]))

    print "generation:", pop.generation, "days:",pop.generation*2,"\t",
    print "N: %.2g" % pop.N, 
    print "C: %.2g" % pop.carrying_capacity, 
    print "%:",round(pop.N/float(pop.carrying_capacity),2), 
    print "part_ratio:",round(pop.participation_ratio,2),
##    print "clones:",pop.number_of_clones,
    print "mutations:",noMutations, # SIratio, 
    print "detectable:",noDMutations, # detSIratio,
    print hist + [noMaxMutations] #,histI

#    alleles = pop.get_allele_frequencies()
#    alleles.sort()
    posAf.sort()
#    print alleles[alleles > 0] #,histI
    print >>outfile, "\t".join(map(str,[pop.generation,pop.N,pop.carrying_capacity] + hist.tolist() + [noMaxMutations]))
    if pop.generation in [14,90]:
        print >>afreqs, "\t".join(map(str,[pop.generation] + posAf.tolist()))

detLimit=0.05

# evolve for 2000 generations and track the allele frequencies
maxgen = 91
popStep=1
#bottlenext (attrition due to lack of RBC)
#nb: no of generation to KILL
BN=args.bottle

maxAF =0.25

#immuneThresh = 1e5 #size when immune system kicks in
immuneThresh = 14 #generation when immune system kicks in
bn=False

while pop.generation < maxgen:
    if pop.generation *2 in maxC and maxC[pop.generation*2] > 1:
        newMax = maxC[long(pop.generation*2)]
#        print "new C = %.2g" % newMax
        if newMax < 1e9:
            pop.carrying_capacity = newMax
            print "setting new max: "+str(newMax)
#    allele_frequencies.append(pop.get_allele_frequencies()) 
    
    if args.bottle > 0 and pop.generation > immuneThresh:
        bn=True
    if bn and random.randint(1, 10)==1:
        print "bottlenecking",pop.N/args.bottle
        pop.bottleneck(pop.N*args.bottle)
    pop.evolve(popStep)
    
    _printGen()


