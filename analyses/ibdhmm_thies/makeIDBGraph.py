#!/usr/bin/python

import sys
import argparse
from string import *
import os.path as path
import os
import gzip
import cPickle as pickle

from graph_tool.all import *

parser = argparse.ArgumentParser(description='make graph from IBD pc output')

parser.add_argument('--ibd','-f', action="store", dest='file', help='IBD results to parse', nargs='?', default=None)

args = parser.parse_args()

gr = Graph(directed=False)
vs = list()
weight = gr.new_edge_property("double")            # Double-precision floating point
sample = gr.new_vertex_property("string")            # Double-precision floating point

i = -1;
for l in open(args.file,"r"):
    if i==-1: i+=1; continue
    fr,to,f,prob,sites = l.split()
    if fr not in vs:
       frV = gr.add_vertex()
       vs += [fr]
       sample[frV]=fr
    frI = vs.index(fr)
    if to not in vs:
        toV = gr.add_vertex()
        vs += [to]
        sample[toV]=to
    toI = vs.index(to)
    e = gr.add_edge(frI,toI)
    weight[e] = float(f)*10
    i+=1
    
print vs
print weight
#graph_draw(gr, vertex_font_size=8,
#    output_size=(200, 200), output=args.file+".png")

lay = fruchterman_reingold_layout(gr, r=500, weight=weight)
#lay = arf_layout(gr, weight=weight)
 
graph_draw(gr, pos=lay, vertex_text=sample, vertex_font_size=14,
    output_size=(600, 600), edge_pen_width=weight, output=args.file+".png")
