import numpy
from numpy import *
from scipy import *
from math import *
from scipy.stats import mode
from scipy.misc.common import factorial
import networkx as nx
import random as rand
from time import time
import operator
import os.path
import sys
import re
import shelve
import argparse
from libsbml import *
from ambient import get_data_tsv
import ambient as amb
import copy

def combine_amb_runs(dat1, dat2, outfile, q_cutoff = 0.05, run1_name = 'Run1', run2_name = 'Run2'):
    """Take the results from two runs of AMBIENT and combine their results,
    outputting a GraphML file containing three modules attribute columns.
    
    Module membership is always done using a significance cutoff controlled by q_cutoff.
    
    These columns are: modules from dat1, modules from dat2 and the third column
    shows how they combine.  The values for the third column are as follows:
    
    0 - not in a module in either
    1 - in a module from dat1
    2 - in a module from dat2
    3 - in modules in both."""
    
    dat1 = shelve.open(dat1)
    ccs1 = dat1['cc_out']
    qvals1 = dat1['qvals']
    dat2 = shelve.open(dat2)
    ccs2 = dat2['cc_out']
    qvals2 = dat2['qvals']
    
    G = dat1['G']
    
    G_class = amb.classify_nodes_single(G, ccs1, run1_name, qvals1, q_cutoff)
    G_class = amb.classify_nodes_single(G_class, ccs2, run2_name, qvals2, q_cutoff)
    
    """Determine memberships for third column."""
    members1 = []
    for idx, mod in enumerate(ccs1):
        if idx < len(qvals1):
            if qvals1[idx] <= q_cutoff:
                for node in mod:
                    members1.append(node)
            
    members2 = []
    for idx, mod in enumerate(ccs2):
        if idx < len(qvals2):
            if qvals2[idx] <= q_cutoff:
                for node in mod:
                    members2.append(node)
    
    members_both = []
    for node1 in members1:
        for node2 in members2:
            if node1 == node2:
                members_both.append(node1)
    
    cc_overlap = []
    cc_overlap.append(members1)
    cc_overlap.append(members2)
    cc_overlap.append(members_both)
    
    G_class = amb.classify_nodes_single(G_class, cc_overlap, 'Overlap')
    
    G_output = amb.gml_export(G_class, amb.cc([outfile, '.graphml']))
    
def compare_cc_f_measure(G, ccs1, ccs2, qvals1 = None, qvals2 = None, q_cutoff = 0.05):
    """Find the similarity between the sets of nodes according to two AMBIENT runs."""
    
    if type(ccs1[0]) != 'list':
        ccs1 = [ccs1]
    if type(ccs2[0]) != 'list':
        ccs2 = [ccs2]
        
    if qvals1 is None:
        qvals1 = [0]*len(ccs1)
    if qvals2 is None:
        qvals2 = [0]*len(ccs2)
    
    qvals1 = tuple(qvals1)
    qvals2 = tuple(qvals2)
    
    all_nodes = G.nodes()
    set_of_nodes = set(all_nodes)
    mod_set_1 = []
    for idx, mod in enumerate(ccs1):
        if idx < len(qvals1):
            if qvals1[idx] <= q_cutoff:
                for node in mod:
                    mod_set_1.append(node)
    mod_set_1 = set(mod_set_1)
    mod_set_2 = []
    for idx, mod in enumerate(ccs2):
        if idx < len(qvals2):
            if qvals2[idx] <= q_cutoff:
                for node in mod:
                    mod_set_2.append(node)
    mod_set_2 = set(mod_set_2)
    mod_set_not_1 = set_of_nodes - mod_set_1
    mod_set_not_2 = set_of_nodes - mod_set_2
    
    f_score = f_measure(mod_set_1, mod_set_not_1, mod_set_2, mod_set_not_2)
    return f_score

def f_measure(pos1, neg1, pos2, neg2):
    """Determine the F-measure for a binary classification against another binary classification."""
    
    pos1 = set(pos1)
    neg1 = set(neg1)
    pos2 = set(pos2)
    neg2 = set(neg2)
    
    tp = len(pos1 & pos2)
    fp = len(pos2 - pos1)
    tn = len(neg1 & neg2)
    fn = len(neg2 - neg1)
    
    try:
        precision = tp / float(tp+fp)
        recall = tp / float(tp + fn)
        f_measure = 2*precision*recall/(precision + recall)
    except:
        f_measure = 0
    return f_measure

def convert_cc_to_module_list(G, ccs, qvals, q_cutoff = 0.05):
    """Convert membership lists of modules to ordered list of module memberships."""
    
    ccs_sig = []
    for idx, cc in enumerate(ccs):
        if idx < len(qvals):
            if qvals[idx] <= q_cutoff:
                ccs_sig.append(cc)

    full_node_list = G.nodes()
    non_member_id = len(ccs_sig)+1
    module_no = 1
    
    module_member_list = [non_member_id] * len(full_node_list)
    for mod_idx, mod in enumerate(ccs_sig):
        mod_ref = mod_idx + 1
        for node in mod:
            module_member_list[full_node_list.index(node)] = mod_ref
    return module_member_list

def sig_module_table(expt_name, outfile, q_cutoff = 0.05):
    """For an experiment create a table of module memberships for statistical analysis."""
    
    #import data and get module lists
    run_no = 1
    modules = []
    while os.path.isfile(amb.cc([expt_name, str(run_no), '.dat'])) | os.path.isfile(amb.cc([expt_name, str(run_no), '.dat.dat'])):
        filename = amb.cc([expt_name, str(run_no), '.dat'])
        p = shelve.open(filename)
        cc_out = p['cc_out']
        qvals = p['qvals']
        G = p['G']
        modlist = convert_cc_to_module_list(G, cc_out, qvals)
        modules.append(modlist)
        no_nodes = len(G.nodes())
        run_no += 1
    
    #Output data into tsv
    
    mod_table = zip(*modules)
    f = open(outfile,'w')
    
    for value in range(1,run_no):
        f.write('Run' + str(value) + '\t')
    f.write('\n')
    for row in mod_table:
        for entry in row:
            f.write(str(entry) + '\t')
        f.write('\n')
    f.close()
    
    
    
    