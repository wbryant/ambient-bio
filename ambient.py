"""
AMBIENT v1.3Y: Active Modules for Bipartite Networks
Copyright 2012 and 2013 William A. Bryant and John W. Pinney

This module undertakes simulated annealing on a metabolic model (arranged as a
bipartite network) to find connected components of the network with a coordinated
characteristic (for instance upregulation in a particular environment).  Edge
toggling is used to sample the space of possible active modules (based on the
work of Ideker et al).

LICENSE
=======
This program is free software, distributed under the terms of the GNU GPL:
you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

INSTALLATION
============
AMBIENT does not require installation, it can be run straight away.  It is
required that the AMBIENT directory is in the path when run from the command
line / that it is in the Python path when imported.

PREREQUISITES
-------------
Due to compatibility issues in NetworkX it is recommended that Python 2.7 be
used, rather than Python 3

Required non-standard Python packages:

NumPy
libSBML
NetworkX

Example use
===========
A set of example data is available with AMBIENT; from the command line AMBIENT can be run thus (starting in the AMBIENT main directory):

C{cd example}

C{python ../ambient.py -m yeast_4.02.xml -s SCE_scores.tsv -e SCE_pos_log_run}

where C{yeast_4.02.xml} is the SBML yeast model, C{SCE_scores.tsv} is the
tab-separated table of reaction scores and C{SCE_pos_log_run} is the name that
is used in naming all of the results files.

Options for all of the parameters of AMBIENT can be set at the command line.  For
a complete list of all options type (from the AMBIENT main directory):

C{python ambient.py -h}

N.B. the default value for C{N} (the number of iterations) is 10,000, which may
take a few minutes to run.  This value is for testing, a more realistic value
for networks of this size might be ~1,000,000, but would depend on the underlying
structure of the network and the distribution of metabolite/reaction scores.

MAIN FUNCTION DESCRIPTION
=========================

C{G, H, scores, cc_out = ambient(expt_name, Q, N = 10000, M = -1, dir = 1, adaptive_interval = 3000, score_change_ratio = 0.2, intervals_cutoff = 4, H_edges = -1, T_init_div = 10, T_chn_factor = 0.8)}

Inputs
------
Required:

C{expt_name} - string - a name for output files to store the annealing results (no extension required)
C{Q} - networkx DiGraph - a bipartite metabolic network as imported by 'import_SBML_to_bipartite'

Optional:

C{N} - int - number of iterations.
C{M} - int - number of module scores to track in simulated annealing.  If M is set to -1, the simulated annealing is done using logarithmic scoring without limiting the number of modules.
C{dir} - int - if dir = -1, the negative of all reaction scores is taken before simulated annealing proceeds.
C{adaptive_interval} - number of steps to take before assessing how the score has changed for adaptive schedule.
C{score_change_ratio} - percentage rate of change per thousand steps below which temperature is dropped.
C{intervals_cutoff} - number of times to go through adaptive_interval number of steps before automatically dropping temperature.
C{H_edges} - list of tuples of 2 ints - a list of edges that if specified is used as the seed set of edges in the simulated annealing process.
C{T_init_div} - factor by which to reduce the initial temperature (tune to prevent local minima and time wasting at initial high temperature).
C{T_chn_factor} - factor determining rate of temperature reduction during annealing.

Outputs
-------
Output variables:

C{G} - original network.
C{H} - networkx Graph consisting of the subnetwork induced by the final selected edge set.
C{scores} - list of floats consisting of the scores of the top 20 modules in the final subnetwork (in order high to low).
C{cc_out} - an ordered list of tuples where each entry contains the nodes in the 20 top scoring modules (in order high to low scoring).

Output files:

'C{expt_name}.log' which contains details of the edge toggling events that affected the top scoring module, which gives an indication of how the algorithm progressed.
'C{expt_name}.dat' which is a shelf file created by the shelve.open() command in the shelve Python module, which contains the three main outputs and the network used for the simulated annealing.
'C{expt_name}.tsv' which is a table of all nodes in all significant modules (q<0.05).

For help on command-line execution type C{python ambient.py -h}.  
"""

import numpy
from numpy import *
from scipy import *
from math import *
import networkx as nx
import random as rand
from time import time
import operator
import sys, re
import shelve, dumbdbm, pickle
from libsbml import SBMLReader
from copy import *

# Simulated annealing algorithm - after Ideker 2002.  See above for full description of use
#@profile
def ambient(expt_name, Q, N = 10000, M = -1, dir = 1,
            adaptive_interval = 3000, score_change_ratio = 0.2, intervals_cutoff = 4,
            H_edges = -1,
            T_init_div = 10, T_chn_factor = 0.8,
            log_score_cutoff = 0.0):
    """Find high scoring modules in a bipartite metabolic network Q."""   
    
    G_orig = Q.copy()
    G = Q.to_undirected()
    q = len(G.edges())/50
    
    print 'Calculating k ...'
    # Calculate k
    s_tot_m = 0
    s_tot_r = 0
    no_r_pos = 0
    no_m = 0
    for node in G.nodes():
        if G.node[node]['type'] == 'metabolite':
            s_tot_m -= G.node[node]['score']
            no_m += 1
        else:
            #print G.node[node]['score']
            if G.node[node]['score'] >= 0:
                s_tot_r += G.node[node]['score']
                no_r_pos += 1
    if no_r_pos == 0:
        print 'No positive scoring reactions, using negative reactions to clculate k ...'
        s_tot_r = 0
        for node in G.nodes():
            if G.node[node]['type'] == 'reaction':
                s_tot_r += abs(G.node[node]['score'])
                no_r_pos += 1
    s_mean_r = s_tot_r/no_r_pos
    s_mean_m = s_tot_m/no_m
    k = s_mean_r/s_mean_m
    
    print 'k calculated, k = %s.' % k
    
    # Normalise all metabolite scores
    for node in G.nodes():
        if G.node[node]['type'] == 'metabolite':
            G.node[node]['score'] = float(k*G.node[node]['score'])
            
    
    # Ensure all node scores are of type 'float'
    for node in G.nodes():
        
        #print G.node[node]['score']
        G.node[node]['score'] = float(G.node[node]['score'])
        #print G.node[node]['score']
    
    # Store all scores in a dictionary for each node
    score_dict = {}
    for node in G.nodes():
        score_dict[node] = G.node[node]['score']
    
    # If M is not given, logarithmic scoring is used
    if M == -1:
        print 'Using logarithmic scoring.'
        chn_test_fcn = compare_scores_total
    else:
        print 'Using limited module number scoring.'
        chn_test_fcn = compare_graph_scores

    # Estimate typical score changes and set beginning temperature
    module_score_dict = {}
    s_H_sample = []
    print 'Determining initial value for T ...'
    for i in range(500):
        H_edges = rand.sample(G.edges(),q)
        H_edges = set(H_edges)
        H = edges_to_graph(G, H_edges)
        s_H, _ = get_graph_scores(H, M, log_score_cutoff, [], score_dict, module_score_dict)
        s_H_sample.append(max(s_H))
    T = max(s_H_sample) - min(s_H_sample)
    T = T/T_init_div
    print 'Initial value of T: %4.4f.' % (T)
    
    # Set initial number of edges to toggle at each step
    n_toggles = len(G.edges())/50
    print 'Initial value of n_toggles: %d.' % (n_toggles)
    
    # Induce new subnetwork with selected edges
    if H_edges == -1:
        H_edges = rand.sample(G.edges(),q)
    H_edges = set(H_edges)
    H = edges_to_graph(G, H_edges)
    H_n = H.copy()
    G_edges = G.edges()
    H_edges = set(H.edges())
    H_n_edges = set(H_n.edges())
    
    # Determine best M s-scores for this network
    s_H, _ = get_graph_scores(H, M, log_score_cutoff, [], score_dict, module_score_dict)
    
    #Initialise variables
    R = len(G.edges())
    intervals_count = 0
    interval_cutoff_reached = 0
    count_attempts = 0
    score_sum = 0.0
    score_win_prev = adaptive_interval*(sumpos(s_H))
    keep_count = 0
    t1 = time()
    count_score_improvements = 0
    count_score_imps_per_round = 0
    print 'Current s_H: %8.2f %8.2f %8.2f' % (sumpos(s_H), s_H[0], s_H[-1])
    
    # Store run details in dictionary expt_details
    expt_details = {}
    expt_details['name'] = expt_name
    expt_details['N'] = N
    expt_details['M'] = M
    expt_details['dir'] = dir
    expt_details['adaptive_interval'] = adaptive_interval
    expt_details['score_change_ratio'] = score_change_ratio
    expt_details['intervals_cutoff'] = intervals_cutoff
    expt_details['init_edges'] = H_edges
    expt_details['init_temp'] = T
    expt_details['init_ntoggles'] = n_toggles
    expt_details['k'] = k
    expt_details['T_init_div'] = T_init_div
    expt_details['T_chn_factor'] = T_chn_factor
    expt_details['log_score_cutoff'] = log_score_cutoff
    
    
    #Run through N steps of this algorithm
    print 'Begin annealing ...'
    for i in range(1, N):
        
        count_attempts += 1
        score_sum += sumpos(s_H)
        
        #Randomly change the network by a set of edges
        toggle_edges = rand.sample(G_edges,n_toggles)
        toggle_edges = set(toggle_edges)
        toggle_subgraph_edges(H_n, H_n_edges, G, toggle_edges)
        
        # Get scores for the new subgraph
        s_H_n, _ = get_graph_scores(H_n, M, log_score_cutoff, [], score_dict, module_score_dict)
        
        #For adaptive annealing, check sum of positive scoring modules and count improvements
        if sumpos(s_H_n) > sumpos(s_H):
            count_score_improvements += 1
            count_score_imps_per_round += 1
        
        #Compare s scores to decide whether to keep the change
        keep_change = chn_test_fcn(s_H_n, s_H, T)
        if keep_change == 1:
            keep_count += 1
            # update scores to new state
            s_H = s_H_n
            toggle_subgraph_edges(H, H_edges, G, toggle_edges)
        else:
            # revert to previous state
            toggle_subgraph_edges(H_n, H_n_edges, G, toggle_edges)
            
        # Keep user informed of progress and approximate time until completion
        if i/adaptive_interval == i/float(adaptive_interval):
            
            t_now = (time()-t1)/60.0
            time_per_round = t_now/(i)
            time_left = time_per_round * (N-i)
            
            try:
                ro_change = 100000*(score_sum-score_win_prev)/(adaptive_interval*score_win_prev)
            except:
                ro_change = 0
            if abs(ro_change) < score_change_ratio or intervals_count == intervals_cutoff:
                T = T*T_chn_factor
                n_toggles = (n_toggles*8)/10
                if n_toggles < 1:
                    n_toggles = 1
                intervals_count = 0
            else:
                intervals_count += 1
            if intervals_count == intervals_cutoff:
                interval_cutoff_reached = 1
            
            score_win_prev = score_sum
            score_sum = 0
            no_modules = 0
            
            # Output status        
            print 'Step %d:' % (i)
            print 'Approx. time left: %4.2f minutes.' % (time_left)
            print 'Current score_sum = %2.4f.' % score_win_prev
            print 'Current percent rate of change (per 1000 steps, avg.): %5.4f.' % (ro_change)
            print 'Number of improving changes this round: %d.' % (count_score_imps_per_round)
            count_score_imps_per_round = 0
            if interval_cutoff_reached == 1:
                print 'Interval cutoff reached, moving to new temperature.'
                interval_cutoff_reached = 0
            print 'Temperature for following: %5.4f.' % (T)
            print 'Number of toggles per step for following: %d.' % (n_toggles)
            print 'Length of module_score_dict: %d.' % (len(module_score_dict))
            print 'Current score of top 20: %5.4f.' % (sumpos(s_H[0:20]))
            print 'Current s_H sum: %4.2f' % (sumpos(s_H))
            print 'Current highest s_H: %4.2f' % (s_H[0])
            print 'Current second highest s_H: %4.2f\n' % (s_H[1])
            
        # If no improvements are made over a period of adaptive_interval*intervals_cutoff then exit algorithm
        if i/(adaptive_interval*intervals_cutoff) == i/float(adaptive_interval*intervals_cutoff):
            if count_score_improvements == 0:
                break
            else:
                count_score_improvements = 0
        
    print 'Annealing complete.'    
    H = H_n.copy()
    scores, cc_out = get_graph_scores(H, M, log_score_cutoff, [], score_dict, module_score_dict)
    t_diff = (time()-t1)/60
    print 'Time taken for %d steps was %4.2f minutes.' % (i, t_diff)
    
    # Output results to backup file using shelve
    d = shelve.open(cc([expt_name, '.dat']))
    d['G'] = Q
    d['H'] = H
    d['scores'] = scores
    d['cc_out'] = cc_out
    d['expt_details'] = expt_details
    d.close()
    print 'Python outputs written to %s' % (cc([expt_name, '.dat']))
    
    return G_orig, H, scores, cc_out

#Find sum of positive values in a list
def sumpos(lst):
    listpos = 0
    for entry in lst:
        if entry > 0:
            listpos += entry
    return listpos
   
# This function scores a bipartite met/rxn graph by summing
# reaction and metabolite scores and returns the sum 
# for all M of the highest scoring (connected) components.
#@profile
def get_graph_scores(G, M, log_score_cutoff = 0.0, G_comps = [], score_dict = {}, module_score_dict = {}):
    """Get top C{M} module scores in network C{G}."""
    
    if len(G_comps) == 0:
        G_comps = nx.connected_components(G)
        G_comps = [G_comp for G_comp in G_comps]
    elif type(G_comps) is set:
        G_comps = list(G_comps)
    else:
        if type(G_comps[0]) is int:
            G_comps = [G_comps]
    G_comps_frozen = []
    #no_comps = len(G_comps)
    s_tot = []
    
#    score_dict = {}
#    for node in G.nodes():
#            score_dict[node] = G.node[node]['score']
    if not score_dict:
        for node in G.nodes():
            score_dict[node] = G.node[node]['score']
    
    # For each connected component get a score
    for component in G_comps:
        if isinstance(component, int):
            component = [component]
        nodeset = frozenset(component)
        G_comps_frozen.append(nodeset)
        if nodeset in module_score_dict:
            # Use known value for this module
            s_tot.append(module_score_dict[nodeset])
        else:
            # Get score for this component
            s_tot_count = 0
            try:
                for node in component:
                    s_tot_count += score_dict[node]
            except:
                s_tot_count += score_dict[component]
            if isnan(s_tot_count) or isinf(s_tot_count):
                s_tot_count = -100
            
            # Multiply module score by log of module size to encourage congregation
            if M == -1:
                log_count = 2
                if log_score_cutoff > 0:
                    for node in component:
                        if score_dict[node] > log_score_cutoff:
                            log_count += 1
                else:
                    log_count += len(component)
                s_tot_count = log(log_count)*s_tot_count
            
            s_tot.append(s_tot_count)    
            
            # Add score and nodeset to module_score_dict
            module_score_dict[nodeset] = s_tot_count
        
    # Run through all scores in score_dict and remove old modules
    for module in module_score_dict.keys():
        if not module in G_comps_frozen:
            del module_score_dict[module]

    
    # Get a list of component scores sorted by score
    ss = sorted(enumerate(s_tot), key=operator.itemgetter(1),reverse = True)
    
    # With this sorted list, reorder the connected components so that the top
    # scorer is first in ordered_comps, and scores are ordered in ordered_s
    ordered_comps = []
    ordered_s = []
    for comp in ss:
        ordered_comps.append(G_comps[comp[0]])
        ordered_s.append(comp[1])
    
    # Record top score, ordered scores and ordered connected component members
    if M != -1:
        if len(ordered_s) > M:
            ordered_s = ordered_s[:M]
            ordered_comps = ordered_comps[:M]
        
    return ordered_s, ordered_comps

    
##Temperature function - last 5% of toggles are for final annealing, with T = 0.
#def T_func(N = 100000, T_start = 1, T_end = 0.01, T_prop = 0.999999):
#    T = zeros(N)
#    T[0] = T_start
#    for i in range(1,N):
#        T[i] = ((T[i-1]-T_end)*(T_prop)) + T_end
#    return T
    
#This function compares two lists of s_scores, and returns '1' if the change
#should be accepted and '0' if it should not.  s1 is the new s_N and s0 is the
#old one.
def compare_graph_scores(s1, s0, T):
    """Compare two sets of scores C{s1} and C{s2} at temperature C{T} and return a decision as to whether to keep the change."""
    k = 0    
    while k < len(s1):
        if k >= len(s0):
            return 1
        else:
            if s1[k] > s0[k]:
                return 1
            elif s1[k] == s0[k]:
                k += 1
            else:
                try:
                    p = exp((s1[k]-s0[k])/T)
                except:
                    p = 0
                r = rand.uniform(0,1)
                #print '\t%d %5.3f %5.3f %5.3f %5.3f %5.3f' % (k, s1[k], s0[k], s1[k] - s0[k], T, p)
                # f_log.write('\t%d %5.3f %5.3f %5.3f %5.3f\n' % (k, p, r, s1[k], s0[k]))
                if r > p:
                    return 0
                else:
                    k += 1
    #What to do if the algorithm doesn't come to a decision?  Accept the change
    return 1

def compare_scores_total(s1, s0, T):
    """Compare the sum of two sets of scores and determine whether a change should be accepted."""
    
    ss1 = sum(s1)
    ss0 = sum(s0)
    
    if ss1 >= ss0:
        return 1
    else:
        p = exp((ss1-ss0)/T)
        r = rand.uniform(0,1)
        if r > p:
            return 0
        else:
            return 1

# SBML import function
# This function imports an SBML model and converts it into a bipartite network
# of metabolites and reactions, noting weights for metabolites.
def import_SBML_to_bipartite(SBML_filename):
    """Import file C{SBML_filename} and convert it into a bipartite metabolic network."""
    
    print '\n\nImporting SBML model ...'
    
    # Import SBML model
    reader = SBMLReader()
    document = reader.readSBMLFromFile(SBML_filename)
    model = document.getModel()
    print 'Model being imported: %s ...' % (model.getId())
    
    # Initialize NetworkX model and populate with metabolite and reaction nodes.
    # At the same time, create an id / node_idx dictionary for use when creating
    # edges.
    G = nx.DiGraph()
    node_idx = 0
    node_id_dictionary = {}
   
    for metabolite in model.getListOfSpecies():
        node_idx += 1
        G.add_node(node_idx)
        G.node[node_idx]['name'] = metabolite.getName()
        G.node[node_idx]['id'] = metabolite.getId()
        G.node[node_idx]['type'] = 'metabolite'
        node_id_dictionary[metabolite.getId()] = node_idx

    for reaction in model.getListOfReactions():
        node_idx += 1
        G.add_node(node_idx)
        G.node[node_idx]['name'] = reaction.getName()
        G.node[node_idx]['id'] = reaction.getId()
        G.node[node_idx]['type'] = 'reaction'
        node_id_dictionary[reaction.getId()] = node_idx
        
        #print node_idx
        #print G.node[node_idx]['name']
        
        notes = reaction.getNotesString()
        
        genelist = []
        genes = re.search('GENE[_ ]ASSOCIATION\:([^<]+)<',notes)
        if genes is not None:
            for gene in re.finditer('([^\s\&\|\(\)]+)', genes.group(1)):
                if not gene.group(1) == 'and' and not gene.group(1) == 'or' and not gene.group(1) == 'none':
                    genelist.append(gene.group(1))   
        G.node[node_idx]['genelist'] = list(set(genelist))
        
        # Cycle through all reactants and products and add edges
        #print 'REACTANTS:'
        for reactant in reaction.getListOfReactants():
            #print reactant.getSpecies()
            reactant_idx = node_id_dictionary[reactant.getSpecies()]
            #print reactant_idx
            G.add_edge(reactant_idx, node_idx)
            #print G.edges(reactant_idx)
        #print 'PRODUCTS:'
        for product in reaction.getListOfProducts():
            #print product.getSpecies()
            G.add_edge(node_idx,node_id_dictionary[product.getSpecies()])
            #print G.edges(node_idx)
        #print '\n'
    # Add degree of each metabolite as 'weight' attribute
    for node in G.nodes():
        if G.node[node]['type'] == 'metabolite':
            G.node[node]['weight'] = float(G.degree(node))
            G.node[node]['score'] = -1*float(G.degree(node))
    print 'Finished model import.'
    
    return G


# This function takes two sets of nodes from a given graph, from positive
# scores and negative scores respectively, and classifies all of the nodes
# according to those sets (0 if not a member of any of them).  This will only
def classify_nodes_ud(G, set_list1, qvals1, set_list2, qvals2, q_cutoff = 0.05, attribute = 'Classification'):
    """Add attribute C{attribute} to network C{G} indicating two sets of module members from C{set_list1} and C{set_list2}."""
    G_out = G.copy()
    for node in G_out.nodes():
            if G_out.node[node]['type'] == 'metabolite':
                    G_out.node[node][attribute] = -0.5
            else:
                    G_out.node[node][attribute] = 0.0
    class_idx = 0.0
    for idx, module in enumerate(set_list1):
        class_idx += 1
        if qvals1[idx] <= q_cutoff:
            for node in module:
                G_out.node[node][attribute] = class_idx
    class_idx = 0.0
    for idx, module in enumerate(set_list2):
        class_idx -= 1
        if qvals2[idx] <= q_cutoff:
            for node in module:
                G_out.node[node][attribute] = class_idx
    return G_out

# This function takes a list of lists of nodes from a given graph and classifies all of
# the nodes according to those sets (0 if not a member of any of them) in node
# attribute 'attribute'.  If q-values are given, only those modules with
# q-value < q-cutoff are classified
def classify_nodes_single(G, set_list, attribute, q_vals = -1, q_cutoff = 0.05):
    """Add attribute C{attribute} to network C{G} indicating sets of module members from C{set_list}."""
    
    sig_mods = []
    if q_vals != -1:
        for idx, q_val in enumerate(q_vals):
            if q_val <= q_cutoff:
                sig_mods.append(idx)
    
    G_out = G.copy()
    for node in G_out.nodes():
            if G_out.node[node]['type'] == 'metabolite':
                    G_out.node[node][attribute] = -0.5
            else:
                    G_out.node[node][attribute] = 0.0
    class_idx = 0.0
    for idx, nodes in enumerate(set_list):
        class_idx += 1
        if idx in sig_mods or q_vals == -1:
            if type(nodes) == int:
                G_out.node[nodes][attribute] = class_idx
            else:
                for node in nodes:
                    G_out.node[node][attribute] = class_idx
    return G_out

# Take a network, a subnetwork and an edge from the network.  If the edge exists
# in the subnetwork, remove it, otherwise add it and nodes connected to it that
# were not originally in the subnetwork
#@profile
def toggle_subgraph_edges(H, H_edges, G, edges):
    """Add/remove C{edges} to/from C{H}, a subgraph of C{G}."""
    
    """C{H_edges} is the full list of edges in C{H}, passed to improve performance."""
    for edge in edges:
        if edge in H_edges:
            H_edges.remove(edge)
            H.remove_edge(edge[0], edge[1])
            for node in edge:
                if H.degree(node) == 0:
                    H.remove_node(node)  
        else:
            H_edges.add(edge)
            H.add_edge(edge[0], edge[1])
            for node in edge:
                for attr in G.node[node]:
                    H.node[node][attr] = G.node[node][attr]            
    

# This function takes a list of 'edges' from a bipartite reaction/metabolite
# graph 'G' and produces the induced subgraph
def edges_to_graph(G, H_edges):
    """Get subgraph of network C{G} induce by edges C{H_edges}."""
    
    # Get nodes from list of edges with relevant attributes
    H = nx.Graph()
    H.add_edges_from(H_edges)
    for edge in H_edges:
        for node in edge:
            for attr in G.node[node]:
                H.node[node][attr] = G.node[node][attr]
    return H


# Calculate p-values for each module in H (being a subgraph of G), or if H is a list of edges, for each module induced by those edges in G
def get_module_pvalues(H, G, P = 10000, M = 1000, log_score_cutoff = 0.0):
    """Find significance values for all modules with a positive score (connected components) in a subgraph of a reaction/metabolite bipartite network."""
    print 'Calculating p-values and q-values for modules ...'
    r_scores = []
    m_scores = []
    # Get complete lists of reaction and metabolite scores
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            r_scores.append(G.node[node]['score'])
        else:
            m_scores.append(G.node[node]['score'])
    
    # For each module, do P samples of reaction/metabolite sets with the same number of each of those in that module
    if type(H) == 'list':
        H = edges_to_graph(G, H)
    H_ccs = nx.connected_components(H)
    H_cc = []
    H_cc.append(H_ccs[0])
    
    pvals_dict = {}
    cc_pvals = []
    no_pos_scores = 0
    for cc in H_ccs:
        #print cc
        # How many reactions/metabolites?  What is module's score?
        no_rxns = 0
        no_mets = 0
        module_score = 0
        for node in cc:
            module_score += G.node[node]['score']    
            if G.node[node]['type'] == 'reaction':
                no_rxns += 1
            else:
                no_mets += 1
        #print 'mets:rxns - %d:%d' % (no_mets, no_rxns)
        #print 'Module score: %f.' % (module_score)
        
        if module_score < 0:
            cc_pval = 1
        else:
            if len(cc) > 100:
                print len(cc)
            no_pos_scores += 1
            print no_pos_scores
            
            # Do sampling
            #print 'Sample scores:\n'
            cc_sample_values = []
            for i in range(P):
                sample_rxn_vals = rand.sample(r_scores, no_rxns)
                sample_met_vals = rand.sample(m_scores, no_mets)
                
                rxn_score = sum(sample_rxn_vals)
                met_score = sum(sample_met_vals)
                tot_score = rxn_score + met_score
                
                #Take into account logarithmic scoring
                no_pos_scoring = 1
                for val in sample_rxn_vals:
                    if val > log_score_cutoff:
                        no_pos_scoring += 1
                for val in sample_met_vals:
                    if val > log_score_cutoff:
                        no_pos_scoring += 1
                
                tot_score = log(no_pos_scoring)*tot_score
                
                #print tot_score
                cc_sample_values.append(tot_score)
            
            #print 'Sample scores:\n'
            #print cc_sample_values
            
            # How many scores are higher than the score of the module?  What is the p-value?
            score_higher_count = 0
            for score in cc_sample_values:
                if score >= module_score:
                    score_higher_count += 1
            cc_pval = score_higher_count / float(P)
        
            cc_pvals.append(cc_pval)
            cc_name = frozenset(cc)
            if len(cc) > 100:
                print cc_name
            pvals_dict[cc_name] = cc_pval
            if len(cc) > 100:
                print pvals_dict[cc_name]
            #print('cc_name: ' + str(cc_name) + '. pval: ' + str(cc_pval) + '. score: ' + str(module_score) + '\n')
            
    #print 'new modules:\n'
    scores, cc_out = get_graph_scores(H, M, {}, {}, log_score_cutoff)
    pvals_ordered = []
    for id, cc_sorted in enumerate(cc_out):
        print id
        print cc_sorted
        print scores[id]
        if scores[id] < 0:
            break
        cc_sorted_set = frozenset(cc_sorted)
        #print('cc_name: ' + str(cc_sorted) + '. pval: ' + str(pvals_dict[cc_sorted_set]) + '\n')
        pvals_ordered.append(pvals_dict[cc_sorted_set])
    
    # Do correction for multiple testing using Benjamini-Hochberg
    pvals_table = sorted(enumerate(pvals_ordered), key=operator.itemgetter(1))
    pval_index, pvals_sorted = zip(*pvals_table)
    qvals_sorted = calc_bh_values(pvals_sorted)
    qval_table = zip(pval_index, qvals_sorted)
    qvals_ordered_table = sorted(qval_table, key=operator.itemgetter(0))
    _, qvals_ordered = zip(*qvals_ordered_table)
    
    print 'P-values and q-values calculated.'
    return cc_out, scores, pvals_ordered, qvals_ordered

# Calculate p-values for each module in H (being a subgraph of G), or if H is a list of edges, for each module induced by those edges in G
def get_module_qvalues(H, G, score_fcn, P = 10000, M = -1, log_score_cutoff = 0.0):
    """Find significance values for all modules with a positive score (connected components) in a subgraph of a reaction/metabolite bipartite network."""
    print 'Calculating p-values and q-values for modules ...'
    rxn_list = []
    met_list = []
    # Get complete lists of reaction and metabolite scores
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            rxn_list.append(node)
        else:
            met_list.append(node)
    
    # For each module, do P samples of reaction/metabolite sets with the same number of each of those in that module
    if type(H) is list:
        if type(H[0]) is list:
            H_nodes = []
            for cc in H:
                for node in cc:
                    H_nodes.append(node)
        else:
            H_nodes = H
        H = G.subgraph(H_nodes)
    H_ccs = nx.connected_components(H)
    H_cc = []
    #? H_cc.append(H_ccs[0])
    
    pvals_dict = {}
    cc_pvals = []
    no_pos_scores = 0
    cc_no = 0
    for cc in H_ccs:
        #cc_no += 1
        #print cc_no
        #print cc
        # How many reactions/metabolites?
        no_rxns = 0
        no_mets = 0
        for node in cc:   
            if G.node[node]['type'] == 'reaction':
                no_rxns += 1
            else:
                no_mets += 1
        
        #print no_rxns
        #print no_mets
        
        module_score, _ = score_fcn(G, M, log_score_cutoff, cc)
        module_score = module_score[0]
        #print module_score
        
        if module_score <= 0:
            cc_pval = 1
        else:
            no_pos_scores += 1
            print no_pos_scores
            
            # Do sampling
            cc_sample_values = []
            for i in range(P):
                
                sample_rxn_list = rand.sample(rxn_list, no_rxns)
                sample_met_list = rand.sample(met_list, no_mets)
                sample_module_list = []
                for rxn in sample_rxn_list:
                    sample_module_list.append(rxn)
                for met in sample_met_list:
                    sample_module_list.append(met)
                
                sample_score, _ = score_fcn(G, M, log_score_cutoff, sample_module_list)
                sample_score = sample_score[0]
                #print tot_score
                cc_sample_values.append(sample_score)
            
            #print 'Sample scores:\n'
            #print cc_sample_values
            
            # How many scores are higher than the score of the module?  What is the p-value?
            score_higher_count = 0
            for score in cc_sample_values:
                if score >= module_score:
                    score_higher_count += 1
            cc_pval = score_higher_count / float(P)
        
            cc_pvals.append(cc_pval)
            cc_name = frozenset(cc)
            pvals_dict[cc_name] = cc_pval
            #print('cc_name: ' + str(cc_name) + '. pval: ' + str(cc_pval) + '. score: ' + str(module_score) + '\n')
            
    scores, cc_out = score_fcn(H, M, log_score_cutoff)
    #print cc_out
    #print scores
    pvals_ordered = []
    for id, cc_sorted in enumerate(cc_out):
        if scores[id] <= 0:
            if len(pvals_ordered) == 0:
                pvals_ordered.append(1)
            break
        cc_sorted_set = frozenset(cc_sorted)
        pvals_ordered.append(pvals_dict[cc_sorted_set])
        #print pvals_ordered
        
    # Do correction for multiple testing using Benjamini-Hochberg
    pvals_table = sorted(enumerate(pvals_ordered), key=operator.itemgetter(1))
    pval_index, pvals_sorted = zip(*pvals_table)
    qvals_sorted = calc_bh_values(pvals_sorted)
    qval_table = zip(pval_index, qvals_sorted)
    qvals_ordered_table = sorted(qval_table, key=operator.itemgetter(0))
    _, qvals_ordered = zip(*qvals_ordered_table)
    
    print 'P-values and q-values calculated.'
    return cc_out, scores, pvals_ordered, qvals_ordered


            
def get_subgraph_reactions(H, G):
    """Get a set of all reactions connected to an edge in subgraph C{H}, and all reactions present in C{G} and not in C{H}."""
    
    H_reactions = []
    for edge in H.edges():
        for node in edge:
            if H.node[node]['type'] == 'reaction':
                H_reactions.append(node)
    
    G_reactions = []
    for edge in G.edges():
        for node in edge:
            if G.node[node]['type'] == 'reaction':
                G_reactions.append(node)
    
    H_reactions = set(H_reactions)
    G_reactions = set(G_reactions)
    G_only_reactions = (G_reactions - H_reactions)
    
    return H_reactions, G_only_reactions

def get_module_reactions(cc_H, G):
    """Get a set of all reactions in modules C{s_H}, and all reactions present in C{G} and not in C{s_H}."""
    
    H_reactions = []
    for module in cc_H:
        for node in module:
            if G.node[node]['type'] == 'reaction':
                H_reactions.append(node)
    
    G_reactions = []
    for edge in G.edges():
        for node in edge:
            if G.node[node]['type'] == 'reaction':
                G_reactions.append(node)
    
    H_reactions = set(H_reactions)
    G_reactions = set(G_reactions)
    G_only_reactions = (G_reactions - H_reactions)
    
    return H_reactions, G_only_reactions

def output_results_table(expt_name, G, ccomps, qvals, q_cutoff = 0.05):
    """Export all significant modules to a flat file tsv table."""
    f = open(cc([expt_name,'.tsv']),'w')
    for idx, ccomp in enumerate(ccomps):
        if idx >= len(qvals):
            qval = 1
        else:
            qval = qvals[idx]
        if qval <= q_cutoff:
            f.write('#Module ' + str(idx+1) + ' - q-value = ' + str(qvals[idx]) + ', no. of nodes: ' + str(len(ccomp)) + '.\n')
            for node in ccomp:
                n = G.node[node]
                if n['type'] == 'reaction':
                    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\t' + n['type'] + '\n')
            for node in ccomp:
                n = G.node[node]
                if n['type'] == 'metabolite':
                    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\t' + n['type'] + '\n')
            f.write('\n')
    f.close()

#def output_rxn_results_table(expt_name, G, ccomps, qvals, q_cutoff = 0.05):
#    """Export all significant modules to a flat file tsv table."""
#    f = open(cc([expt_name,'_rxns.tsv']),'w')
#    for idx, ccomp in enumerate(ccomps):
#        if qvals[idx] <= q_cutoff:
#            f.write('\tModule ' + str(idx+1) + ' - q-value = ' + str(qvals[idx]) + ', no. of nodes: ' + str(len(ccomp)) + '.\n')
#            for node in ccomp:
#                n = G.node[node]
#                if n['type'] == 'reaction':
#                    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\n')
#    f.close()
#
#    
#def output_met_results_table(expt_name, G, ccomps, qvals, q_cutoff = 0.05):
#    """Export all significant modules to a flat file tsv table."""
#    f = open(cc([expt_name,'_mets.tsv']),'w')
#    for idx, ccomp in enumerate(ccomps):
#        if qvals[idx] <= q_cutoff:
#            f.write('\tModule ' + str(idx+1) + ' - q-value = ' + str(qvals[idx]) + ', no. of nodes: ' + str(len(ccomp)) + '.\n')
#            for node in ccomp:
#                n = G.node[node]
#                if n['type'] == 'metabolite':
#                    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\n')
#    f.close()

#Export a given graph K to a GraphML file with name 'filename', with the
#following provisos:
#
#1)  write_graphml does not support attributes as lists, so these are flattened for export,
#2)  write_graphml does not export int type attributes, so ensure that all numbered attributes are floats.
def gml_export(K, filename, Cyto = -1):
    """Export bipartite graph C{K} in GraphML format to file C{filename}."""
    print 'Exporting graph to GraphML ...'
    G = K.copy()
    for node in G.nodes():
        for attr in G.node[node]:
            if type(G.node[node][attr]) is list:
                if len(G.node[node][attr]) > 1:
                    attr_list = G.node[node][attr][:]
                    G.node[node][attr] = ''
                    for element in attr_list:
                        G.node[node][attr] += str(element)
                        G.node[node][attr] += ', '
                elif len(G.node[node][attr]) == 1:
                    G.node[node][attr] = G.node[node][attr][0]
                else:
                    G.node[node][attr] = 'none'
            if type(G.node[node][attr]) is numpy.float64:
                    G.node[node][attr] = numpy.asscalar(G.node[node][attr])

    nx.write_graphml(G, cc([filename,'.tmp']))
    
    # Write GraphML, when imported into Cytoscape, doesn't use attribute names,
    # but just attribute IDs for identification.  Therefore to make it easier to
    # use in Cytoscape, the IDs in the GraphML ile are converted to the names
    # of the attributes.
    
    # Create dictionary of id/name pairs from the file
    id_name = {}
    f_in = open(cc([filename,'.tmp']), 'r')
    while True:
        current_line = f_in.readline()
        if len(current_line) == 0:
            break
        id_name_line = re.search('(^  <key.+)',current_line)
        if id_name_line is not None:
            id = re.search('id="([^"]+)"', current_line)
            id = id.group(1)
            name = re.search('attr.name="([^"]+)"', current_line)
            name = name.group(1)
            id_name[id] = name
    f_in.close()
    
    # For each line in the temporary GraphML check for ID mentions and replace them with names
    f_in = open(cc([filename,'.tmp']), 'r')
    f_out = open(filename, 'w')
    while True:
        current_line = f_in.readline()
        if len(current_line) == 0:
            break
        id_name_line = re.search('(^  <key.+)',current_line)
        if id_name_line is not None:
            rep_line = re.search('(^.+id=")([^"]+)(".+$)',current_line)
            new_name = id_name[rep_line.group(2)]
            f_out.write(rep_line.group(1) + new_name + rep_line.group(3) + '\n')
        else:
            rep_line = re.search('(^.+key=")([^"]+)(".+$)', current_line)
            if rep_line is not None:
                new_name = id_name[rep_line.group(2)]
                f_out.write(rep_line.group(1) + new_name + rep_line.group(3) + '\n')
            else:
                f_out.write(current_line)
    f_in.close()
    f_out.close()
    
    print 'Graph exported.'
    return G


# This function concatenates the strings in list strings into a single string
# and returns it
def cc(strings):
    """Concatenate strings."""
    str_out = ''
    for string in strings:
        str_out += string
    return str_out
 

# Read scores for reactions and add values to relevant nodes.  Assign the
# median value of changes to each reaction for which there is no score.
def read_rxn_scores(G, filename, unknown_to_zero = -1):
    """Import scores for reactions in C{G} from file C{filename}."""
    print 'Reading reaction scores ...'
    table = get_data_tsv(filename)
    
    if unknown_to_zero == -1:
        # All unassigned reactions get assigned the median score of the scoreset
        print "Unassigned reaction scores set to median."
        ids, scores = zip(*table)
        scores_med = []
        for score in scores:
            scores_med.append(float(score))
        score_default = median(scores_med)
        score_default = float(score_default)
    else:
        print "Unassigned reaction scores set to zero."
        score_default = 0
    
    # Assign scores from table
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            G.node[node]['score'] = float(score_default)
            G.node[node]['no_data'] = True
            for row in table:
                #print str(G.node[node]['id'])
                #print str(row[0])
                if G.node[node]['id'] == row[0]:
                    G.node[node]['score'] = float(row[1])
                    G.node[node]['no_data'] = False
                    break
    print 'Reaction scores read.'
    

    return G


# This function takes gene transcription values and maps them (using simple
# averages) onto reactions in network G (according to gene associations in
# attribute - default is 'genelist')
def genescore2rxnscore(G, gene_scores, gene_attr = 'genelist', unknown_to_zero = -1):
    """Use gene scores to create scores for reactions in C{G}."""
    all_scores = []
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            if G.node[node][gene_attr] is not empty:
                new_genelist = []
                for gene in G.node[node][gene_attr]:
                    gene = gene.lower()
                    gene = gene.rstrip('c')
                    new_genelist.append(gene)
                G.node[node][gene_attr] = new_genelist
                scores = []
                for gene_score in gene_scores:
                    gene_name = gene_score[0]
                    gene_name = gene_name.lower()
                    if gene_name in G.node[node][gene_attr]:
                        scores.append(gene_score[1])
                if len(scores) > 0:
                    mean_score = 0
                    for score in scores:
                        mean_score += float(score)
                    mean_score = mean_score/len(scores)
                    G.node[node]['no_data'] = False
                    all_scores.append(mean_score)
                else:
                    G.node[node]['no_data'] = True
                    mean_score = 0
            else:
                G.node[node]['no_data'] = True
                mean_score = 0
            G.node[node]['score'] = float(mean_score)
    
    #Assign scores to all unknown reactions
    if unknown_to_zero == -1:
        # Assign median score to all reactions without scores    
        score_default = float(median(all_scores))
        print "Unassigned reaction scores set to median."
    else:
        print "Unassigned reaction scores set to zero."
        score_default = 0
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            if G.node[node]['no_data'] == True:
                G.node[node]['score'] = score_default
            G.node[node]['score'] = float(G.node[node]['score'])
            
    return G

def move_median_to_zero(G):
    """Shift all reaction scores by the value of the median, to put the median at zero."""
    print "Median reaction score set to zero."
    # Re-zero reaction scores
    rxn_score_list = []
    for i in G.nodes():
        node = G.node[i]
        if node['type'] == 'reaction':
            rxn_score_list.append(node['score'])
    median_score = median(rxn_score_list)
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            G.node[node]['score'] -= median_score
    return G

def move_median_to_offset(G, offset):
    """Shift all reaction scores by the value of the median, to put the median at zero, then offset median by offset value."""
    
    # Re-zero reaction scores
    rxn_score_list = []
    for i in G.nodes():
        node = G.node[i]
        if node['type'] == 'reaction':
            rxn_score_list.append(node['score'])
    median_score = median(rxn_score_list)
    
    # Move median to offset.
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            G.node[node]['score'] -= (median_score-offset)
    
    print "Median reaction score set to %d." % offset
    return G


#Get any file that is a tab-separated table (without headers) and import it as
#a list of tuples, one for each line of the table.
def get_data_tsv(filename):
    """Import a tab separated file C{filename} into a list of tuples."""
    infile = open(filename,'r')
    table = []
    for line in infile:
        if not list:
            break
        line = line.rstrip()
        data = line.split('\t')
        table.append(data)
    infile.close()
    return table

def export_attributes(G, attribute_list, file_out):
    """Export a tsv with node IDs and selected attributes."""
    outfile = open(file_out,'w')
    for node in G.nodes():
        row_vals = []
        row_vals.append(str(node))
        for attribute in attribute_list:
            if attribute in G.node[node]:
                attr_vals = G.node[node][attribute]
                if type(attr_vals) is str:
                    new_val = attr_vals
                else:
                    if len(attr_vals) > 1:
                        new_val = ''
                        for val_idx, val in enumerate(attr_vals):
                            new_val += val
                            if val_idx+1 < len(attr_vals):
                                new_val += ', '
                    else:
                        new_val = str(attr_vals).strip("'[]")
                row_vals.append(new_val)
            else:
                row_vals.append(G.node[node]['name'])
        for idx, entry in enumerate(row_vals):
            if idx + 1 == len(row_vals):
                outfile.write('{}\n'.format(entry))
            else:
                outfile.write('{}\t'.format(entry))
    outfile.close()


"""
import_channels_genepix
-----------------------

Code for importing and averaging TBDB excel spreadsheets in the GenePix Pro format.  N.B. 'log_score_table' can be used to score reactions in a network by using the 'genescore2rxnscore' function below.

Format notes:

1)     Number of title lines = 21
2)     Column of Channel 1 net intensity (median) = 34
3)     Column of Channel 2 net intensity (median) = 42
4)     Column of identifiers = 8

Function notes:

'files_in' should be a list of names of files that contain experimental replicates
In principle, if the format of any transcription file is known, the optional parameters can be modified to input that file / set of files.

"""

# This function imports a number of GenePix Pro Excel spreadsheets and determines the estimated log fold-change of each transcript, based on their identifier
def import_channels_genepix(files_in, file_out, ch1 = 34, ch2 = 42, id_col = 8, title_lines = 21):
    """Import genepix microarray data from files C{files_in} and return calculated log fold-changes of genes."""
    # Open all files in files_in and lump all data together to get mean changes
    table = []
    for file_in in files_in:        
        # Open file and import complete table, ignoring title lines
        infile = open(file_in,'r')
        no_lines = 1
        for line in infile:
            if no_lines > title_lines:
                if not list:
                    break
                line = line.rstrip()
                data = line.split('\t')
                table.append(data)
            no_lines += 1
        infile.close()
    
    # For each line in the table create a tuple of ID, ch1 net and ch2 net
    #
    # Also, filter out negative values and samples without IDs
    score_table = []
    for line in table:
        if line[id_col-1] is not '' and float(line[ch1-1]) > 0 and float(line[ch2-1]) > 0:
            score_table.append([line[id_col-1], float(line[ch1-1]), float(line[ch2-1])])
    ids, ch1_vals, ch2_vals = zip(*score_table)
    
    # Get all unique identifiers
    unique_ids = []
    for id in ids:
        if id not in unique_ids:
            unique_ids.append(id)
    
    # For each identifier get the mean ch1 and ch2, then get the fold-change
    mean_scores_ch1 = []
    mean_scores_ch2 = []
    mean_scores_rel = []
    for id in unique_ids:
        scores_ch1 = []
        scores_ch2 = []
        for line in score_table:
            if line[0] == id:
                scores_ch1.append(float(line[1]))
                scores_ch2.append(float(line[2]))
        
        mean_ch1 = sum(scores_ch1)/float(len(scores_ch1))
        mean_ch2 = sum(scores_ch2)/float(len(scores_ch2))        
        mean_scores_ch1.append(mean_ch1)
        mean_scores_ch2.append(mean_ch2)
        mean_scores_rel.append(mean_ch2/mean_ch1)

    # Get log_2 fold-change for each identifier
    mean_scores_log2 = []
    for score in mean_scores_rel:
        mean_scores_log2.append(log(score, 2))   
    
    unique_ids_new = []
    for id in unique_ids:
        id = id.lower()
        id = id.rstrip('c')
        unique_ids_new.append(id)
    
    log_score_table = zip(unique_ids_new, mean_scores_log2)
    
    outfile = open(file_out,'w')
    for gene in log_score_table:
        outfile.write(gene[0])
        outfile.write('\t')
        outfile.write(str(gene[1]))
        outfile.write('\n')
    outfile.close()
    
    return log_score_table

def calc_bh_values(p_values, num_total_tests = -1):
    """
    Calculates the Benjamini-Hochberg correction for multiple hypothesis
    testing from a list of p-values *sorted in ascending order*.

    See
    http://en.wikipedia.org/wiki/False_discovery_rate#Independent_tests
    for more detail on the theory behind the correction.

    Parameters:
        - `p_values`: a list or iterable of p-values sorted in ascending order
        - `num_total_tests`: the total number of tests (p-values)

    """
    if num_total_tests == -1:
        num_total_tests = len(p_values)
    prev_bh_value = 0
    bh_values = []
    for i, p_value in enumerate(p_values):
        bh_value = p_value * num_total_tests / float(i + 1)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(bh_value, 1)

        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        bh_values.append(bh_value)
    return bh_values

def inv_transform(G_old, mbc_table):
    """Transform all known metabolite scores by S_met/(1+|mbc_met|)."""
    G = G_old.copy()
    for node in G.nodes():
        if G.node[node]['type'] == 'metabolite':
            for met in mbc_table:
                G.node[node]['score_bak'] = deepcopy(G.node[node]['score'])
                if G.node[node]['name'] == met[0]:
                    G.node[node]['score'] = G.node[node]['score']/(1+abs(float(met[1])))
                    break
    return G       

# Script function - inputs required are:
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Execute AMBIENT.')
    parser.add_argument('-m', action='store', dest='model_file', required=True, help='input metabolic network in SBML format')
    parser.add_argument('-s', action='store', dest='score_file', help='input scores from a tsv, using reaction IDs')
    parser.add_argument('-g', action='store', dest='gene_score_file', help='input gene scores from a tsv to score reactions')
    parser.add_argument('-e', action='store', dest='expt_name', required=True, help='output name for results files')
    parser.add_argument('-N', action='store', dest='N', type=int, default = 10000, help='number of edge toggles')
    parser.add_argument('-M', action='store', dest='M', type=int, default = -1, help='number modules to track')
    parser.add_argument('-d', action='store', dest='d', type=int, default = 1, help='direction of search')
    parser.add_argument('-P', action='store', dest='P', type=int, default = 10000, help='number of tests for empirical significance testing')
    
    #adaptive_interval
    parser.add_argument('-i', action='store', dest='adaptive_interval', type=int, default = 3000, help='number of steps before testing score change for adaptive annealing')
    
    #score_change_ratio
    parser.add_argument('-r', action='store', dest='score_change_ratio', type=int, default = 0.2, help='percentage cutoff for adaptive temperature change')
    
    #intervals_cutoff
    parser.add_argument('-c', action='store', dest='intervals_cutoff', type=int, default = 6, help='number of step intervals before automatic temperature reduction')
    
    #Metabolomics values
    parser.add_argument('-t', action='store', dest='metabolomics_file', help='input metabolomics data from a tsv to score metabolites')
    
    #T_init_div
    parser.add_argument('-T', action='store', dest='T_div', type=float, default=10, help='Temperature divisor to fix initial looseness of negative score acceptance')
    
    #T_chn_factor
    parser.add_argument('-U', action='store', dest='T_mult', type=float, default=0.8, help='Temperature multiplier to determine rate of cooling')
    
    #Set median scores to zero
    parser.add_argument('-Y', action='store', dest='median_to_offset', default=None, help='Set flag to move reaction median score to offset value.')
    
    #Set unknowns to zero
    parser.add_argument('-Z', action='store_true', dest='unknown_to_zero', help='Set flag to zero reactions without assigned scores.')
    
    #Set score cutoff ratio for encouraging coagulation of modules
    parser.add_argument('-V', action='store', dest='log_score_cutoff', type=float, default=0.0, help='Set to a value between 0 and 1 to reduce inflatory effects of low scoring nodes to modules.')
    
    ##Use network in this shelve file if it is specified, instead of the model file provided
    #parser.add_argument('-s', action='store', dest='score_file', help='input scores from a tsv, using reaction IDs')
    #
    args = parser.parse_args()
    
    print("AMBIENT v1.3 running ...")
    
    # Import model and data
    if args.model_file[-3:] == "dat":
        ## File is a storage file and G should be extracted straight from it - network must be called 'G' if Shelve has been used.
        
        ## Try Pickle first
        try:
            f_in = open(args.model_file, "r")
            G = pickle.load(f_in)
            print("Pickle file '%s' being used for model ..." % args.model_file)
        except:
            try:
                d = shelve.Shelf(dumbdbm.open(args.model_file, 'c'))
                G = d['G']
                print("Shelve file '%s' being used for model ..." % args.model_file)
            except:
                print("Unable to load metabolic network, please check input file '%s' ..." % args.model_file)
    else:
        ## File is an SBML model
        G = import_SBML_to_bipartite(args.model_file)
        print("SBML file '%s' imported ..." % args.model_file)
    
    #Determine whether to set unknowns to zero
    if args.unknown_to_zero and args.median_to_offset is None:
        unknown_to_zero = 1
    else:
        unknown_to_zero = -1        
    
    # Depending on whether a reaction score file or a gene score file is selected, import the relevant scores and integrate them with the model
    if not args.score_file is None:
        print 'Reaction scores given.'
        G = read_rxn_scores(G, args.score_file, unknown_to_zero)
    elif not args.gene_score_file is None:
        print 'Gene scores given.'
        gene_score_table = get_data_tsv(args.gene_score_file)
        G = genescore2rxnscore(G, gene_score_table, unknown_to_zero = unknown_to_zero)
    else:
        print 'No input score file was specified.  Please specify a score file with the -s or -g argument.'
        sys.exit()
    expt = args.expt_name
    
    if args.d == -1:
        # Determine direction for reaction scores
        for node in G.nodes():
            if G.node[node]['type'] == 'reaction':
                G.node[node]['score'] = -1*G.node[node]['score']
    
    
    #If flag is set, move reaction median score to zero.
    if not args.median_to_offset is None:

        G = move_median_to_offset(G, float(args.median_to_offset))

                
    if not args.metabolomics_file is None:
        print 'Metabolomics data given, integrating with weight scores.'
        mbc_table = get_data_tsv(args.metabolomics_file)
        G = inv_transform(G, mbc_table)
    else:
        print 'No metabolomics data given.'
    
#    for node in G.nodes():
#        if G.node[node]['type'] == 'metabolite':
#            print G.node[node]['score']
#        
    # Execute simulated annealing
    G, H, _, ccomps = ambient(expt, G, args.N, args.M, args.d, args.adaptive_interval, args.score_change_ratio, args.intervals_cutoff, -1, args.T_div, args.T_mult, args.log_score_cutoff)
    
#    for node in G.nodes():
#        if G.node[node]['type'] == 'metabolite':
#            print G.node[node]['score']
#    
    
    #Get significance for each module
    _, _, _, qvals = get_module_qvalues(H, G, get_graph_scores, args.P, args.M, args.log_score_cutoff)
    
    #Output flat file of significant results
    output_results_table(expt, G, ccomps, qvals)
    #output_met_results_table(expt, G, ccomps, qvals)
    
    # Classify nodes in the modules found and export a GraphML file including
    # those data.  Also store intermediate networks in the .dat file.
    G_class = classify_nodes_single(G, ccomps, 'modules')
    G_output = gml_export(G_class, cc([expt, '.graphml']))
    print 'Output graph written to %s' % (cc([expt, '.graphml']))
    
    d = shelve.open(cc([expt, '.dat']))
    
    #Update experiment details
    expt_details = d['expt_details']
    expt_details['T_div'] = args.T_div
    expt_details['T_mult'] = args.T_mult
    expt_details['median_to_offset'] = args.median_to_offset
    if args.unknown_to_zero:
        expt_details['unknown_to_zero'] = True
    else:
        expt_details['unknown_to_zero'] = False
    expt_details['log_score_cutoff'] = args.log_score_cutoff
    d['expt_details'] = expt_details
    
    #Add results analysis
    d['qvals'] = qvals
    d['G_class'] = G_class
    d['G_output'] = G_output
    
    d.close()
    