"""
AMBIENT v0.6.3: Active Modules for Bipartite Networks
Copyright 2012 William A. Bryant and John W. Pinney

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
Two sets of example data are available with AMBIENT; presented here is the set of commands
required to run the analysis of transcriptomic data for Mycobacterium
tuberculosis invading the Human Macrophage - see the README for reference details
(assuming you are in the C{examples/MTU} directory and C{ambient.py} is in your path):

C{########
# Import ambient
import ambient as amb}

C{# Load metabolic network
G = amb.import_SBML_to_bipartite('inj661.xml')}

C{# Read reaction scores into the network
G = amb.read_rxn_scores(G, 'MTU_scores.txt')}

C{# Execute simulated annealing for positive and negative scores.
H_up, scores_up, cc_up = amb.ambient('MTU_example_up', G)
H_down, scores_down, cc_down = amb.ambient('MTU_example_down', G)}

C{# If a network export is required (for instance for visualisation) currently
# GraphML export is implemented}

C{# Map modules onto the metabolic network
G_class = amb.classify_nodes_ud(G, cc_up, cc_down, 'modules')}

C{# Export network to GraphML format.  GraphML cannot handle lists as attribute
# values, so the network attributes are all flattened.  The final output
# network is G_output
G_output = amb.gml_export(G_class, 'MTU_example.graphml')
########}

N.B. the default value for C{N} (the number of iterations) is 10000, which may
take a few minutes to run.  This value is for testing, a more realistic value
for networks of this size would be ~1,000,000, but would depend on the underlying
structure of the network and the distribution of metabolite/reaction scores.

MAIN FUNCTION DESCRIPTION
=========================

C{G, H, scores, cc = ambient(expt_name, Q, N = 10000, M = 20, dir = 1, adaptive_interval = 3000, score_change_ratio = 0.2, intervals_cutoff = 6, H_edges = -1)}

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
from scipy.stats import mode
from scipy.misc.common import factorial
import networkx as nx
import random as rand
from time import time
import operator
import sys
import re
import shelve
import argparse
from libsbml import *

# Simulated annealing algorithm - after Ideker 2002.  See above for full description of use
#@profile
def ambient(expt_name, Q, N = 10000, M = -1, dir = 1, adaptive_interval = 3000, score_change_ratio = 0.2, intervals_cutoff = 6, H_edges = -1):
    """Find high scoring modules in a bipartite metabolic network Q."""

    G = Q.to_undirected()
    q = len(G.edges())/50
    
    # Determine direction for reaction scores
    if dir == -1:
        for node in G.nodes():
            if G.node[node]['type'] == 'reaction':
                G.node[node]['score'] = -1*G.node[node]['score']
    
    # Calculate k
    s_tot_m = 0
    s_tot_r = 0
    no_r_pos = 0
    no_m = 0
    for node in G.nodes():
        if G.node[node]['type'] == 'metabolite':
            s_tot_m += G.node[node]['weight']
            no_m += 1
        else:
            G.node[node]
            if G.node[node]['score'] > 0:
                s_tot_r += G.node[node]['score']
                no_r_pos += 1
    s_mean_r = s_tot_r/no_r_pos
    s_mean_m = s_tot_m/no_m
    k = s_mean_r/s_mean_m
    
    # Make metabolite 'score' = -k*weight and sum all scores
    for node in G.nodes():
	if G.node[node]['type'] == 'metabolite':
	    G.node[node]['score'] = -k*G.node[node]['weight']
    
    # Store all scores in a dictionary for each node
    score_dict = {}
    for node in G.nodes():
	score_dict[node] = G.node[node]['score']
    
    # Estimate typical score changes and set beginning temperature
    module_score_dict = {}
    s_H_sample = []
    for i in range(100):
	H_edges = rand.sample(G.edges(),q)
	H_edges = set(H_edges)
	H = edges_to_graph(G, H_edges)
	s_H, _ = get_graph_scores(H, M, score_dict, module_score_dict)
	s_H_sample.append(max(s_H))
    T = max(s_H_sample) - min(s_H_sample)
    T = T/5.0
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
    s_H, _ = get_graph_scores(H, M, score_dict, module_score_dict)
    
    #Initialise variables
    R = len(G.edges())
    intervals_count = 0
    interval_cutoff_reached = 0
    count_attempts = 0
    score_sum = 0.0
    score_win_prev = adaptive_interval*(sum(s_H))
    keep_count = 0
    t1 = time()
    count_score_improvements = 0
    count_score_imps_per_round = 0
    print 'Current s_H: %8.2f %8.2f %8.2f' % (sum(s_H), s_H[0], s_H[-1])
    
    #Run through N steps of this algorithm
    print 'Begin annealing ...'
    for i in range(1, N):
	
	count_attempts += 1
	score_sum += sum(s_H)
	
        #Randomly change the network by a set of edges
        toggle_edges = rand.sample(G_edges,n_toggles)
	toggle_edges = set(toggle_edges)
        toggle_subgraph_edges(H_n, H_n_edges, G, toggle_edges)
        
        # Get scores for the new subgraph
        s_H_n, _ = get_graph_scores(H_n, M, score_dict, module_score_dict)
	
	#Compare s scores to decide whether to keep the change
        keep_change = compare_graph_scores(s_H_n, s_H, T)
        if sum(s_H_n) > sum(s_H):
	    count_score_improvements += 1
	    count_score_imps_per_round += 1
	
        if keep_change == 1:
	    keep_count += 1
            # update scores
            s_H = s_H_n
            toggle_subgraph_edges(H, H_edges, G, toggle_edges)
        else:
            toggle_subgraph_edges(H_n, H_n_edges, G, toggle_edges)
	    
	
	
        # Keep user informed of progress and approximate time until completion
        if i/adaptive_interval == i/float(adaptive_interval):
	    
            t_now = (time()-t1)/60.0
            time_per_round = t_now/(i)
            time_left = time_per_round * (N-i)
            
	    ro_change = 100000*(score_sum-score_win_prev)/(adaptive_interval*score_win_prev)
	    if abs(ro_change) < score_change_ratio or intervals_count == intervals_cutoff:
		T = T*0.8
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
	    print 'Current score of top 20: %5.4f.' % (sum(s_H[0:20]))
	    print 'Current s_H sum: %4.2f' % (sum(s_H))
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
    scores, cc_out = get_graph_scores(H, M, score_dict, module_score_dict)
    t_diff = (time()-t1)/60
    print 'Time taken for %d steps was %4.2f minutes.' % (i, t_diff)
    
    # Output results to backup file using shelve
    d = shelve.open(cc([expt_name, '.dat']))
    d['G'] = G
    d['H'] = H
    d['scores'] = scores
    d['cc_out'] = cc_out
    d.close()
    print 'Python outputs written to %s' % (cc([expt_name, '.dat']))
    
    return G, H, scores, cc_out



   
# This function scores samples from a bipartite met/rxn graph by summing
# reaction scores and summing metabolite weights and returns the sum 
# Z = z_r - (k * z_m) for all M of the highest scoring (connected) components.
# For improved annealing the top M s-scores must be returned from s_score.
#@profile
def get_graph_scores(G, M, score_dict = {}, module_score_dict = {}, score = 'score'):
    """Get top C{M} module scores in network C{G}."""
    G_comps = nx.connected_components(G)
    G_comps_frozen = []
    no_comps = len(G_comps)
    s_tot = []
    
#    score_dict = {}
#    for node in G.nodes():
#	    score_dict[node] = G.node[node]['score']
    if not score_dict:
	for node in G.nodes():
	    score_dict[node] = G.node[node]['score']
    
    # For each connected component get a score
    for component in G_comps:
	nodeset = frozenset(component)
	G_comps_frozen.append(nodeset)
	if nodeset in module_score_dict:
	    # Use known value for this module
	    s_tot.append(module_score_dict[nodeset])
	else:
	    # Get score for this component
	    s_tot_count = 0
	    for node in component:
		s_tot_count += score_dict[node]
	    if isnan(s_tot_count) or isinf(s_tot_count):
		s_tot_count = -100
	    
	    # Multiply module score by log of module size to encourage congregation
	    if M == -1:
		s_tot_count = log(len(component))*s_tot_count
	    
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
		p = exp((s1[k]-s0[k])/T)
		r = rand.uniform(0,1)
		#print '\t%d %5.3f %5.3f %5.3f %5.3f %5.3f' % (k, s1[k], s0[k], s1[k] - s0[k], T, p)
		# f_log.write('\t%d %5.3f %5.3f %5.3f %5.3f\n' % (k, p, r, s1[k], s0[k]))
		if r > p:
		    return 0
		else:
		    k += 1
    #What to do if the algorithm doesn't come to a decision?  Accept the change
    return 1


# SBML import function
# This function imports an SBML model and converts it into a bipartite network
# of metabolites and reactions, noting weights for metabolites.
def import_SBML_to_bipartite(SBML_filename):
    """Import file C{SBML_filename} and convert it into a bipartite metabolic network."""
    
    # Import SBML model
    reader = SBMLReader()
    document = reader.readSBMLFromFile(SBML_filename)
    model = document.getModel()
    print 'Model being imported: %s.' % (model.getId())
    
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
    print 'Finished model import.'
    return G


# This function takes two sets of nodes from a given graph, from positive
# scores and negative scores respectively, and classifies all of the nodes
# according to those sets (0 if not a member of any of them).
def classify_nodes_ud(G, set_list1, set_list2, attribute):
    """Add attribute C{attribute} to network C{G} indicating two sets of module members from C{set_list1} and C{set_list2}."""
    G_out = G.copy()
    for node in G_out.nodes():
    	if G_out.node[node]['type'] == 'metabolite':
    		G_out.node[node][attribute] = -0.5
    	else:
    		G_out.node[node][attribute] = 0.0
    class_idx = 0.0
    for module in set_list1:
        class_idx += 1
        for node in module:
            G_out.node[node][attribute] = class_idx
    class_idx = 0.0
    for module in set_list2:
        class_idx -= 1
        for node in module:
            G_out.node[node][attribute] = class_idx
    return G_out


# This function takes a set of nodes from a given graph and classifies all of
# the nodes according to those sets (0 if not a member of any of them) in node
# attribute 'attribute'.
def classify_nodes_single(G, set_list, attribute):
    """Add attribute C{attribute} to network C{G} indicating sets of module members from C{set_list}."""
    G_out = G.copy()
    for node in G_out.nodes():
    	if G_out.node[node]['type'] == 'metabolite':
    		G_out.node[node][attribute] = -0.5
    	else:
    		G_out.node[node][attribute] = 0.0
    class_idx = 0.0
    for nodes in set_list:
        class_idx += 1
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
def get_module_pvalues(H, G, M = 1000, P = 10000):
    """Find significance values for all modules (connected components) in a subgraph of a reaction/metabolite bipartite network."""
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
    for cc in H_ccs:
	if len(cc) < 400:
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
	    
	    # Do sampling
	    #print 'Sample scores:\n'
	    cc_sample_values = []
	    for i in range(P):
		rxn_score = sum(rand.sample(r_scores, no_rxns))
		met_score = sum(rand.sample(m_scores, no_mets))
		tot_score = rxn_score + met_score
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
	else:
	    cc_pval = 1
	#print 'p-value: %f\n' % (cc_pval)
        cc_pvals.append(cc_pval)
	cc_name = frozenset(cc)
        pvals_dict[cc_name] = cc_pval
	
	    
    scores, cc_out = get_graph_scores(H, M)
    pvals_ordered = []
    for cc_sorted in cc_out:
	cc_sorted_set = frozenset(cc_sorted)
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
	if qvals[idx] <= q_cutoff:
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
#	if qvals[idx] <= q_cutoff:
#	    f.write('\tModule ' + str(idx+1) + ' - q-value = ' + str(qvals[idx]) + ', no. of nodes: ' + str(len(ccomp)) + '.\n')
#	    for node in ccomp:
#		n = G.node[node]
#		if n['type'] == 'reaction':
#		    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\n')
#    f.close()
#
#    
#def output_met_results_table(expt_name, G, ccomps, qvals, q_cutoff = 0.05):
#    """Export all significant modules to a flat file tsv table."""
#    f = open(cc([expt_name,'_mets.tsv']),'w')
#    for idx, ccomp in enumerate(ccomps):
#	if qvals[idx] <= q_cutoff:
#	    f.write('\tModule ' + str(idx+1) + ' - q-value = ' + str(qvals[idx]) + ', no. of nodes: ' + str(len(ccomp)) + '.\n')
#	    for node in ccomp:
#		n = G.node[node]
#		if n['type'] == 'metabolite':
#		    f.write(n['id'] + '\t' + n['name'] + '\t' + str(n['score']) + '\n')
#    f.close()

#Export a given graph K to a GraphML file with name 'filename', with the
#following provisos:
#
#1)  write_graphml does not support attributes as lists, so these are flattened for export,
#2)  write_graphml does not export int type attributes, so ensure that all numbered attributes are floats.
def gml_export(K, filename):
    """Export bipartite graph C{K} in GraphML format to file C{filename}."""
    print 'Exporting graph to gml ...'
    G = K.copy()
    for node in G.nodes():
        for attr in G.node[node]:
            if type(G.node[node][attr]) is list:
                if len(G.node[node][attr]) > 1:
		    attr_list = G.node[node][attr][:]
		    G.node[node][attr] = ''
                    for element in attr_list:
			G.node[node][attr] += element
			G.node[node][attr] += ', '
                elif len(G.node[node][attr]) == 1:
                    G.node[node][attr] = G.node[node][attr][0]
                else:
                    G.node[node][attr] = 'none'

    nx.write_graphml(G, filename)
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
def read_rxn_scores(G, filename):
    """Import scores for reactions in C{G} from file C{filename}."""
    print 'Reading reaction scores ...'
    table = get_data_tsv(filename)
    
    # All unassigned reactions get assigned the median score of the scoreset
    ids, scores = zip(*table)
    scores_med = []
    for score in scores:
	scores_med.append(float(score))
    score_default = median(scores_med)
    score_default = float(score_default)
    
    # Assign scores from table
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            G.node[node]['score'] = float(score_default)
            for row in table:
		#print str(G.node[node]['id'])
		#print str(row[0])
                if G.node[node]['id'] == row[0]:
                    G.node[node]['score'] = float(row[1])
		    break
    print 'Reaction scores read.'
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

# This function takes gene transcription values and maps them (using simple
# averages) onto reactions in network G (according to gene associations in
# attribute - default is 'genelist')
def genescore2rxnscore(G, gene_scores, gene_attr = 'genelist'):
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
	
    # Assign median score to all reactions without scores    
    median_score = float(median(all_scores))
    for node in G.nodes():
        if G.node[node]['type'] == 'reaction':
            if G.node[node]['no_data'] == True:
		G.node[node]['score'] = median_score
	    G.node[node]['score'] = float(G.node[node]['score'])
    	
    return G


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
    
    
    args = parser.parse_args()
    
    # Import model and data
    G = import_SBML_to_bipartite(args.model_file)
    
    # Depnding on whether a reaction score file or a gene score file is selected, import the relevant scores and integrate them with the model
    if not args.score_file is None:
	print 'Reaction scores given.'
	G = read_rxn_scores(G, args.score_file)
    elif not args.gene_score_file is None:
	print 'Gene scores given.'
	gene_score_table = get_data_tsv(args.gene_score_file)
	G = genescore2rxnscore(G, gene_score_table)
    else:
	print 'No input score file was specified.  Please specify a score file with the -s or -g argument.'
	sys.exit()
    expt = args.expt_name
    
    # If M is not given, logarithmic scoring is used (and M is set to 1000)
    if args.M == -1:
	print 'Using logarithmic scoring.'
    
    # Execute simulated annealing
    G, H, _, ccomps = ambient(expt, G, args.N, args.M, args.d, args.adaptive_interval, args.score_change_ratio, args.intervals_cutoff)
    
    #Get significance for each module
    _, _, _, qvals = get_module_pvalues(H, G, args.M, args.P)
    
    #Output flat file of significant results
    output_results_table(expt, G, ccomps, qvals)
    #output_met_results_table(expt, G, ccomps, qvals)
    
    # Classify nodes in the modules found and export a GraphML file including
    # those data.  Also store intermediate networks in the .dat file.
    G_class = classify_nodes_single(G, ccomps, 'modules')
    G_output = gml_export(G_class, cc([expt, '.graphml']))
    print 'Output graph written to %s' % (cc([expt, '.graphml']))
    
    d = shelve.open(cc([expt, '.dat']))
    d['qvals'] = qvals
    d['G_class'] = G_class
    d['G_output'] = G_output
    d.close()
    