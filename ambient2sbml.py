"""Export a set of nodes from an AMBIENT network to SBML."""

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

def export_ambient_to_SBML(G, infile, outfile, nodes = -1):
    """Export a set of nodes from an AMBIENT network to SBML."""
    if nodes == -1:
        nodes = G.nodes()
    
    #Establish nodes to be exported
    export_nodes = nodes
    for nodeid in nodes:
        node = G.node[nodeid]
        if node['type'] == 'reaction':
            for nb in G.neighbors(nodeid):
                if nb not in export_nodes:
                    export_nodes.append(nb)
    
    #Get ID for all nodes to be exported
    export_ids = []
    for nodeid in export_nodes:
        export_ids.append(G.node[nodeid]['id'])
    
    #Get original model
    reader = SBMLReader()
    document = reader.readSBMLFromFile(infile)
    model = document.getModel()
    
    ##For each node, if it is not in the export list, delete it
    #for metabolite in model.getListOfSpecies():
    #    sbml_id = metabolite.getId()
    #    if sbml_id not in export_ids:
    #        model.removeSpecies(sbml_id)
    #
    #for reaction in model.getListOfReactions():
    #    sbml_id = reaction.getId()
    #    if sbml_id not in export_ids:
    #        model.removeReaction(sbml_id)
    #
    SBMLWriter.writeSBMLToFile(model, outfile)
    
    return export_ids


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
        genes = re.search('GENE_ASSOCIATION\:([^<]+)<',notes)
        if genes is not None:
            for gene in re.finditer('([^\s\&\|\(\)]+)', genes.group(1)):
                if not gene.group(1) == 'and' and not gene.group(1) == 'or':
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
