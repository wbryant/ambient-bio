AMBIENT v0.4

June 2012.  This is a beta release

GETTING STARTED

Language: Python
Non-standard requirements: libSBML and NetworkX
AMBIENT has been tested and works with Python v2.7, libSBML v5.3.0, NetworkX v1.6

Installation requires only importing the module into Python: 'import ambient'.  Details on use of the functions in the module can be found in the module file, including a short tutorial and information on customising the simulated annealing algorithm.  The program can also be run without installation from the command line (as long as ambient.py is in the path).
 
EXAMPLE

The example can be found in /example for Saccharomyces cerevisiae diauxie (baker's yeast).  Details on running the example can be found in the module file.  References for these data and models can be found at the end of this document.

Basic command line usage, from the /example directory type:

'python ../ambient.py -m yeast_4.02.xml -s SCE_scores.txt -e SCE_pos_log_run -N 1000000 -P 10000'

which runs the algorithm for up to 1,000,000 steps and bases the empirical significance values for the found modules on 10,000 random samples per module.  This may take several hours, depending on the processor speed of the computer used.

FULL COMMAND LINE USAGE

usage: ambient.py [-h] -m MODEL_FILE [-s SCORE_FILE] [-g GENE_SCORE_FILE] -e
                  EXPT_NAME [-N N] [-M M] [-d D] [-P P] [-i ADAPTIVE_INTERVAL]
                  [-r SCORE_CHANGE_RATIO] [-c INTERVALS_CUTOFF]

Execute AMBIENT.

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL_FILE         input metabolic network in SBML format
  -s SCORE_FILE         input scores from a tsv, using reaction IDs
  -g GENE_SCORE_FILE    input gene scores from a tsv to score reactions
  -e EXPT_NAME          output name for results files
  -N N                  number of edge toggles
  -M M                  number modules to track
  -d D                  direction of search
  -P P                  number of tests for empirical significance testing
  -i ADAPTIVE_INTERVAL  number of steps before testing score change for
                        adaptive annealing
  -r SCORE_CHANGE_RATIO
                        percentage cutoff for adaptive temperature change
  -c INTERVALS_CUTOFF   number of step intervals before automatic temperature
                        reduction

OUTPUTS

EXPT_NAME.graphml which is a GraphML file for visualisation.

EXPT_NAME.DAT which is a shelf file created by the shelve.open() command in the shelve module, which contains the three main outputs and the network used for the simulated annealing.

EXPT_NAME.TSV which is a table of all nodes in all significant modules (q<0.05).

N.B. The Python shelve module may behave slightly differently depending on which system it is run on.  In certain cases it will produce 3 files, with the extensions '.dat', '.dir' and '.bak'.  The file with '.dat' (or '.dat.dat') is the full shelve file and can be used with shelve.

GENEPIX DATA

If transcriptional data are in the GenePix Pro Excel format, they can be imported automatically and scores for each reaction in the relevant organism's metabolic network can be inferred.  This is achieved using the 'import_channels_genepix' function, which can take an arbitrary number of results files (representing, say, Biological replicates) and get mean log-fold-changes of transcription for each gene and reaction.  

CONTACT

If you have any comments or queries, please contact Dr William Bryant, w.bryant@imperial.ac.uk.

LICENSE

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

CREDITS

Code and concept are by Dr William A. Bryant and Dr John W. Pinney.

REFERENCES

Yeast Model: Herrgård, M. J., Swainston, N., Dobson, P., Dunn, W. B., Arga, K. Y., Arvas, M., Blüthgen, N., et al. (2008). A consensus yeast metabolic network reconstruction obtained from a community approach to systems biology. Nature Biotechnology, 26(10), 1155–1160. doi:10.1038/nbt1492

Yeast data: DeRisi, J. L. (1997). Exploring the Metabolic and Genetic Control of Gene Expression on a Genomic Scale. Science, 278(5338), 680–686. doi:10.1126/science.278.5338.680

