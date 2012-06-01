AMBIENT v0.1

June 2012.  This is a beta release

GETTING STARTED

Language: Python
Non-standard requirements: libSBML and NetworkX

Installation requires only importing the module into Python: 'import ambient'.  Details on use of the functions in the module can be found in the module file, including a short tutorial and information on customising the simulated annealing algorithm.  The program can also be run without installation from the command line (as long as ambient.py is in the path).
 
EXAMPLES

The example can be found in /examples/SCE for Saccharomyces cerevisiae diauxie (baker's yeast).  Details on running the example can be found in the module file.  References for these data and models can be found at the end of this document.

If transcriptional data are in the GenePix Pro Excel format, they can be imported automatically and scores for each reaction in the relevant organism's metabolic network can be inferred.  This is achieved by the 'import_channels_genepix' function, which can take an arbitrary number of results files (representing, say, Biological replicates) and get mean log-fold-changes of transcription for each reaction.  

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

