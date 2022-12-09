# Dominating Sets
This repository demonstrates the use of iterative projection methods to solve a classic discrete problem: finding dominating sets. For a given network, a dominating set is a subset of the network such that every node in the complete network is either in the dominating subset or shares an edge with a member of the subset. The algorithm in `ds.c` is one way of finding dominating sets. See arXiv:2106.05206 for details about this algorithm (and similar algorithms for a few other problems).  

To compile: `gcc -O2 ds.c -lm -o ds`  
To run: `./ds netfile D epsilon beta maxiter iterstride stoperr trials id`  
- `netfile` : name of the file containing the network connections (see below for formatting info)  
- `D` : the maximum cardinality of the dominating set to look for  
- `epsilon` : determines how quickly to tune the metric parameter; generally epsilon≤100/maxiter is reliable (assuming maxiter>>100); set epsilon=0 for no tuning  
- `beta` : a hyperparameter of the iteration (see arXiv:2106.05206); values close to 0 are slow but reliable, values close to 1 can be fast but might get stuck in limit cycles of other misbehavior; start with 0.5  
- `maxiter` : the maximum number of iterations for the algorithm to make before giving up  
- `iterstride` : the number of iterations between each print to the .err file (the errors provide a diagnostic of the algorithm's progress)  
- `stoperr` : the error at which the algorithm will stop and declare it has found a solution; try a small nonzero value like 0.001 or 0.000001  
- `name` : the name for the output files (`name.err`, `name.sol`, `name.stats`)  

An example run: `./ds queens/q8 5 0. .5 10000 100 .000001 10 test`  
For a reminder of these parameters, simply run `./ds` without any arguments. The output files ending `.err`, `.sol`, and `.stats` contain, respectively, a log of the errors and variable states throughout the search, the vertex numbers in the dominating set (if found), and statistics on the performance of the algorithm.  

## Network File Format
The net file should have the following format:  
`number of nodes`  
(empty line)  
`0` (SPACE) `number of nodes connected to node 0`  
(list of nodes connected to node 0, separated by SPACE)  
`1` (SPACE) `number of nodes connected to node 1`  
(list of nodes connected to node 1, separated by SPACE)  
...  

As an example, the `queens` folder contains the net files for several sizes of the queen's graph: e.g., `q8` is the 8x8 queen's graph (i.e., a normal chessboard). The domination numbers (smallest possible size for a dominating set) of the n x n queen's graph are known for n up to 25 (see sequence A075458 on the OEIS). For n=8 the domination number is 5, meaning you can find a dominating set of size 5 but not any smaller. Try out the algorithm and see if you can find one! Or if you're feeling ambitious, try one of the larger graphs for which the domination number is not known and see how small a dominating set you can find.  
