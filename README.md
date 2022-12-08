# Dominating Sets
This repository demonstrates the use of iterative projection methods to solve a classic discrete problem: finding dominating sets. For a given network, a dominating set is a subset of the network such that every node in the complete network is either in the dominating subset or shares an edge with a member of the subset. The algorithm in `ds.c` is one way of finding dominating sets. See arXiv:2106.05206 for details.  

To compile: `gcc -O2 ds.c -lm -o ds`  
To run: `./ds netfile D epsilon beta maxiter iterstride stoperr trials id`  
- `netfile`: the name of the file containing the network connectivity  
- `D`: the maximum cardinality of the dominating set to look for  
- `epsilon`: determines how quickly to tune the metric parameter; generally epsilon≤100/maxiter is safe (assuming maxiter>>100); set epsilon=0 for no tuning  
- `beta`: another hyperparameter; values less than 1 work best; start with 0.5  
- `maxiter`: the maximum number of iterations for the algorithm to make before giving up  
- `iterstride`: the number of iterations to skip between each print to the .err file (must be at least 1 and at most equal to maxiter)  
- `stoperr`: the error at which the algorithm will stop and declare it has found a solution; try a small nonzero value like 0.001  
- `name`: the name for the output files (name.err, name.sol, name.stats)  

An example run: `./ds queens/q8 5 0. .5 10000 100 .001 10 test`  
For a reminder of these parameters, simply run `./ds` without any arguments. The program generates output files ending .err, .sol, and .stats that contain, respectively, a log of the errors and variable states throughout the search, the vertex numbers in the dominating set (if found), and statistics on the performance of the algorithm.  

## Network File Format
The net file should have the following format:  
`number of nodes`  
(empty line)  
`0` (SPACE) `number of nodes connected to node 0`  
(list of nodes connected to node 0, separated by SPACE)  
`1` (SPACE) `number of nodes connected to node 1`  
(list of nodes connected to node 1, separated by SPACE)  
...  

As an example, `q8` in the `queens` folder is the net file for the 8x8 queen's graph.
