# chromnum

#### Authors: Abhirup Bandyopadhyay, Amit Kumar Dhar and Sankar Basu *
#### * Author for all correspondence

#### (c) Indian Copyright filed, Diary Number. 4955/2018-CO/SW, Filed. April, 2018

#### Funding: The work was supported by the Department of Science and Technology – Science and Engineering Research Board (DST-SERB research grant PDF/2015/001079)

#### Reference:   https://link.springer.com/article/10.1007%2Fs00500-019-04278-8

##### Cite: Bandyopadhyay, A.; Dhar, A.K. Basu, S.; Graph coloring: a novel heuristic based on trailing path — properties, perspective and applications in structured networks. Soft Computing, 2019, https://doi.org/10.1007/s00500-019-04278-8


Requires: `Octave` or `MATLAB (version R2016a or higher)`  

Installation
```
git clone https://github.com/nemo8130/Chromnum
cd Chromnum
run the script chromnum.m / chromnum_octave.m in MATLAB command window / chromnum_octave.m in Octave as a function with a single input argument
which is the filename (with full path / locally from the current directory) containing the adjacency matrix
```

```
% Chromnum computes the chromatic number of a graph 
% Chromatic Number of a graph is the minimum number of colors by whcih 
% All nodes of the graph could be exhaustively mapped so that no two adjacent nodes share the same color. 
% Chromnum implements a novel trailing path meta-heuristic approximation algorithm
% The function takes the 1-0 adjacency matrix stored in a text file as its sole input argument
% And returns the chromatic number of the corresponding graph
```

**Usage: crn = chromnum('example/adj.inp')**
**Usage: crn = chromnum_octave('example/adj.inp')**

**RUN: help chromnum or help chromnum_octave to get familiarize with the different options:: The default mode is Graphical. But the user can also turn-off the graphics (in case of say, batch jobs) and also controll the number of finite iterations (particularly helpful for higher dimensional regular graphs) 


```
% It also generates one color map for the graph and tabulates the same 
% Please note that this color map could be degenerate but the chromatic
% number should be identical 
% Furthermore it also returns the order at which the nodes have been
% colored. Again this may have degenerate solutions.
%
% Example adjacency matrix file :
% (Each row should contain the adjacencies of one node in the graph)
% (The matrix entries should be separated by a single white-space)
%
% 0 1 1 0 1
% 1 0 0 1 1
% 1 0 0 0 1
% 0 1 0 0 1
% 1 1 1 1 0
% 
% OUTPUT:
% Number of Nodes in the given undirected graph : 5
% Number of edges : 7 
% Chromatic Number : 3 

% Proposed Colormap:
% ------------------
% Node-1 : color-3
% Node-2 : color-2
% Node-3 : color-2
% Node-4 : color-3
% Node-5 : color-1
%


