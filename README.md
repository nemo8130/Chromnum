# chromnum
Requires: `Octave` or `MATLAB (version R2016a or higher)`  

#### The program is not yet made public due to security reasons while it is being reviewed. Please send an email request to nemo8130@gmail.com for the time being if you wish to download and use it. Please mention your identity and describe your interest. We will get back to you with a licence agreement and the source code.

Installation
```
git clone https://github.com/nemo8130/Chromnum
cd Chromnum
run the script chromnum.m / chromnum_octave.m in MATLAB command window / chromnum_octvae.m in Octave as a function with a single input argument
which is the filename (with full path / locally from the current directory) containing the adjacency matrix
```

```
% Chromnum computes the chromatic number of a graph 
% Chromatic Number of a graph is the minimum number of colors by whcih 
% All nodes of the graph could be exhaustively mapped so that no two
% adjacent nodes share the same color. 
% Chromnum implements the "trailpathSA" algorithm
% The function takes the 1-0 adjacency matrix stored in a text file as its sole input argument
% And returns the chromatic number of the corresponding graph
```

**Usage: crn = chromnum('adj.inp')**
**Usage: crn = chromnum_octave('adj.inp')**

```
% It also generates one color map for the graph and tabulates the same 
% Please note that this color map could be degenerate but the chromatic
% number should be identical 
% Furthermore it also returns the order at which the nodes have been
% colored. Again this may have degenerate solutions.
%
% Example adjacency matrix file :
% (Each row should contain the adjacencies of one node in the graph)
% (The matrix entries should be seperated by a single white-space)
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

