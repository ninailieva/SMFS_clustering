# SMFS_clustering

We provide the source code of our method for clustering heterogeneous AFM-SMFS data from pulling experiments performed on proteins together with the files neccesary to run a toy example.


The files included in this repository are:
* **cluster_traces.cpp** - the source code
* **input.txt** - contains the input parameters of the algorithm
* **Dataset_Mixed.txt** - test dataset containing four groups of traces corresponding to the unfolding of four different proteins
* **reference_output.txt** - the correct output results for Dataset Mixed

The code has been tested with an Intel compiler: intel/18.0.3 with the following command:

icc -std=c++0x -O3 cluster_traces.cpp -fopenmp -o cluster_traces.x

To run the code:

./cluster_traces.x input.txt > output.txt

The generated output file contains a list of the used parameter values followed by the actual results divided in columns containing the following information:

col1=*Trace number*;

col2=*Length*;

col3=*Quality score*;

col4=*Number of peaks*;

col5=*Cluster number*;

col6=*Distance from the cluster center*;

In the very end of the output file, information about the clusters is printed consistent with the one in Table 2 in our reference paper.

The values of the input parameters in **input.txt** are the default ones.

**Results**

Dataset Mixed contains 4 groups of traces corresponding to the unfolding of 4 different proteins:
* Group 1 - unfolding of the CNGA1 channel (Maity et al, Nat Commun, 2015)
* Group 2 - unfolding of a tandem globular polyprotein (Alpha3D + 4xNUG2) (Heenan and Perkins, Biophys J, 2018)
* Group 3 and Group 4 correspond to the unfolding of two unidentified proteins from AFM-SMFS experiments in the plasma membrane of the rod outer segment.

The total number of clusters you should obtain for Dataset Mixed is 7. Two of them: cluster number 5 and cluster number 7 contain one member only: the cluster center. We exclude these clusters from the analysis and 5 clusters remain.

The traces from the 4 groups in the dataset are distributed in our clusters as follows:
* cluster 1 - contains the traces from Group 2
* cluster 2 - contains the traces from Group 1
* cluster 3 - contains the traces from Group 4
* cluster 4 - contains the traces from Group 3
* cluster 5 - also contains traces from Group 3 

Both clusters 4 and 5 bear the same force pattern. Since the core of cluster 4 is more robust we consider cluster 4 to be the cluster matching the traces from Group 3. 
