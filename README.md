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

col1="tr ";

col2=*Trace number*;

col3=*Length*;

col4=*Quality score*;

col5=*Number of peaks*;

col6=*Cluster number*;

col7=*Cluster number without filters*;

col8=*Distance from the cluster center*;

In the very end of the output file, information about the clusters is printed. Each cluster is descirbed in terms of cluster center, size and members (cluster core traces).

The values of the input parameters in **input.txt** are the default ones.

**Results**

Dataset Mixed contains 4 groups of traces corresponding to the unfolding of 4 different proteins:
* Group 1 - unfolding of the CNGA1 channel (Maity et al, Nat Commun, 2015)
* Group 2 - unfolding of a tandem globular polyprotein (Alpha3D + 4xNUG2) (Heenan and Perkins, Biophys J, 2018)
* Group 3 and Group 4 correspond to the unfolding of two unidentified proteins from AFM-SMFS experiments in the plasma membrane of the rod outer segment.

The total number of clusters you should obtain for Dataset Mixed is 5. 

The traces from the 4 groups in the dataset are distributed in our clusters as follows:
* cluster 1 - contains the traces from Group 2
* cluster 2 - contains the traces from Group 1
* cluster 3, 4, 5 - contain the traces from Group 3

The traces from group 4 do not belong to the core of any of the five clusters indicating that they are not similar enough according to the criteria implemented in the code. 

For more information you can check our reference paper.
