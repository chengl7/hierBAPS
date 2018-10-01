Instructions to run hierBAPS with its bundles in Windows (Windows 8 tested)
Lu Cheng & Jukka Corander
15.8.2012

1. Launch command prompt, go to the folder containing hierBAPS files and type 

exData.exe seqs.fa fasta (or exData.exe seqs.xls xls)

This will produce an input file called "seqs.mat" from you alignment file. The example file seqs.fa and seqs.xls are provided with this readme file.

2. Type 

hierBAPS.exe seqs.mat L maxK results  

This will launch hierBAPS. Clustering is performed with L levels in the hierarchy and maxK is the prior upper bound for number of clusters. As in BAPS, hierBAPS will estimate the maximum a posteriori partition (MAP) with the number of clusters in the interval 1 to maxK. hierBAPS will save an output file named results.mat (binary format) and a partition file "results.partition.txt", where each column represents the MAP partion of that layer.

3. (optional) Assume you did use

hierBAPS.exe seqs.mat 2 20 results 

In the previous step. You can continue the clustering from the previous result to a deeper level (here L=4) by typing 

hierBAPS results.mat 4 20 results2nd


4. Type 

drawSnpMat results.mat shuffle

This will draw a SNP matrix with rows and columns shuffled as shown in BAPS 6.0 manual. Without the "shuffle" parameter, the columns will not be shuffled. Consequtive rows between horizontal black lines represent a 1st layer cluster.

A tab delimited file "figInfo.txt" will also be produced.


