Instructions to run hierBAPS with its bundles in Linux (Ubuntu 14.04.1 LTS tested)
Lu Cheng & Jukka Corander
16.01.2015

1. Use a text editor (gedit) to open the file "hierBAPS.sh"
set MCRROOT in line 16, that is the directory where you install MCR

2. Launch terminal, go to the folder containing hierBAPS files and type 

chmod u+x hierBAPS.sh
chmod u+x exData
chmod u+x hierBAPS
chmod u+x drawSnpMat
./hierBAPS.sh exData seqs.fa fasta

This will produce an input file called "seqs.mat" from you alignment file. The example file seqs.fa and seqs.xls are provided with this readme file.

3. Type 

./hierBAPS.sh hierBAPS seqs.mat L maxK results  

This will launch hierBAPS. Clustering is performed with L levels in the hierarchy and maxK is the prior upper bound for number of clusters. As in BAPS, hierBAPS will estimate the maximum a posteriori partition (MAP) with the number of clusters in the interval 1 to maxK. hierBAPS will save an output file named results.mat (binary format) and a partition file "results.partition.txt", where each column represents the MAP partion of that layer.

4. (optional) Assume you did use

./hierBAPS.sh hierBAPS seqs.mat 2 20 results 

In the previous step. You can continue the clustering from the previous result to a deeper level (here L=4) by typing 

./hierBAPS.sh hierBAPS results.mat 4 20 results2nd


5. Type 

./hierBAPS.sh drawSnpMat results.mat shuffle

This will draw a SNP matrix with rows and columns shuffled as shown in BAPS 6.0 manual. Without the "shuffle" parameter, the columns will not be shuffled. Consecutive rows between horizontal black lines represent a 1st layer cluster.

A tab delimited file "figInfo.txt" will also be produced, where the third column "SampleID" is the index in the input fasta file.



