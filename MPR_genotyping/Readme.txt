The tar file MPR_0.1.tar.gz (or zip file MPR_0.1.zip) is an R package named 
MPR. 

The package is built for implementing algorithms based on the principle of 
maximum parsimony of recombination (MPR) in a mapping population to infer 
parental genotypes and genotype the population in a parent-independent manner.

To install the package, please install R software before. Additional R 
packages genefilter and qtl are needed.

For linux systems, type:

R CMD INSTALL MPR_0.1.tar.gz

Or install the package (MPR_0.1.zip) from menu in windows.

Main algorithms/functions implemented in the MPR package are:

1. localMPR. The core function to infer parental genotypes (MPR inference) in 
a local chromosome region by minimizing the number of recombination events in 
the population. Several factors may affect the accuracy of MPR inference: 
(a). The number of SNPs processed each time (window size);
(b). The density of putative SNPs (the distance between SNP sites);
(c). The maximum step size of the heuristic perturbation (the parameter of 
maxNStep);
(d). The number of RILs.


2. globalMPRByMarkers. The function to do MPR inference in whole chromosome 
by using localMPR to infer parental genotypes in hundreds of local regions 
and assemble them aiding with low-coverage sequences of one parent or known 
markers.

3. globalMPRRefine. The function to refine SNPs by resampling and Bayesian 
inference based on results of globalMPRByMarkers. 

SNP alleles at all putative SNP sites identified from the population (snpData
), recovering SNP alleles from low-depth sequences of one parent (markerData)
, postitive (checkData) and negative (fSNP) SNP datasets are stored in data 
directory and can be loaded into R.
 
See MPR_example.R for usages of the package.

Thanks for your interest in our works.

Weibo Xie <xwbcn@webmail.hzau.edu.cn>


