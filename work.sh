echo "use a new dir for all MPR-related package"
export R_LIBS=$R_LIBS:"/rhome/cjinfeng/software/tools/R_package_MPR/"

echo "install Biocgeneric"
/rhome/cjinfeng/software/tools/R-2.15.3/bin/R CMD INSTALL BiocGenerics_0.10.0.tar.gz

echo "install Biobase"
/rhome/cjinfeng/software/tools/R-2.15.3/bin/R
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")

echo "install MPR"
/rhome/cjinfeng/software/tools/R-2.15.3/bin/R
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("annotate")
biocLite("AnnotationDbi")

cd MPR_genotyping
/rhome/cjinfeng/software/tools/R-2.15.3/bin/R CMD INSTALL MPR_0.1.tar.gz
qsub MPR.sh


