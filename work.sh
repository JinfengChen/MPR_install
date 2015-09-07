echo "use a new dir for all MPR-related package"
export R_LIBS=$R_LIBS:"/bigdata/stajichlab/cjinfeng/software/R_package_MPR/"

echo "install Biocgeneric"
#R.3.2.0
R CMD INSTALL BiocGenerics_0.10.0.tar.gz

echo "install Biobase"
#R.3.2.0, update all other package required, yes
R
source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")

echo "install MPR"
#R.3.2.0
cd MPR_genotyping
R CMD INSTALL MPR_NAMESPACE_0.1.tar.gz
qsub MPR.sh


