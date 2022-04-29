### this script uses the DSS package to determine DMRs with Bismark coverage files as input ###

# Jilke De Wilde
# 04-2022

usage <- function() {
  cat('Usage: Rscript DMRs_DSS.R  [options]  -g1 <groupList1> -g2 <groupList2> \
                                      -n1 <groupName1> -n2 <groupName2> \    
                                      -o <output>  \
                
    -g1 -g2 <input files>   comma-separated list of coverage files
    -n1 -n2 <group names>   names for the groups that should be compared (same order as input files)
    -o <output>   Output file listing DMR results
')
  q()
}

# get CL args
args <- commandArgs(trailingOnly=T)
i <- 1
while (i <= length(args)) {
  if (args[i] == '-h' || args[i] == '--help') {
      usage()
    } else if (i < length(args)) {
      if (args[i] == '-g1') {
        g1 <- strsplit(gsub("'","",args[i + 1])," ") 
      } else if (args[i] == '-g2') {
        g2 <- strsplit(gsub("'","",args[i + 1])," ")
      } else if (args[i] == '-n1') {
        n1 <- args[i + 1]
      } else if (args[i] == '-n2') {
        n2 <- args[i + 1]
      } else if (args[i] == '-o') {
        output <- args[i + 1]
      } else {
        cat('Error! Unknown parameter:', args[i], '\n')
        usage()
      }
      i <- i + 1
    } else {
      cat('Error! Unknown parameter with no arg:', args[i], '\n')
      usage()
    }
  i <- i + 1
}

# load libraries
library(DSS)
require(bsseq)

# load files and convert to correct input for DSS
prep_samples <- function(file){
  s <- read.table(file, col.names = c('chr','pos','stop','beta_perc','X','unmeth'))
  s$N <- s$X + s$unmeth
  s <- s[c('chr','pos','N','X')]
  return(s)
}

print("prepping samples")
g1.samples <- lapply(g1[[1]], prep_samples)
g2.samples <- lapply(g2[[1]], prep_samples)
grouped.samples <- c(g1.samples,g2.samples)

# name samples
name_samples <- function(file){
  n <- strsplit(file,"/")
  n <- n[[1]][length(n[[1]])]
  n <- gsub("_.*","",n)
  return(n)
}

names(g1.samples) <- lapply(g1[[1]], name_samples)
names(g2.samples) <- lapply(g2[[1]], name_samples)
name.lst <- c(names(g1.samples),names(g2.samples))
  
# make BSseq object
print(paste("making a BSseq object of",n1,"against",n2))

BSobj <- makeBSseqData(grouped.samples,name.lst)

# define amount of cores for parallelisation
mParam = MulticoreParam(workers=9, progressbar=TRUE)

# call DMLs
print(paste("calculating DMLs of",n1,"against",n2))

dmlTest <- DMLtest(BSobj, names(g1.samples), names(g2.samples), BPPARAM = mParam)

# call DMRs
print(paste("calculating DMRs of",n1,"against",n2))

dmrs <- callDMR(dmlTest, p.threshold = 0.05, delta = 0.05)

write.table(dmrs, file = output, sep='\t', col.names = TRUE)
