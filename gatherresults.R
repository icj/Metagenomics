##====Gather Results====
## make result directories in git repository
system("mkdir ~/my_research/Metagenomics/betabinsim/")
system("mkdir ~/my_research/Metagenomics/negbinsim/")

## Gather all files from the HPC and put in my git repository
system("scp -r icj@login.hpc.arizona.edu:~icj/Metagenomics/betabinsim/\"*\" ~/my_research/Metagenomics/betabinsim/")
system("scp -r icj@login.hpc.arizona.edu:~icj/Metagenomics/negbinsim/\"*\" ~/my_research/Metagenomics/negbinsim/")

## Load RData files
rm(list = ls())
bb.sim.ids <- paste("betabin", 1:36, sep = "")
for(w in 1:length(bb.sim.ids)){
    load(paste(paste(paste("~/my_research/Metagenomics/betabinsim/", bb.sim.ids[w], sep = ""), 
                     bb.sim.ids[w], sep = "/"), "RData", sep = "."))
}
nb.sim.ids <- paste("negbin", 1:36, sep = "")
for(w in 1:length(nb.sim.ids)){
    load(paste(paste(paste("~/my_research/Metagenomics/negbinsim/", nb.sim.ids[w], sep = ""), 
                     nb.sim.ids[w], sep = "/"), "RData", sep = "."))
}

## Load Pvalue files
for(w in 1:length(bb.sim.ids)){
    aa <- read.delim(paste(paste(paste("~/my_research/Metagenomics/betabinsim/", bb.sim.ids[w], sep = ""), 
                                 bb.sim.ids[w], sep = "/"), "pvalues.txt", sep = ""))
    assign(paste(bb.sim.ids[w], "pvalues", sep = "."), as.matrix(aa))
}
for(w in 1:length(nb.sim.ids)){
    aa <- read.delim(paste(paste(paste("~/my_research/Metagenomics/negbinsim/", nb.sim.ids[w], sep = ""), 
                                 nb.sim.ids[w], sep = "/"), "pvalues.txt", sep = ""))
    assign(paste(nb.sim.ids[w], "pvalues", sep = "."), as.matrix(aa))
}

rm(list = c("aa", "w"))

## Save workspace containing all results
save.image("~/my_research/Metagenomics/bbResults.RData")

