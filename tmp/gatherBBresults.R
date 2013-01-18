##====Gather and load Results====
## run locally

sim.ids <- paste("betabin", 1:36, sep = "")

## IF option 1 already ran:
load("~/Dropbox/Metagenomics/bbResults.RData")

##----Option 1: Get whole directory structure from HPC----
## Gather all files from the HPC and put in my Dropbox
system("scp -r icj@login.hpc.arizona.edu:~icj/Metagenomics/betabinsim/"*" ~/Dropbox/Metagenomics/betabinsim/")

## Load RData files from Dropbox
for(w in 1:length(sim.ids)){
    load(paste(paste(paste("~/Dropbox/Metagenomics/betabinsim/", sim.ids[w], sep = ""), 
                     sim.ids[w], sep = "/"), "RData", sep = "."))
}

## Load Pvalue files
for(w in 1:length(sim.ids)){
    aa <- read.delim(paste(paste(paste("~/Dropbox/Metagenomics/betabinsim/", sim.ids[w], sep = ""), sim.ids[w], sep = "/"), "pvalues.txt", sep = ""))
    assign(paste(sim.ids[w], "pvalues", sep = "."), as.matrix(aa))
}

rm(list = c("aa", "w"))

## Save workspace containing all BetaBinomial results
save.image("~/Dropbox/Metagenomics/bbResults.RData")

##----Option 2: copy results from HPC only (Option 1 is better)
system("mkdir ~/Dropbox/Metagenomics/betabinsim/results/")
## RData files
for(w in 1:length(sim.ids)){
    system(paste(
        paste(
            paste(
                paste(
                    paste("scp icj@login.hpc.arizona.edu:~icj/Metagenomics/betabinsim/", sim.ids[w], sep = ""), "/", sep =""),
                        sim.ids[w], sep = ""),
                            "RData", sep = "."),
                                "~/Dropbox/Metagenomics/betabinsim/results/", sep = " "
        )
    )
}

## Pvalue text files
for(w in 1:length(sim.ids)){
    system(paste(
        paste(
            paste(
                paste(
                    paste("scp icj@login.hpc.arizona.edu:~icj/Metagenomics/betabinsim/", sim.ids[w], sep = ""), "/", sep =""),
                        sim.ids[w], sep = ""),
                            "pvalues.txt", sep = ""),
                                "~/Dropbox/Metagenomics/betabinsim/results/", sep = " "
        )
    )
}


## Load RData files
for(w in 1:length(sim.ids)){
    load(paste(paste("~/Dropbox/Metagenomics/betabinsim/results/", sim.ids[w], sep = ""), "RData", sep = "."))
}

## Load Pvalue files
for(w in 1:length(sim.ids)){
    aa <- read.delim(paste(paste("~/Dropbox/Metagenomics/betabinsim/results/", sim.ids[w], sep = ""), "pvalues.txt", sep = ""))
    assign(paste(sim.ids[w], "pvalues", sep = "."), as.matrix(aa))
}

rm(list = c("aa", "w"))
       