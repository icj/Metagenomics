## Move from git directory to simulation working directory on HPC

system("cd ~/my_research/Metagenomics/")
system("mkdir ~/Metagenomics/")
system("mkdir ~/Metagenomics/betabinsim/")
system("mkdir ~/Metagenomics/negbinsim")
system("cp betabinsim.R ~/Metagenomics/betabinsim/")
system("cp negbinsim.R ~/Metagenomics/negbinsim")
system("cp metastats.R ~/Metagenomics/")

##====Beta Binomial Simulations====

sim.ids <- paste("betabin", 1:36, sep = "")

##----Create Directories----
for(h in 1:length(sim.ids)) system(paste("mkdir ~/Metagenomics/betabinsim/", sim.ids[h], sep = ""))


##----Simulation Settings Files----
n <- c(5, 10, 25)
mu <- c(.2, .1, .05, .01)
a <- c(4, 2, 1.5)
l <- 0
for(i in n){
    for(j in mu){
        for(k in a){
            l <- l + 1
            cat("##====Simulation Settings====", "\n",
                "a  <- ", k, "\n", 
                "mu <- ", j, "\n", 
                "n  <- ", i, "\n",
                "sim.id <- ", "\"", paste("betabin", l, sep = ""), "\"", "\n", 
                sep = "", 
                file = paste(paste("~/Metagenomics/betabinsim/", sim.ids[l], sep = ""),
                             paste(paste("/betabin", l, sep = ""), "settings.R", sep = ""), sep = "")
            )
        }
    }
}


##----CSH Files----
m <- 1
for(m in 1:length(sim.ids)){
    cat("#!/bin/csh", "\n",
        paste("#PBS -N", sim.ids[m], sep = " "), "\n",
        "#PBS -m bea", "\n",
        "#PBS -M icj@email.arizona.edu", "\n",
        "#PBS -W group_list=anling", "\n",
        "#PBS -q standard", "\n",
        "#PBS -l jobtype=serial", "\n",
        "#PBS -l select=1:ncpus=1:mem=1gb", "\n",
        "#PBS -l place=pack:shared", "\n",
        "#PBS -l walltime=05:00:00", "\n",
        "#PBS -l cput=05:00:00", "\n",
        "\n",
        "### Load R", "\n",
        "source /usr/share/Modules/init/csh", "\n",
        "module load R/2.15.2", "\n",
        "\n",
        paste("cd ~icj/Metagenomics/betabinsim/", sim.ids[m], sep = ""), "\n",
        "### clean the workspace if needed", "\n",
        "rm -f .RDATA", "\n",
        "\n",
        "### combine the sim settings with the generic sim file", "\n",
        paste(paste(paste("cat", paste(sim.ids[m], "settings.R", sep = ""), sep = " "),
                    "~/Metagenomics/betabinsim/betabinsim.R >", sep = " "), 
              paste(sim.ids[m], "sim.R", sep = ""), sep = " "), "\n",
        "\n",
        "### run simulation", "\n",
        "date", "\n",
        paste("R CMD BATCH", paste(sim.ids[m], "sim.R", sep = ""), sep = " "), "\n",
        "date", "\n",
        
        sep = "", 
        file = paste(paste("~/Metagenomics/betabinsim/", sim.ids[m], sep = ""),
                     paste(paste("/betabin", m, sep = ""), ".csh", sep = ""), sep = "")
    )
}

##====Run Simulations====
v <- 1
for(v in 1:length(sim.ids)){
    system(paste(paste("cd", paste("~/Metagenomics/betabinsim/", sim.ids[v], sep = ""), sep = " "),
                 paste("qsub", paste(sim.ids[v], "csh", sep = "."), sep = " "), sep = " ; ")
    )
}

rm(list = ls()) # cleanup

##====Negative Binomial Simulations====

sim.ids <- paste("negbin", 1:36, sep = "")

##----Create Directories----
for(h in 1:length(sim.ids)) system(paste("mkdir ~/Metagenomics/negbinsim/", sim.ids[h], sep = ""))


##----Simulation Settings Files----
n <- c(5, 10, 25)
mu <- c(.2, .1, .05, .01)
a <- c(4, 2, 1.5)
l <- 0
for(i in n){
    for(j in mu){
        for(k in a){
            l <- l + 1
            cat("##====Simulation Settings====", "\n",
                "a  <- ", k, "\n", 
                "mu <- ", j, "\n", 
                "n  <- ", i, "\n",
                "sim.id <- ", "\"", paste("negbin", l, sep = ""), "\"", "\n", 
                sep = "", 
                file = paste(paste("~/Metagenomics/negbinsim/", sim.ids[l], sep = ""),
                             paste(paste("/negbin", l, sep = ""), "settings.R", sep = ""), sep = "")
            )
        }
    }
}


##----CSH Files----
for(m in 1:length(sim.ids)){
    cat("#!/bin/csh", "\n",
        paste("#PBS -N", sim.ids[m], sep = " "), "\n",
        "#PBS -m bea", "\n",
        "#PBS -M icj@email.arizona.edu", "\n",
        "#PBS -W group_list=anling", "\n",
        "#PBS -q standard", "\n",
        "#PBS -l jobtype=serial", "\n",
        "#PBS -l select=1:ncpus=1:mem=1gb", "\n",
        "#PBS -l place=pack:shared", "\n",
        "#PBS -l walltime=05:00:00", "\n",
        "#PBS -l cput=05:00:00", "\n",
        "\n",
        "### Load R", "\n",
        "source /usr/share/Modules/init/csh", "\n",
        "module load R/2.15.2", "\n",
        "\n",
        paste("cd ~icj/Metagenomics/negbinsim/", sim.ids[m], sep = ""), "\n",
        "### clean the workspace if needed", "\n",
        "rm -f .RDATA", "\n",
        "\n",
        "### combine the sim settings with the generic sim file", "\n",
        paste(paste(paste("cat", paste(sim.ids[m], "settings.R", sep = ""), sep = " "),
                    "~/Metagenomics/negbinsim/negbinsim.R >", sep = " "), 
              paste(sim.ids[m], "sim.R", sep = ""), sep = " "), "\n",
        "\n",
        "### run simulation", "\n",
        "date", "\n",
        paste("R CMD BATCH", paste(sim.ids[m], "sim.R", sep = ""), sep = " "), "\n",
        "date", "\n",
        
        sep = "", 
        file = paste(paste("~/Metagenomics/negbinsim/", sim.ids[m], sep = ""),
                     paste(paste("/negbin", m, sep = ""), ".csh", sep = ""), sep = "")
    )
}

##====Run Simulations====
for(v in 1:length(sim.ids)){
    system(paste(paste("cd", paste("~/Metagenomics/negbinsim/", sim.ids[v], sep = ""), sep = " "),
                 paste("qsub", paste(sim.ids[v], "csh", sep = "."), sep = " "), sep = " ; ")
    )
}

