
#INSTALL bGMYC bGMYC_1.0.2.tar.gz:: FROM ARCHIVE (I.E. FOLDER) AND NOT CRAN REPOS.
bGMYC



library(ape)
library("bGMYC", lib.loc="~/R/win-library/3.0")
setwd("C:\\Users\\engelbrechth\\Documents\\Since May 2016_17-06-2016\\RAXML\\CROT\\NEW ALIGNMENTS\\crot5_all_cod_Dasypel-16S-EXCL_21-02-2017_1000reps")


raxml.tree = read.nexus("nexustree.tre") #
plot(raxml.tree)
raxml.tree

raxml.tree=ladderize(raxml.tree, right = TRUE)
plot(raxml.tree)

tip=c("T_pul_CAS220642","D_fas_YPM13202", "T_sem", "B_forst", "D_unicolo")# this is for trimming unwanted tips (in my case the outgroup)
tip


PRUNE=drop.tip(raxml.tree, tip, trim.internal = TRUE, subtree = FALSE,
              root.edge = 0, rooted = is.rooted(raxml.tree)) # goes with the above...
plot(PRUNE)
PRUNE


is.ultrametric(PRUNE)
utree = chronos(PRUNE, lambda = 0, model = "correlated")

is.ultrametric(utree)
plot(utree)


library(phytools)
write.tree(NOOG, "pruned NOOG.tre")# not needed really..

bgmyc.dataprep(utree)


#bgmyc.singlephy(raxml.tree[[1]], mcmc=50000, burnin=1, thinning=10, t1=2, t2=100, start=c(1,1,25))->result.single #from "instructions" text file and will need to change according your tree

bgmyc.singlephy(utree, mcmc=1000000, burnin=1, thinning=10, t1=2, t2=87, start=c(1,1,25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)->result.single #adapted from "instructions" text file as well as pacakge cmd line
#second try mcmc=1mil instead of 50k and burnin also adjusted

bgmyc.singlephy(phylo, mcmc, burnin, thinning, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 50), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)

plot(result.single)# from website

spec.probmat(result.single)->result.probmat

plot(result.probmat, utree)

bgmycrates=checkrates(result.single)

plot(bgmycrates)

