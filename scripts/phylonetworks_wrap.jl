# DESCRIPTION
# Wrapper script to execute PhyloNetworks analysis in julia programming language
#
# USAGE
# > julia phylonetworks_wrap.jl consensus.nwk raxmltrees.nwk outgroup_taxon
# where consensus.nwk is a single newick tree with the species topology,
# raxmltrees.nwk is a set of multiple (usually discordant) local newick topologies, e.g. windowed RAxML output
# and outgroup_taxon is a single outgroup taxon present in the sample set
# addprocs(35) may need to be tuned to the number of processors available on the compute node
#
# AUTHOR
# Jonas Lescroart
# Created 02DEC2021
# Added slope heuristic 27FEB23

# Multithreading and packages
using Distributed
addprocs(35)
@everywhere using PhyloNetworks
@everywhere using PhyloPlots
@everywhere using RCall

# Load input
consensus = readTopology(ARGS[1])
raxmlCF = readTrees2CF(ARGS[2], writeTab = false, writeSummary = false)
#outgroup = ARGS[3]

# Prefixes for output filenames
prefix_net0 = string(ARGS[2][begin:end-4], "_net0")
prefix_net1 = string(ARGS[2][begin:end-4], "_net1")
prefix_net2 = string(ARGS[2][begin:end-4], "_net2")
prefix_net3 = string(ARGS[2][begin:end-4], "_net3")

# PDF filenames
net0_pdf = string(prefix_net0, ".pdf")
net1_pdf = string(prefix_net1, ".pdf")
net2_pdf = string(prefix_net2, ".pdf")
net3_pdf = string(prefix_net3, ".pdf")
slope_pdf = string(ARGS[2][begin:end-4], "_slope_heuristic.pdf")

# Network estimation and print PDFs
net0 = snaq!(consensus, raxmlCF, hmax = 0, filename = prefix_net0, seed = 8, runs = 50)
#rootatnode!(net0, outgroup)

R"pdf"(net0_pdf, width = 20, height = 10);
plot(net0, :R, showGamma=true);
R"dev.off()";

net1 = snaq!(net0, raxmlCF, hmax = 1, filename = prefix_net1, seed = 9, runs = 50)
#rootatnode!(net1, outgroup)

R"pdf"(net1_pdf, width = 20, height = 10);
plot(net1, :R, showGamma=true);
R"dev.off()";

net2 = snaq!(net1, raxmlCF, hmax = 2, filename = prefix_net2, seed = 10, runs = 50)
#rootatnode!(net2, outgroup)

R"pdf"(net2_pdf, width = 20, height = 10);
plot(net2, :R, showGamma=true);
R"dev.off()";

net3 = snaq!(net2, raxmlCF, hmax = 3, filename = prefix_net3, seed = 11, runs = 50)
#rootatnode!(net3, outgroup)

R"pdf"(net3_pdf, width = 20, height = 10);
plot(net3, :R, showGamma=true);
R"dev.off()";

# Slope heuristic
scores = [net0.loglik, net1.loglik, net2.loglik, net3.loglik]
hmax = collect(0:3)

R"pdf"(slope_pdf, width = 5, height = 5);
R"plot"(hmax, scores, type="b", main="Slope heuristic", ylab="Negative log-likelihood", xlab="Hybrid nodes", col="blue");
R"dev.off()";

