#!/usr/bin/env Rscript
# packages ####
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

# dev packages for flipping x axis of mutant
library(scales)
library(facetscales)

# constants ####
mRNA = "mRNA"
prot = "protein"
phos = "active protein"
wildtype = "Wildtype"
mutant = "KP1 Knockout"

KP_colors = c("#af70b6", "#654169", "#9f3b5a", "#9400d3")
TF_colors = c("#87ac33", "#50661e", "#2fa540", "#8FBC8F")
O_colors  = c("#79712f", "#DAA520")

# functions ####
readfunc = function(fnames) {
    
    meas = list()
    for (i in names(fnames)) {
        meas[[i]] = read.table(fnames[[i]], header=F, sep=" ")
    }
    
    # time along axis=1, node along axis=2
    for (i in c("r", "p", "psi")) {
        meas[[i]] = t(meas[[i]])
    }
    
    KP = cbind(meas[["r"]][,1:np], meas[["p"]][,1:np], meas[["psi"]][,1:np])
    TF = cbind(meas[["r"]][,np+1:nt], meas[["p"]][,np+1:nt], meas[["psi"]][,np+1:nt])
    if(no > 0) O = meas[["r"]][,np+nt+1:no]
    
    KP = data.frame(KP)
    TF = data.frame(TF)
    if(no > 0) O  = data.frame(O)
    
    colnames(KP) = c(paste0("KP",1:np," mRNA"), paste0("KP",1:np," prot"), paste0("KP",1:np," phos"))
    colnames(TF) = c(paste0("TF",1:nt," mRNA"), paste0("TF",1:nt," prot"), paste0("TF",1:nt," phos"))
    if(no > 0) colnames(O) = paste0("O",1:no," mRNA")
    
    # split as mRNA, prot, phos
    KP_ = list(mRNA=KP[,1:np], prot=KP[,np+1:np], phos=KP[,2*np+1:np])
    TF_ = list(mRNA=TF[,1:nt], prot=TF[,nt+1:nt], phos=TF[,2*nt+1:nt])
    
    # to remove curves always at zero
    rmz = function(df) {df[colSums(df) > 0]}
    
    # add time and remove curves that are always at zero
    for (i in c("mRNA", "prot", "phos")) {
        KP_[[i]] = rmz(KP_[[i]])
        TF_[[i]] = rmz(TF_[[i]])
        KP_[[i]]["time"] = meas[["t"]]
        TF_[[i]]["time"] = meas[["t"]]
    }
    if(no > 0) {
        O = rmz(O)
        O["time"] = meas[["t"]]
    }
    
    # melt
    mlt = function(df) {melt(df, id="time")}
    
    
    KP_mRNA=mlt(KP_[["mRNA"]]); KP_mRNA["measure"] = mRNA; KP_mRNA["type"] = "KP"
    KP_prot=mlt(KP_[["prot"]]); KP_prot["measure"] = prot; KP_prot["type"] = "KP"
    KP_phos=mlt(KP_[["phos"]]); KP_phos["measure"] = phos; KP_phos["type"] = "KP"
    TF_mRNA=mlt(TF_[["mRNA"]]); TF_mRNA["measure"] = mRNA; TF_mRNA["type"] = "TF"
    TF_prot=mlt(TF_[["prot"]]); TF_prot["measure"] = prot; TF_prot["type"] = "TF"
    TF_phos=mlt(TF_[["phos"]]); TF_phos["measure"] = phos; TF_phos["type"] = "TF"
    if(no > 0) {
        O_mRNA=mlt(O); O_mRNA["measure"] = mRNA; O_mRNA["type"] = "O"
        return(rbind(KP_mRNA, KP_prot, KP_phos, TF_mRNA, TF_prot, TF_phos, O_mRNA))
    }
    else return(rbind(KP_mRNA, KP_prot, KP_phos, TF_mRNA, TF_prot, TF_phos))
}

# settings ####
# mut = 1; outfname = "simulation_KP1.pdf"
# setwd("~/cwd/data/testdata/data_cas"); nt = 3; np = 3; no = 1
# setwd("~/cwd/data/testdata/pres18c"); nt = 3; np = 3; no = 0
# setwd("~/cwd/data/testdata/simi"); nt = 2; np = 3; no = 1
# setwd("~/cwd/data/toy/simulation"); nt = 2; np = 3; no = 1
args = commandArgs(T)
if(length(args) >= 3) {
    np = args[1]
    nt = args[2]
    no = args[3]
}
if(length(args) == 4) {
    mut = args[4]
    outfname = paste0("simulation_", mut, ".pdf")
} else {
    outfname = "simulation.pdf"
}


n = np+nt+no
wt_fnames = list(
    r = "sim_r.mat",
    p = "sim_p.mat",
    psi = "sim_psi.mat",
    t = "sim_t.mat"
)
mut = 1
mut_fnames = list(
    r = paste0("sim_r_", mut, ".mat"),
    p = paste0("sim_p_", mut, ".mat"),
    psi = paste0("sim_psi_", mut, ".mat"),
    t = paste0("sim_t_", mut, ".mat")
)
# process ####

wt  = readfunc(wt_fnames);   wt["experiment"] = wildtype
mut = readfunc(mut_fnames); mut["experiment"] = mutant
data = rbind(wt, mut)
# explicit order using factors
data$experiment = factor(data$experiment, levels=c(wildtype, mutant))
data$measure = factor(data$measure, levels=c(mRNA, prot, phos))
data$node = gsub(" .*", "", as.character(data$variable))
# be they are already sorted correctly so unique will provide them in correct KP->TF->O order
data$node = factor(data$node, levels=unique(data$node))
data$type = factor(data$type, levels=unique(data$type))

# plot ####

p = ggplot(data, aes(x=time, y=value, group=variable, color=node)) +
    geom_line(aes(linetype=measure, size=measure)) +
    scale_linetype_manual(values=c("solid", "solid", "dashed")) +
    scale_size_manual(values=c(1.2, .5, .5)) +
    scale_color_manual(values=c(P_colors[1:np], T_colors[1:nt], O_colors[1:no])) +
    scale_x_continuous(expand=c(0,0)) +
    xlab("time [min]") +
    ylab("nondimensionalized concentration") + 
    ggtitle("Steady state simulation") +
    guides(color = guide_legend(override.aes = list(size=8))) +
    theme_linedraw() +
    theme(panel.grid.major = element_line(colour = "gray"), 
          panel.grid.minor = element_line(colour = "lightgray"))


# without flipping x-axis of mutant
p = p + facet_grid(vars(type), vars(experiment), scales="free_x")
# or with flipping x-axis of mutant
# x_scales = list(Wildtype=scale_x_continuous(), Mutant=scale_x_reverse())
# p = p + facet_grid_sc(vars(type), vars(experiment), scales=list(x=x_scales))

ggsave("~/cwd/data/toy/simulation/simulation_KP1.pdf", p, width=7, height=5)


#####



