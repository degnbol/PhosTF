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
phos = "phos. protein"
wildtype = "Wildtype"
mutant = "Mutant"

P_colors = c("#af70b6", "#654169", "#9f3b5a")
T_colors = c("#87ac33", "#50661e", "#2fa540")
X_colors = c("#79712f")

# functions ####
readfunc = function(fnames) {
    
    meas = list()
    for (i in names(fnames)) {
        meas[[i]] = read.table(fnames[[i]], header=F, sep=" ")
    }
    
    # time along axis=1, node along axis=2
    for (i in c("r", "p", "phi")) {
        meas[[i]] = t(meas[[i]])
    }
    
    P  = cbind(meas[["r"]][,1:np], meas[["p"]][,1:np], meas[["phi"]][,1:np])
    TF = cbind(meas[["r"]][,np+1:nt], meas[["p"]][,np+1:nt], meas[["phi"]][,np+1:nt])
    if(nx > 0) X = meas[["r"]][,np+nt+1:nx]
    
    P  = data.frame(P)
    TF = data.frame(TF)
    if(nx > 0) X  = data.frame(X)
    
    colnames(P)  = c(paste0("P",1:np," mRNA"), paste0("P",1:np," prot"), paste0("P",1:np," phos"))
    colnames(TF) = c(paste0("T",1:nt," mRNA"), paste0("T",1:nt," prot"), paste0("T",1:nt," phos"))
    if(nx > 0) colnames(X) = paste0("X",1:nx," mRNA")
    
    # split as mRNA, prot, phos
    P_ = list(mRNA= P[,1:np], prot= P[,np+1:np], phos= P[,2*np+1:np])
    T_ = list(mRNA=TF[,1:nt], prot=TF[,nt+1:nt], phos=TF[,2*nt+1:nt])
    
    # to remove curves always at zero
    rmz = function(df) {df[colSums(df) > 0]}
    
    # add time and remove curves that are always at zero
    for (i in c("mRNA", "prot", "phos")) {
        P_[[i]] = rmz(P_[[i]])
        T_[[i]] = rmz(T_[[i]])
        P_[[i]]["time"] = meas[["t"]]
        T_[[i]]["time"] = meas[["t"]]
    }
    if(nx > 0) {
        X = rmz(X)
        X["time"] = meas[["t"]]
    }
    
    # melt
    mlt = function(df) {melt(df, id="time")}
    
    
    P_mRNA=mlt(P_[["mRNA"]]); P_mRNA["measure"] = mRNA; P_mRNA["type"] = "PK/PP"
    P_prot=mlt(P_[["prot"]]); P_prot["measure"] = prot; P_prot["type"] = "PK/PP"
    P_phos=mlt(P_[["phos"]]); P_phos["measure"] = phos; P_phos["type"] = "PK/PP"
    T_mRNA=mlt(T_[["mRNA"]]); T_mRNA["measure"] = mRNA; T_mRNA["type"] = "TF"
    T_prot=mlt(T_[["prot"]]); T_prot["measure"] = prot; T_prot["type"] = "TF"
    T_phos=mlt(T_[["phos"]]); T_phos["measure"] = phos; T_phos["type"] = "TF"
    if(nx > 0) {
        X_mRNA=mlt(X); X_mRNA["measure"] = mRNA; X_mRNA["type"] = "X"
        return(rbind(P_mRNA, P_prot, P_phos, T_mRNA, T_prot, T_phos, X_mRNA))
    }
    else return(rbind(P_mRNA, P_prot, P_phos, T_mRNA, T_prot, T_phos))
}

# settings ####
# setwd("/Users/christian/GoogleDrev/PKTFX/testdata/data_cas")
setwd("/Users/christian/GoogleDrev/PKTFX/testdata/pres18c")
nt = 3; np = 3; nx = 0; n = nt+np+nx
wt_fnames = list(
    r = "sim_r.mat",
    p = "sim_p.mat",
    phi = "sim_phi.mat",
    t = "sim_t.mat"
)
mut = 2
mut_fnames = list(
    r = paste0("sim_r_", mut, ".mat"),
    p = paste0("sim_p_", mut, ".mat"),
    phi = paste0("sim_phi_", mut, ".mat"),
    t = paste0("sim_t_", mut, ".mat")
)
# process ####

wt  = readfunc(wt_fnames);   wt["experiment"] = wildtype
mut = readfunc(mut_fnames); mut["experiment"] = mutant
data = rbind(wt, mut)
# explicit order using factors
data$experiment = factor(data$experiment, levels=c(wildtype, mutant))
data$measure = factor(data$measure, levels=c(mRNA, prot, phos))
data$node = sapply(as.character(data$variable), function(x) strsplit(x, " ")[[1]][1])


# plot ####

p = ggplot(data, aes(x=time, y=value, group=variable, color=node)) +
    geom_line(aes(linetype=measure, size=measure)) +
    scale_linetype_manual(values=c("solid", "solid", "dashed")) +
    scale_size_manual(values=c(.5, 1.2, .5)) +
    scale_color_manual(values=c(P_colors[1:np], T_colors[1:nt], X_colors[1:nx])) +
    xlab("time [min]") +
    ylab("nondimensionalized concentration") + 
    ggtitle("Steady state simulation") +
    guides(color = guide_legend(override.aes = list(size=8))) +
    theme_linedraw() +
    theme(panel.grid.major = element_line(colour = "gray"), 
          panel.grid.minor = element_line(colour = "lightgray"))


# without flipping x-axis of mutant
# p + facet_grid(vars(type), vars(experiment), scales="free_x")
# or with flipping x-axis of mutant
x_scales = list(Wildtype=scale_x_continuous(), Mutant=scale_x_reverse())
p + facet_grid_sc(vars(type), vars(experiment), scales=list(x=x_scales))




#####



