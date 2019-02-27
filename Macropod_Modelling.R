library(RPANDA)
library(dplyr); library(gridExtra)
library(geiger); library(phytools)
library(parallel); library(pspline)
library(mvMORPH); library(cladeMode)
library(awtools); library(ggplot2); library(ggridges) # show_col(mpalette); or pie(rep(1, length(mpalette)), col = mpalette, labels = mpalette)


macro.data <- read.csv("/PATH/CrownHeight_Data.csv", header=T)
min.tree <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MinAges_Constrained_PF/Macro_MinAges_HKY_Constrained_CON.tre")
mean.tree <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MeanAges_Constrained/Macro_MeanAges_HKY_Constrained_CON.tre")
max.tree <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MaxAges_Constrained_PF/Macro_MaxAges_HKY_Constrained_CON_fixed.tre")
cp.min <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPMinAges_Constrained_PF/Macro_CP_MinAges_CON.tre")
cp.mean <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPEstAges_Constrained_PF/Macro_CP_EstAges_CON.tre")
cp.max <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPMaxAges_Constrained_PF/Macro_CP_MaxAges_CON.tre")

enviro.data <- read.csv("/Users/Ian/Desktop/Macropod_Dating/Andrea_S1.csv", header=T)
flux.data <- read.csv("/Users/Ian/Desktop/Macropod_Dating/Aeolian_Flux.csv", header=T)

#sp.means <- read.csv("/Users/Ian/Desktop/Macropod_Dating/Couzens&Prideaux2018_dryad_package/data/CrownHeight_spMEANS.csv", header=T)
sp.means <- read.csv("/Users/Ian/Desktop/Macropod_Dating/Couzens&Prideaux2018_dryad_package/data/CrownHeight_Macropodinae_spMEANS.csv", header=T)

# Trim tree and data down to overlapping taxa
overlaps <- intersect(min.tree$tip.label, unique(macro.data$Taxon))
trim.tree <- drop.tip(mean.tree, setdiff(mean.tree$tip.label, overlaps))
trim.data <- filter(macro.data, Taxon %in% overlaps)

# OR trim tree and data down to just Macropodinae
macros <- filter(macro.data, Higher_tax == "Macropodinae")
overlaps <- intersect(min.tree$tip.label, unique(macros$Taxon))
trim.tree <- drop.tip(mean.tree, setdiff(mean.tree$tip.label, overlaps))
trim.data <- filter(macros, Taxon %in% overlaps)


# create a tibble to get the species means (ONLY DO THIS ONCE!)
#   sp.means <- trim.data %>%
#     group_by(Taxon) %>%
#     summarise_at(vars(H_HYPCD), mean)
#   # make sure to drop riparius, kaindensis, nubicola_2!
#   write.csv(sp.means, row.names=FALSE,
#             file="/Users/Ian/Desktop/Macropod_Dating/Couzens&Prideaux2018_dryad_package/data/CrownHeight_Macropodinae_spMEANS.csv")

species.means <- as.data.frame(sp.means$H_HYPCD); rownames(species.means) <- sp.means$Taxon
macro.means <- species.means[,1]; names(macro.means) <- rownames(species.means); name.check(trim.tree, macro.means)
data(InfTemp)


# define function that returns the SSE
multiFit <- function(x){
  resfit <- fit_t_env(trim.tree, macro.means, env_data=InfTemp, df=x, scale=F, plot=F, model="EnvExp")
  return(resfit$LH)
}

# Run optim to find span that gives min SSE, starting at 0.5
optim(par=c(0), multiFit, method="Nelder-Mead")


head(enviro.data)


# make a function to extrapolate data from a mean and confidence intervals (we're making extra points from the grass data)
time.estimates <- function (some.data, reps){
  time.col <- NULL; recon.col <- NULL
  for(j in 1:nrow(enviro.data)) {
    
    time.sd <- sd(c(enviro.data$Age[j], (enviro.data$Age[j] + enviro.data$Age_Error[j]), (enviro.data$Age[j] - enviro.data$Age_Error[j])))
    time.samples <- rnorm(reps, enviro.data$Age[j], time.sd)
    #time.samples <- runif(reps, min=(enviro.data$Age[j] - enviro.data$Age_Error[j]), max=(enviro.data$Age[j] + enviro.data$Age_Error[j]))
    time.col <- rbind(time.col, as.data.frame(time.samples))
    
    recon.sd <- sd(c(enviro.data$C4_recon_mean[j], enviro.data$C4_recon_lower[j], enviro.data$C4_recon_upper[j]))
    recon.samples <- rnorm(reps, enviro.data$C4_recon_mean[j], recon.sd/2)
    #recon.samples <- runif(reps, min=enviro.data$C4_recon_lower[j], max=enviro.data$C4_recon_upper[j])
    recon.col <- rbind(recon.col, as.data.frame(recon.samples))
  }
  time.est.samples <- cbind(time.col, recon.col)
  time.est.samples[time.est.samples<0] <- 0
  time.est.samples[time.est.samples>100] <- 100
  time.est.samples <- time.est.samples[order(time.est.samples$time.samples),]
  return(time.est.samples)
}

testo <- time.estimates(enviro.data, 100)
plot(testo)
testo[1,] <- c(0.1, 70)
extrap.data <- matrix(nrow=7, ncol=2, c(10,11,12,13,14,15,16,
                                        20,10,5,2.5,1.25,0.7,0)); colnames(extrap.data) <- colnames(testo)
testo <- rbind(testo, extrap.data)
plot(testo)

raw.data <- enviro.data[,c("Age", "C4_recon_mean")]

l1 <- loess(recon.samples ~ time.samples, data=testo, span=0.5)
l1.pred <- predict(l1)
plot(l1.pred)

lo <- loess(recon.samples ~ time.samples, data=testo, control=loess.control(surface="direct"))
l1.predo <- predict(lo, newdata=seq(max(testo$time.samples), max(nodeHeights(trim.tree)), 0.1))
plot(l1.predo)

test.pred <- append(l1.pred, l1.predo)
plot(test.pred)

seq(max(testo$time.samples), max(nodeHeights(trim.tree)), 0.5)

# Visualize the data in a few different ways to get an idea of what's going on
phenogram(trim.tree, macro.means)
plotTree.wBars(trim.tree, macro.means, args.plotTree=list(border=F, fsize=0.5)) # type="fan"
plotTree.barplot(trim.tree, macro.means, args.barplot=list(beside=TRUE, border=F))
dotTree(trim.tree, macro.means, ftype="i")

# I'm not sure what this function was for, but I can probably get rid of it
search.DF <- function(model, n.iter = 10, traits, n.proc = 8) {
  beginning <- Sys.time()
  
  init.params <- lapply(1:n.iter, function(x) {
    c(m0 = mean(traits), 
      v0 = runif(1, min = 1e-10, max = 0.01),
      d1 = runif(1, min=1e-10, max=0.1),
      d2 = -runif(1, min=1e-10, max=0.1),
      S = rnorm(1, 0, 0.25),
      sigma = runif(1, min=1e-10, max=0.1))
  })
  
  res.list <- mclapply(1:n.iter, function(x) {
    fitTipData(model, traits, GLSstyle=T, params0 = init.params[[x]])}, mc.cores = n.proc)
  
  for (j in 1:length())
  
  res.values <- unlist(lapply(res.list, function(x) x$value)) # make a vector of the values, so we can get the index number of the best
  res.values[which(abs(res.values) > abs(2 * median(res.values)))] <- NA # remove any nonsensical model fits from the options
  res.values[which((res.values) < median(res.values)/2)] <- NA
  best.res <- res.list[[which.min(res.values)]]
  
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time to fit the", model@name, "model from", n.iter, "starting points:", duration))
  
  return(list(all.results = res.list, best.result = best.res))
}

# Create a function that searches for the optimum smoothness of the trend by fitting a set of values
best.smoothing <- function (phy, trait.data, time.data=InfTemp, degrees=c(0,10,20,30,40,50), model="EnvExp", cores=6) {
  res.list <- mclapply(1:length(degrees), function(x) {
    fit_t_env(phy, trait.data, env_data=time.data, df=degrees[x], scale=F, plot=T, model=model)}, mc.cores = cores)
  for(i in 1:length(res.list)){res.list[[i]]$df <- degrees[i]}
  res.values <- unlist(lapply(res.list, function(x) x$aicc)) # make a vector of the values, so we can get the index number of the best
  best.res <- res.list[[which.min(res.values)]]
  
  plot(best.res, main=paste(model, "; AICc = ", round(best.res$aicc,2)), sub=paste("sigma = ",round(best.res$param[1],2), " beta = ",round(best.res$param[2],2), " df = ", best.res$df), col="red")
  
  return(list(all.results=res.list, best.result=best.res, best.df=best.res$df))
}

# Fit a number of models to the data (ENV, GRASS, BM, EB, Trend, Drift)
ENVexp <- best.smoothing(trim.tree, macro.means, time.data=InfTemp, degrees=c(10,20,30,40,50), model="EnvExp", cores=5)
ENVlin <- best.smoothing(trim.tree, macro.means, time.data=InfTemp, degrees=c(10,20,30,40,50), model="EnvLin", cores=5)

GRASSexp <- best.smoothing(trim.tree, macro.means, time.data=testo, degrees=c(40,50), model="EnvExp", cores=3)
#GRASSlin <- best.smoothing(trim.tree, macro.means, time.data=testo, degrees=c(30,40,50), model="EnvLin", cores=3)

FLUXexp <- best.smoothing(trim.tree, macro.means, time.data=flux.data, degrees=c(10,20,30,40,50), model="EnvExp", cores=5)
FLUXlin <- best.smoothing(trim.tree, macro.means, time.data=flux.data, degrees=c(10,20,30,40,50), model="EnvLin", cores=5)


BM_res <- fitContinuous(trim.tree, species.means, model="BM")
trend_res <- fitContinuous(trim.tree, species.means, model="trend")
EB_res <- fitContinuous(trim.tree, species.means, model="EB")
#drift_res <- fitContinuous(trim.tree, species.means, model="drift")


# Compare the models with AICc, and check differences across the trees
min_tree_FIT <- c(ENVexp$best.result$aicc, ENVlin$best.result$aicc, 
                  GRASSexp$best.result$aicc, FLUXexp$best.result$aicc, FLUXlin$best.result$aicc, #GRASSlin$best.result$aicc, 
                  BM_res$opt$aicc, trend_res$opt$aicc, EB_res$opt$aicc); names(min_tree_FIT) <- c("ENVexp", "ENVlin", "GRASSexp", "FLUXexp", "FLUXlin", "BM", "Trend", "EB")
aic.w(min_tree_FIT)
min.aic <- as.data.frame(as.vector(aic.w(min_tree_FIT))); min.aic$model <- names(min_tree_FIT); colnames(min.aic) <- c("aiccw", "model"); min.aic$age <- "Min_Ages"

mean_tree_FIT <- c(ENVexp$best.result$aicc, ENVlin$best.result$aicc, 
                  GRASSexp$best.result$aicc, FLUXexp$best.result$aicc, FLUXlin$best.result$aicc, #GRASSlin$best.result$aicc, 
                  BM_res$opt$aicc, trend_res$opt$aicc, EB_res$opt$aicc); names(mean_tree_FIT) <- c("ENVexp", "ENVlin", "GRASSexp", "FLUXexp", "FLUXlin", "BM", "Trend", "EB")
aic.w(mean_tree_FIT)
mean.aic <- as.data.frame(as.vector(aic.w(mean_tree_FIT))); mean.aic$model <- names(mean_tree_FIT); colnames(mean.aic) <- c("aiccw", "model"); mean.aic$age <- "Mean_Ages"

max_tree_FIT <- c(ENVexp$best.result$aicc, ENVlin$best.result$aicc, 
                  GRASSexp$best.result$aicc, FLUXexp$best.result$aicc, FLUXlin$best.result$aicc, #GRASSlin$best.result$aicc,
                  BM_res$opt$aicc, trend_res$opt$aicc, EB_res$opt$aicc); names(max_tree_FIT) <- c("ENVexp", "ENVlin", "GRASSexp", "FLUXexp", "FLUXlin", "BM", "Trend", "EB")
max.aic <- as.data.frame(as.vector(aic.w(max_tree_FIT))); max.aic$model <- names(max_tree_FIT); colnames(max.aic) <- c("aiccw", "model"); max.aic$age <- "Max_Ages"

all.aic <- rbind.data.frame(min.aic, mean.aic, max.aic); 
all.aic$age <- factor(all.aic$age, levels=c("Min_Ages", "Mean_Ages", "Max_Ages"))

all.aic$model.type <- c("ENV", "ENV", "GRASS", "FLUX", "FLUX", "BM", "Trend", "EB",
                        "ENV", "ENV", "GRASS", "FLUX", "FLUX", "BM", "Trend", "EB",
                        "ENV", "ENV", "GRASS", "FLUX", "FLUX", "BM", "Trend", "EB")
all.aic$model.type <- factor(all.aic$model.type, levels=c("ENV", "FLUX", "GRASS", "BM", "EB", "Trend"))

(ggplot(all.aic)
  + geom_bar(aes(y=aiccw, x=age, fill=model.type), stat="identity")
  + theme(axis.text.x=element_text(angle=25, hjust=1), panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual( values=wes_palette("Zissou1", 6, "continuous")))

# Plot what the different smoothnesses did to our data
grass.spline0  <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=0)
grass.spline10 <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=10)
grass.spline20 <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=20)
grass.spline30 <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=30)
grass.spline40 <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=40)
grass.spline50 <- sm.spline(x=testo$time.samples, y=testo$recon.samples, df=50)

plot(testo)
lines(grass.spline0, col="red", lwd=4)
lines(grass.spline10, col="green", lwd=4)
lines(grass.spline20, col="blue", lwd=4)
lines(grass.spline30, col="yellow", lwd=4)
lines(grass.spline40, col="violet", lwd=4)
lines(grass.spline50, col="gold", lwd=4)


env.spline0 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=0)
env.spline10 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=10)
env.spline20 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=20)
env.spline30 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=30)
env.spline40 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=40)
env.spline50 <- sm.spline(x=InfTemp$Age, y=InfTemp$Temperature, df=50)

plot(InfTemp)
lines(env.spline0, col="red")
lines(env.spline10, col="green")
lines(env.spline20, col="blue")
lines(env.spline30, col="yellow")
lines(env.spline40, col="magenta")
lines(env.spline50, col="cyan")

flux.spline0  <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=0)
flux.spline10 <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=10)
flux.spline20 <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=20)
flux.spline30 <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=30)
flux.spline40 <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=40)
flux.spline50 <- sm.spline(x=flux.data$Age, y=flux.data$A_Flux, df=50)

plot(flux.data)
lines(flux.spline0, col="red", lwd=4)
lines(flux.spline10, col="green", lwd=4)
lines(flux.spline20, col="blue", lwd=4)
lines(flux.spline30, col="yellow", lwd=4)
lines(flux.spline40, col="magenta", lwd=4)
lines(flux.spline50, col="cyan", lwd=4)


# Fit and plot Disparity Through Time to our data
dmin <- dtt(trim.tree, species.means, nsim=10000, plot=T)
dmean <- dtt(trim.tree, species.means, nsim=10000, plot=T)
dmax <- dtt(trim.tree, species.means, nsim=10000, plot=T)


# We can simulate data under the preferred model (GRASSexp or GRASSlin)
simGRASSexp <- sim_t_env(trim.tree, GRASSexp$best.result$param, model="EnvExp", env_data=testo,
          root.value=GRASSexp$best.result$root, plot=T)
simGRASSlin <- sim_t_env(trim.tree, GRASSlin$best.result$param, model="EnvLin", env_data=testo,
                         root.value=GRASSlin$best.result$root, plot=T)
#phenogram(trim.tree, macro.means)
#phenogram(trim.tree, simGRASSexp)
dtt(trim.tree, simGRASSexp, nsim=10000, plot=T) # looks similar for the youngest tree, very different for the oldest
dtt(trim.tree, simGRASSlin, nsim=100, plot=T) # looks similar for the youngest tree, very different for the oldest


best.shifttime <- function (phy, trait.data, time.data=InfTemp, degrees=c(0,10,20,30,40,50), model="EnvExp", cores=6) {
  res.list <- mclapply(1:length(degrees), function(x) {
    fit_t_env(phy, trait.data, env_data=time.data, df=degrees[x], scale=F, plot=T, model=model)}, mc.cores = cores)
  for(i in 1:length(res.list)){res.list[[i]]$df <- degrees[i]}
  res.values <- unlist(lapply(res.list, function(x) x$aicc)) # make a vector of the values, so we can get the index number of the best
  best.res <- res.list[[which.min(res.values)]]
  
  plot(best.res, main=paste(model, "; AICc = ", round(best.res$aicc,2)), sub=paste("sigma = ",round(best.res$param[1],2), " beta = ",round(best.res$param[2],2), " df = ", best.res$df), col="red")
  
  return(list(all.results=res.list, best.result=best.res, best.df=best.res$df))
}


# Now for the model adequacy bit:
simGRASS <- sim_t_env(trim.tree, param=GRASSexp$best.result$param, model="EnvExp", root.value=GRASSexp$best.result$root, env_data=testo)
      # phenogram(trim.tree, simGRASS)
arbutus(trim.tree, simGRASSexp)

simBM <- fastBM(trim.tree, sig2=BM_res$opt$sigsq)
arbutus(trim.tree, simBM)







temp.slice <- max(nodeHeights(trim.tree)) - 12.5
split.tree <- make.era.map(trim.tree, c(0,temp.slice))
OUBM_res <- mvSHIFT(split.tree, macro.means, model="OUBM", method="sparse", diagnostic=F, echo=F)


mclapply(1:length(degrees), function(x) {
  fit_t_env(trim.tree, macro.means, env_data=time.data, df=degrees[x], scale=F, plot=T, model=model)}, mc.cores = cores)

fit_t_env(trim.tree, macro.means, env_data=InfTemp, df=degrees[1], scale=F, plot=T, model="EnvExp")

arbutus()




# Plot the trees against Geological Time Scales
library(strap)
trim.tree$root.time<-max(nodeHeights(trim.tree))
geoscalePhylo(tree=ladderize(trim.tree, right=FALSE), label.offset=0.2, cex.age=0.6, cex.ts=0.8, cex.tip=0.8)

# Set up the trees we'll use
min.macro <- drop.tip(min.tree, setdiff(min.tree$tip.label, overlaps)); max(nodeHeights(min.macro))
mean.macro <- drop.tip(mean.tree, setdiff(mean.tree$tip.label, overlaps)); max(nodeHeights(mean.macro))
max.macro <- drop.tip(max.tree, setdiff(max.tree$tip.label, overlaps)); max(nodeHeights(max.macro))

#par(mfrow=c(1,3)); plot(min.macro, cex=0.3); plot(mean.macro, cex=0.3); plot(max.macro, cex=0.3)

min.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MinAges_Constrained_PF/Min_20.tre")
min.macros <- lapply(min.trees, drop.tip, tip=setdiff(min.trees[[1]]$tip.label, overlaps)); class(min.macros) <- "multiPhylo"; min.macros <- append(min.macros, min.macro)

mean.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MeanAges_Constrained_PF/Mean_20.tre")
mean.macros <- lapply(mean.trees, drop.tip, tip=setdiff(mean.trees[[1]]$tip.label, overlaps)); class(mean.macros) <- "multiPhylo"; mean.macros <- append(mean.macros, mean.macro)

max.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MaxAges_Constrained_PF/Max_20.tre")
max.macros <- lapply(max.trees, drop.tip, tip=setdiff(max.trees[[1]]$tip.label, overlaps)); class(max.macros) <- "multiPhylo"; max.macros <- append(max.macros, max.macro)

rescaled.trees <- append(min.macros, append(mean.macros, max.macros))

tree.span <- read.tree("/Users/Ian/Desktop/Macropod_Dating/Operators/Tree_Span.tre")
tree.span <- lapply(tree.span, drop.tip, tip=setdiff(tree.span[[79]]$tip.label, overlaps)); class(tree.span) <- "multiPhylo";
write.tree(tree.span, file="/Users/Ian/Desktop/Macropod_Dating/Operators/Tree_Span.tre")

# Function to rescale a number of trees 
rescale.series <- function(phy, min, max, interval){
  new.trees <- NULL
  ages <- seq(min, max, by = interval)
  for (j in 1:length(ages)){
    new.trees[[j]] <- geiger::rescale(phy, model="depth", ages[j])
  }
  class(new.trees) <- "multiPhylo"
  return(new.trees)
}

all.aics <- NULL; all.results <- NULL
for (k in 1:length(tree.span)){
  beginning <- Sys.time()
  int.results <- NULL
  
  # Fit a number of models to the data (ENV, GRASS, BM, EB, Trend, Drift)
  ENVexp <- best.smoothing(tree.span[[k]], macro.means, time.data=InfTemp, degrees=c(10,20,30,40,50), model="EnvExp", cores=5); int.results[["ENVexp"]] <- ENVexp$best.result; 
  ENVlin <- best.smoothing(tree.span[[k]], macro.means, time.data=InfTemp, degrees=c(10,20,30,40,50), model="EnvLin", cores=5); int.results[["ENVlin"]] <- ENVlin$best.result; 
  
  GRASSexp <- best.smoothing(tree.span[[k]], macro.means, time.data=testo, degrees=c(30,50), model="EnvExp", cores=3); int.results[["GRASSexp"]] <- GRASSexp$best.result; 
  #GRASSlin <- best.smoothing(trim.tree, macro.means, time.data=testo, degrees=c(30,40,50), model="EnvLin", cores=3)
  
  FLUXexp <- best.smoothing(tree.span[[k]], macro.means, time.data=flux.data, degrees=c(10,20,30,40,50), model="EnvExp", cores=5); int.results[["FLUXexp"]] <- FLUXexp$best.result; 
  FLUXlin <- best.smoothing(tree.span[[k]], macro.means, time.data=flux.data, degrees=c(10,20,30,40,50), model="EnvLin", cores=5); int.results[["FLUXlin"]] <- FLUXlin$best.result; 
  
  BM_res    <- fitContinuous(tree.span[[k]], species.means, model="BM"); int.results[["BM"]] <- BM_res
  trend_res <- fitContinuous(tree.span[[k]], species.means, model="trend"); int.results[["Trend"]] <- trend_res
  EB_res    <- fitContinuous(tree.span[[k]], species.means, model="EB"); int.results[["EB"]] <- EB_res
  
  curr_tree_FIT <- c(ENVexp$best.result$aicc, ENVlin$best.result$aicc, 
                    GRASSexp$best.result$aicc, FLUXexp$best.result$aicc, FLUXlin$best.result$aicc, #GRASSlin$best.result$aicc,
                    BM_res$opt$aicc, trend_res$opt$aicc, EB_res$opt$aicc); names(curr_tree_FIT) <- c("ENVexp", "ENVlin", "GRASSexp", "FLUXexp", "FLUXlin", "BM", "Trend", "EB")
  curr.aic <- as.data.frame(as.vector(aic.w(curr_tree_FIT))); 
        curr.aic$model <- names(curr_tree_FIT); colnames(curr.aic) <- c("aiccw", "model"); 
            curr.aic$age <- round(max(nodeHeights(tree.span[[k]])), 3)
                curr.aic$tree <- k
  all.results[[k]] <- int.results
  
  curr.aic$model.type <- c("ENV", "ENV", "GRASS", "FLUX", "FLUX", "BM", "Trend", "EB")
  
  # ENVexp <- ENVexp$best.result;  ENVlin <- ENVlin$best.result; GRASSexp <- GRASSexp$best.result; FLUXexp <- FLUXexp$best.result;  FLUXlin <- FLUXlin$best.result
  # curr.aic[which(curr.aic$aiccw == max(curr.aic$aiccw)),"model"]
  
  all.aics <- rbind.data.frame(all.aics, curr.aic); 
  
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time :", duration))
}
saveRDS(all.aics, file="/Users/Ian/Desktop/Macropod_Dating/Model_Fitting_AICCs.RDS")
all.results


#all.aic$age <- factor(all.aic$age, levels=c("Min_Ages", "Mean_Ages", "Max_Ages"))
all.aics$model.type <- factor(all.aics$model.type, levels=c("ENV", "FLUX", "GRASS", "BM", "EB", "Trend"))

test.aics$age <- as.factor(all.aics$age)

(ggplot(test.aics)
  + geom_bar(aes(y=aiccw, x=age, fill=model.type), stat="identity")
  + theme(axis.text.x=element_text(angle=90, hjust=1), panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual( values=wes_palette("Zissou1", 6, "continuous")))


max(nodeHeights(rescaled.trees[[1]])) - findMRCA(rescaled.trees[[1]], c("Dorcopsis_veterum", "Macropus_irma"), type="height")

lapply(rescaled.trees[[]], function(x) (max(nodeHeights(rescaled.trees[[x]])) - findMRCA(rescaled.trees[[x]], c("Dorcopsis_veterum", "Macropus_irma"), type="height")))

min.macros <- lapply(min.trees, drop.tip, tip=setdiff(min.trees[[1]]$tip.label, overlaps))

lapply(rescaled.trees, findMRCA, tips=c("Dorcopsis_veterum", "Macropus_irma"), type="height")



#############################################################################################
## Exercise to compare node ages across the different tree dating schemes
#############################################################################################

# Quick function to get the DEPTH of a node (from present), instead of the HEIGHT (from root)
MRCA.depth <- function(phy){max(nodeHeights(phy)) - findMRCA(phy, tips=c("Setonix_brachyurus", "Wallabia_bicolor"), type="height")}

#macros <- filter(macro.data, Higher_tax == "Macropodinae" | Higher_tax == "Sthenurinae" | Higher_tax == "Lagostrophinae" | Higher_tax == "Hypsiprymnodontidae")
#overlaps <- intersect(min.tree$tip.label, unique(macros$Taxon))

# Prepare trees for comparison of ages
max.age.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MaxAges_Constrained_PF/Newick_MAX_ages.tree"); max.ages <- max.age.trees[436:636]
#max.ages <- lapply(max.age.trees, drop.tip, tip=setdiff(max.age.trees[[1]]$tip.label, overlaps)); class(max.ages) <- "multiPhylo"; 
age.max <- as.data.frame(unlist(lapply(max.ages, MRCA.depth))); age.max$tree <- "strat.max"; colnames(age.max) <- c("age", "tree")

mean.age.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MeanAges_Constrained_PF/Newick_MEAN_species.trees"); mean.ages <- mean.age.trees[435:635]
#mean.ages <- lapply(mean.age.trees, drop.tip, tip=setdiff(mean.age.trees[[1]]$tip.label, overlaps)); class(mean.ages) <- "multiPhylo"
age.mean <- as.data.frame(unlist(lapply(mean.ages, MRCA.depth))); age.mean$tree <- "strat.mean"; colnames(age.mean) <- c("age", "tree")

min.age.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_MinAges_Constrained_PF/Newick_MIN_species.trees"); min.ages <- min.age.trees[433:633]
#min.ages <- lapply(min.age.trees, drop.tip, tip=setdiff(min.age.trees[[1]]$tip.label, overlaps)); class(min.ages) <- "multiPhylo"
age.min <- as.data.frame(unlist(lapply(min.ages, MRCA.depth))); age.min$tree <- "strat.min"; colnames(age.min) <- c("age", "tree")

age.cp.min  <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPMinAges_Constrained_PF/Macro_CP_MinAges_NEWICK.trees"); age.cp.min <- age.cp.min[(length(age.cp.min)-200):length(age.cp.min)]
cp.min <- as.data.frame(unlist(lapply(age.cp.min, MRCA.depth))); cp.min$tree <- "cp.min"; colnames(cp.min) <- c("age", "tree")

age.cp.mean <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPEstAges_Constrained_PF/Macro_CP_EstAges_NEWICK.trees"); age.cp.mean <- age.cp.mean[(length(age.cp.mean)-200):length(age.cp.mean)]
cp.est <- as.data.frame(unlist(lapply(age.cp.mean, MRCA.depth))); cp.est$tree <- "cp.est"; colnames(cp.est) <- c("age", "tree")

age.cp.max <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPMaxAges_Constrained_PF/Macro_CP_MaxAges_NEWICK.trees"); age.cp.max <- age.cp.max[(length(age.cp.max)-200):length(age.cp.max)]
cp.max <- as.data.frame(unlist(lapply(age.cp.max, MRCA.depth))); cp.max$tree <- "cp.max"; colnames(cp.max) <- c("age", "tree")

age.range.strict <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_AgesRanges_LinksCorrected_STRICT/Macro_AgesRanges_STRICT_NEWICK.trees"); age.range.strict <- age.range.strict[(length(age.range.strict)-200):length(age.range.strict)]
range.strict <- as.data.frame(unlist(lapply(age.range.strict, MRCA.depth))); range.strict$tree <- "range.strict"; colnames(range.strict) <- c("age", "tree")

#age.range.relax <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_AgesRanges_LinksCorrected_UCLN/Macro_AgesRanges_UCLN_NEWICK.trees"); age.range.relax <- age.range.relax[(length(age.range.relax)-200):length(age.range.relax)]
#range.relax <- as.data.frame(unlist(lapply(age.range.relax, MRCA.depth))); range.relax$tree <- "range.ucln"; colnames(range.relax) <- c("age", "tree")

age.all <- rbind(age.min, cp.min, cp.est, cp.max, age.max, range.strict)
age.all$tree <- factor(age.all$tree, levels=c("strat.max", "cp.max", "range.strict", "cp.est", "cp.min", "strat.min"))


macro_lago.plot <- (ggplot(age.all, aes(x=age, fill=tree))
  + geom_density(alpha=0.75, adjust=1.5)
  + theme(axis.text.x=element_text(angle=0, hjust=1), panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Zissou1", type="continuous", 6))
  + scale_x_reverse(lim=c(32,10)))
  #+ scale_fill_manual(values=c("#F21A00", "#EBCC2A", "#3B9AB2")))

macropodini.plot; macropodinae.plot; macro_sthen.plot; macro_lago.plot; root.plot

grid.arrange(macropodini.plot,
             macropodinae.plot,
             macro_sthen.plot, 
             #sthenurinae.plot,
             macro_lago.plot,
             root.plot,
             nrow=5)

cp.mean$root.time<-max(nodeHeights(cp.mean))
geoscalePhylo(tree=ladderize(cp.mean, right=FALSE), label.offset=0.2, cex.age=0.6, cex.ts=0.8, cex.tip=0.8)


#### To get the Bayes Factors for Fossil Taxa that are putative Sampled Ancestors
####################################################################################

prior.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Priors_Only/Prior_Only_NEWICK.trees"); prior.trees <- prior.trees[(length(prior.trees)-200):length(prior.trees)]
post.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_CPEstAges_Constrained_PF/Macro_CP_EstAges_NEWICK.trees"); post.trees <- post.trees[(length(post.trees)-200):length(post.trees)]
post.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_AgesRanges_LinksCorrected_STRICT/Macro_AgesRanges_STRICT_NEWICK.trees"); post.trees <- post.trees[(length(post.trees)-200):length(post.trees)]
post.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_AgesRanges_LinksCorrected_UCLN/Macro_AgesRanges_UCLN_NEWICK.trees"); post.trees <- post.trees[(length(post.trees)-200):length(post.trees)]
post.trees <- read.nexus("/Users/Ian/Desktop/Macropod_Dating/Operators/Macro_AgesRanges_Constrained_PF/Macro_AgesRanges_Constrained_NEWICK.trees"); post.trees <- post.trees[(length(post.trees)-200):length(post.trees)]


## Choose which tips you want information for:
fossil_taxa <- c("Baringa_nelsonensis", "Bohra_illuminata", "Congruus_congruus",
                 "Dorcopsoides_fossilis", "Ganguroo_bilamina","Hadronomas_puckridgi",
                 "Kurrabi_mahoneyi","Macropus_pavana","Ngamaroo_archeri",
                 "Bulungamaya_delicata","Prionotemnus_palankarinnicus",
                 "Procoptodon_goliah","Protemnodon_anak","Simosthenurus_occidentalis",
                 "Sthenurus_andersoni","Troposodon_minor","Wanburoo_hilarus")

## Get the node numbers of the tips
nodes <- sapply(fossil_taxa,function(x,y) which(y==x),y=tree$tip.label)

## then get the edge lengths for those nodes
edge.lengths <- setNames(tree$edge.length[sapply(nodes,
                                               function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))

## The faster way is to make a function to do this:
get_terminal_branchlengths <- function(phy, tipnames){
  ## Get the node numbers of the tips
  nodes <- sapply(tipnames,function(x,y) which(y==x),y=phy$tip.label)
  ## Then get the edge lengths for those nodes
  edge.lengths <- setNames(phy$edge.length[sapply(nodes,
                                                   function(x,y) which(y==x),y=phy$edge[,2])],names(nodes))
  return(edge.lengths)
}
#get_terminal_branchlengths(test.tree, tipnames=c("Macropus_pavana", "Protemnodon_anak"))

## Now that we've got the tips and branch lengths, we can compare the posterior to the prior
BFSA <- function(prior.phy, posterior.phy, tips){
  post <- lapply(posterior.phy, get_terminal_branchlengths, tipnames=tips); names(post) <- NULL; post <- unlist(post)
  prior <- lapply(prior.phy,    get_terminal_branchlengths, tipnames=tips); names(prior)<- NULL; prior <- unlist(prior)
  
  BFs <- NULL
  for (j in 1:length(tips)){
    curr.tip <- subset(post, names(post)==tips[j]); 
    probSA <- sum(curr.tip<=0); probTIP <- length(curr.tip)-probSA;
    
    curr.tip <- subset(prior, names(prior)==tips[j]);
    priorSA <- sum(curr.tip<=0); priorTIP <- length(curr.tip)-priorSA;
    
    curr.BF <- log((probSA * priorTIP) / (probTIP * priorSA))
    if(is.na(curr.BF)){curr.BF <- 0}
    
    #curr.BF <- log(probSA/(length(curr.tip)-probSA))
    names(curr.BF) <- tips[j]; curr.BF <- round(curr.BF, 2)
    BFs <- append(BFs, curr.BF)
  }
  return(BFs)
}

#macro_BFs <- BFSA(age.cp.mean, tips=fossil_taxa)
macro_BFs <- BFSA(prior.trees, post.trees, tips=fossil_taxa)

macro_BFs <- as.data.frame(macro_BFs) # make the vector a data frame
macro_BFs[which(macro_BFs$macro_BFs > 5),] <- 5 # change any really big (INF) numbers to 5
macro_BFs[which(macro_BFs$macro_BFs < -5),] <- -5 # change any really small (-INF) numbers to -5

macro_BFs$taxa <- rownames(macro_BFs); # create column with the taxon names
macro_BFs <- macro_BFs[order(macro_BFs$macro_BFs),] # reorder by BF values
# set colors for plotting
macro_BFs$color <- "black";
macro_BFs[which(macro_BFs$macro_BFs > 1),]$color <- "#3d98d3"
macro_BFs[which(macro_BFs$macro_BFs < -1),]$color <- "#FF7175"

macro_BFs$taxa <- factor(macro_BFs$taxa, levels=c(macro_BFs$taxa)) # set factors for plotting (NOT NECESSARY)


BF.strict <- ggplot(macro_BFs, aes(x=taxa, y=macro_BFs, label=macro_BFs)) + 
  geom_ribbon(ymin=-1, ymax=+1) +
  geom_point(stat='identity', size=8, color=macro_BFs$color)  +
  # geom_segment(aes(y = 0, 
  #                  x = taxa, 
  #                  yend = macro_BFs, 
  #                  xend = taxa), 
  #              color = macro_BFs$color) +
  geom_text(color="white", size=2) +
  #labs(title="Bayes Factor Support", subtitle="for Fossil Taxa as Sampled Ancestors") + 
  #ylim(-5, 5) +
  scale_y_continuous(name="log Bayes Factors", limits=c(-5,5), breaks=c(-5:5)) +
  #theme(panel.background=element_blank()) +
  geom_hline(yintercept=-1) +
  geom_hline(yintercept=1) +
  xlab("Fossil Taxa") +
  #ylab("log Bayes Factors") +
  theme_classic() +
  coord_flip()

#plot_grid(BF.ucln, BF.strict, BF.CPest)
plot_grid(BF.strict, BF.CPest)


#### Plot the C4 and Flux data on a geological time scale
library(deeptime)

pp <- ggplot(enviro.data, aes(Age)) +
  geom_ribbon(aes(ymin = C4_recon_lower, ymax = C4_recon_upper), fill = "pink") +
  geom_line(aes(y = C4_recon_mean), color="red") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0,80), expand = FALSE) 

qq <- gggeo_scale(p, dat="epochs")


rr <- ggplot(flux.data, aes(Age)) +
  geom_ribbon(aes(ymin = A_Flux-35, ymax = A_Flux+35), fill = "light blue") +
  geom_line(aes(y = A_Flux), color="blue") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 13), ylim = c(0,150), expand = FALSE) 

ss <- gggeo_scale(rr, dat="epochs")

grid.arrange(qq, ss, nrow=1)



## Let's plot the estimated ages of each fossil taxon
####################################################################################

# Create a function to pull the ages of each fossil taxon estimated
get.fossil.ages <- function(fossil.tips, trees){
  tree.tables <- lapply(trees, print.tree)
  fossil.tables <- lapply(1:length(tree.tables), function(x) {
    subset(tree.tables[[x]], tree.tables[[x]]$label %in% fossil.tips)
  })
  fossil.ages <- lapply(1:length(fossil.tables), function(x) {
    select(fossil.tables[[x]], label, time_bp)
  })
  final <- bind_rows(fossil.ages)
}

my.test <- get.fossil.ages(fossil.tips = fossil_taxa, trees = post.trees)
my.test$label <- factor(my.test$label, levels=c("Bulungamaya_delicata",
                                                "Protemnodon_anak",
                                                "Sthenurus_andersoni",
                                                "Troposodon_minor",
                                                "Baringa_nelsonensis",
                                                "Prionotemnus_palankarinnicus",
                                                "Congruus_congruus",
                                                "Bohra_illuminata",
                                                "Kurrabi_mahoneyi",
                                                "Ngamaroo_archeri",
                                                "Dorcopsoides_fossilis",
                                                "Wanburoo_hilarus",
                                                "Ganguroo_bilamina",
                                                "Hadronomas_puckridgi",
                                                "Simosthenurus_occidentalis",
                                                "Macropus_pavana",
                                                "Procoptodon_goliah"))

(ggplot(my.test, aes(x=time_bp, y=label, fill=..x..)) 
  + scale_fill_gradientn(colours=wes_palette("Zissou1"))
  + geom_density_ridges_gradient(scale=2)
  + scale_x_reverse()
  + theme_classic())

## We can also compare the mean fossil age estimates from different analyses  
cpest.test <- my.test
range.test <- my.test

range.means <- range.test %>%
  group_by(label) %>%
  summarise_at(vars(time_bp), mean)
range.means <- as.data.frame(range.means)

cpest.means <- cpest.test %>%
  group_by(label) %>%
  summarise_at(vars(time_bp), mean)
cpest.means <- as.data.frame(cpest.means)

range.means <- range.means[order(range.means$label),]; names(range.means) <- c("label", "time_bp_range")
cpest.means <- cpest.means[order(cpest.means$label),]; names(cpest.means) <- c("label", "time_bp_cpest")

both.means <- cbind(cpest.means, range.means$time_bp_range); names(both.means) <- c("label", "time_bp_cpest", "time_bp_range")
both.means[,2] <- round(both.means[,2],2); both.means[,3] <- round(both.means[,3],2)


(ggplot(both.means, aes(x=time_bp_cpest, xend=time_bp_range, y=label, group=label))
        + geom_dumbbell(color="light green", colour_x="light green", colour_xend="darkGreen", size_x=6, size_xend=6, lwd=1)
        + theme_classic()
        + scale_x_reverse()
        + geom_text(color="light green", size=2, hjust=-2,
                    aes(x=time_bp_cpest, label=time_bp_cpest))
        + geom_text(aes(x=time_bp_range, label=time_bp_range), 
                    color="darkGreen", size=2, hjust=2))


