### BMC Bioinformatics 
### Designs for the simultaneous inference of concentration-response curves
### Authors: Leonie Schürmeyer, Kirsten Schorning and Jörg Rahnenführer 
### 2023

# This R-Code provides the main R-Code and ideas used in the analysis of the paper
# "Designs for the simultaneous inference of concentration-response curves". 
# The central functions for running the code can be found in 
# "BMCBioinformatics_SchuermeyerSchorningRahnenfuehrer_CentralFunctions.R"

# Corresponding datasets or intermediate results can also be found on the GitHub page 
# https://github.com/schuermeyer/Designs-for-the-simultaneous-inference-of-concentration-response-curves


####################################################################################################

# Load packages 
library(DoseFinding)
library(MASS)
library(Metrics)
library(sfsmisc)

####################################################################################################

# Preprocessing of the data

# Load data
load("VPA_Dataset.RData")
load("VPADataset_matrixVersion.RData")

doses.single <- c(0, 25, 150, 350, 550, 800, 1000)
doses <- c(rep(0, 6), rep(c(25, 150, 350, 450, 550, 800, 1000), each=3))


# First step: Find genes with biological activity using MCPMod-procedure to extract genes with 
# convenient sigmoidal fit


# Define candidate models for MCPMod-procedure, for increasing (up) and decreasing (down) fit
candMods.up <- Mods(sigEmax = c(450, 5.117523),
                    doses=doses.single,
                    maxEff=1)

candMods.down <- Mods(sigEmax=c(450, 5.117523),
                      doses=doses.single,
                      maxEff=-1)


# Extract p-values for all genes with MCPMod-procedure while using sigEmax model both for 
# up and down directed models (ATTENTION: long runtime)

#### DO NOT RUN
# pvalues.up <- apply(ev.VPA.without650, 1, function(x){
#   dataframe <- data.frame(dose=doses, resp=x)
#   pval.up <- attr(suppressWarnings(MCTtest(dose=dose, resp=resp, 
#                                            data=dataframe, models=candMods.up))$tStat, which = "pVal")
#   pval.up
# })
# 
# pvalues.down <- apply(ev.VPA.without650, 1, function(x){
#   dataframe <- data.frame(dose = doses, resp=x)
#   pval.down <- attr(suppressWarnings(MCTtest(dose=dose, resp=resp, 
#                                              data=dataframe, models=candMods.down))$tStat, which = "pVal")
#   pval.down
# })
# 
# 
# all.pvalues <- cbind(pvalues.down, pvalues.up)

# Because of the long runtime, results can also be found in: 
load("PoCPvalues.VPA.RData")

# Overview, which genes show biologic activity (extract those with resulting pvalue < 0.01)
pvalues <- as.data.frame(all.pvalues)
pval <- 0.01
bactiv <- which(pvalues$pvalues.down < pval | pvalues$pvalues.up < pval) # 20791 genes show biologic activity

# Extract gene names
Genenames_all <- rownames(pvalues)
Genenames_bactiv <- Genenames_all[bactiv]

# Only extract genes with biologic activity
ind <- which(all_data$gene %in% Genenames_bactiv)
bdata <- all_data[ind,]



#### DO NOT RUN
# # Fit sigmoid Emax model for every gene separately
# Models_sigEmax <- vector(mode="list", length = length(Genenames_bactiv))
# 
# # Set maximal Dose
# mD <- 1000
# 
# # Fix boundary for sigmoid Emax to c(0.05,100) for h and c(0,1500) for EC50
# for(i in 1:length(Genenames_bactiv)){
#   Models_sigEmax[[i]] <- fitMod(dose,resp, data=bdata[bdata$gene==Genenames_bactiv[i],], 
#                                  model="sigEmax", 
#                                  bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 100), 2))
# }
# # Extract all parameters
# Par_bactiv <- data.frame(e0=numeric(length(Genenames_bactiv)),
#                          eMax=numeric(length(Genenames_bactiv)),
#                          EC50=numeric(length(Genenames_bactiv)),
#                          h=numeric(length(Genenames_bactiv)))
# 
# for(i in 1:length(Genenames_bactiv)){
#   Par_bactiv[i,1] <- Models_sigEmax[[i]]$coefs[1]
#   Par_bactiv[i,2] <- Models_sigEmax[[i]]$coefs[2]
#   Par_bactiv[i,3] <- Models_sigEmax[[i]]$coefs[3]
#   Par_bactiv[i,4] <- Models_sigEmax[[i]]$coefs[4]
# }
# 
# # Now fit each model with a higher steepness h>10 again with constraint to h<=10,
# # because higher steepness parameters lead to numerical problems in further analysis. Furthermore 
# # a steepness over 10 leads to nearly no change in the model compared to a steepness of 10.
# Genes_over10 <- Genenames_bactiv[Par_bactiv$h>10]
# 
# for(i in which(Par_bactiv$h>10)){
#   Models_sigEmax[[i]] <- fitMod(dose,resp, data=bdata[bdata$gene==Genenames_bactiv[i],], 
#                                  model="sigEmax", 
#                                  bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
# }
# 
# # Extract parameters
# Par_bactiv <- data.frame(e0=numeric(length(Genenames_bactiv)),
#                          eMax=numeric(length(Genenames_bactiv)),
#                          EC50=numeric(length(Genenames_bactiv)),
#                          h=numeric(length(Genenames_bactiv)))
# 
# for(i in 1:length(Genenames_bactiv)){
#   Par_bactiv[i,1] <- Models_sigEmax[[i]]$coefs[1]
#   Par_bactiv[i,2] <- Models_sigEmax[[i]]$coefs[2]
#   Par_bactiv[i,3] <- Models_sigEmax[[i]]$coefs[3]
#   Par_bactiv[i,4] <- Models_sigEmax[[i]]$coefs[4]
# }
# 
# # Also exclude genes with an EC50 not available inside the design space of [0,1000]
# ind_ec1000 <- which(Par_bactiv$EC50>999)
# relGene <- Genenames_bactiv[-ind_ec1000]
# 
# # Only extract relevant parameters and models and store them
# relPar <- Par_bactiv[-ind_ec1000,]
# ind <- which(all_data$gene %in% relGene)
# relData <- all_data[ind,]
# Models_sigEmax_rel <- Models_sigEmax[-ind_ec1000]

# Parameters can be found in 
load("relPar.RData")
load("Models.RData")

####################################################################################################


# Calculate locally D-optimal designs 

# Calculate locally D-optimal designs for each gene using Particle Swarm Optimization (PSO)
# Therefore, load own functions psoOptDesign_fixed used for locally D-optimal designs with fixed weights 
# and boundary support points and other functions like "mygrad", "eq_func" 
# from "BMCBioinformatics_SchuermeyerSchorningRahnenfuehrer_CentralFunctions.R"


# Set design space
Lb0 = 0
Ub0 = 1000

# Initialize object to store variables 
optDesigns <- vector(mode="list", length=length(relPar$e0))


#### DO NOT RUN:
# # Use PSO to calculate locally optimal designs for every gene (ATTENTION: high runtime)
# for(i in 1:length(relPar$e0)){
#   # fix parameter values in each step
#   mytheta = as.numeric(relPar[i,])
#   
#   # Run PSO with fixed weights and include boundary support points 
#   opt_des = psoOptDesign_fixed(crit=Dcrit, control = list(numIt = 100, numPart = 300,
#                                                     setProgressBar= FALSE),
#                          nPoints=4, wFixed=rep(1/4, 4), Lb= Lb0, Ub =Ub0, 
#                          gradient= mygrad, theta= mytheta, boundarySup = TRUE)
#   
#   optDesigns[[i]] <- opt_des
# }

# Sometimes PSO needs to be run multiple times or with higher iterations and particles to
# find the optimal design


# Load stored optimal designs and check if they are all locally D-optimal
load("optDesigns.RData")

# Use Equivalence theorem to check D-optimality 

# Initialize sequence over design space 
xseq =seq(Lb0, Ub0, length= 1001)

# Attention long runtime:
yseqDesigns <- vector(mode="list", length=length(optDesigns))
for(i in 1:length(optDesigns)){
  yseq1 <- try(vapply(xseq, FUN= eq_func, FUN.VALUE=1, design=optDesigns[[i]],
                     gradient=mygrad,  theta=as.numeric(relPar[i,])), silent=TRUE)

  if(inherits(yseq1, "try-error")) { yseqDesigns[[i]] <- NA }

  else { yseqDesigns[[i]] <- yseq1 }
}

EquValue <- rep(NA,length(optDesigns))
for(i in 1:length(optDesigns)){
  EquValue[i] <- sum(which(yseqDesigns[[i]]>0.01))
} # No values over 0.01 => all optimal, 
# except for 85 genes, they lead to singular matrices, remove those and reduce dataset to 15233 genes

# Check NAs (singular info matrix)
checkNA <- rep(0,length(optDesigns))
for(i in 1:length(optDesigns)){
  checkNA[i] <- sum(is.na(yseqDesigns[[i]]))
}
withNA <- which(checkNA>0)
length(withNA)#85 singular matrices, numerical problems 

# Only use the optimal designs without singular info matrix  
length(relPar$e0) - length(withNA)#15233
opt_new <- 1:length(optDesigns)
opt_new <- opt_new[-withNA]

relPar_opt <- relPar[opt_new,]
Models_opt <- Models[opt_new]

####################################################################################################

# Develop kmeans design

# Initialize list 
sortsup <- vector(mode="list", length=length(optDesigns))

# Only look at genes in opt_new which have a corresponding locally D-optimal design 
for(i in opt_new){
  sortsup[[i]] <- sort(optDesigns[[i]]$supPoints)
}

# Remove 1st and 4th support point (because C0={{0},{1000}})
midsup <- vector(mode="list", length=length(sortsup))
for(i in 1:length(sortsup)){
  midsup[[i]] <- sortsup[[i]][2:3]
}
midsup <- unlist(midsup)
midsup <- unname(midsup)
midsuppoints <- round(midsup,digits=0)


# More robust approach by developing kmeans design 500 times and get the mean value per support point

kmeans500 <- data.frame(T1=numeric(500),T2=numeric(500),T3=numeric(500),
                        T4=numeric(500),T5=numeric(500),
                        T6=numeric(500),T7=numeric(500), 
                        w1=numeric(500),w2=numeric(500),w3=numeric(500),
                        w4=numeric(500),w5=numeric(500),
                        w6=numeric(500),w7=numeric(500))

set.seed(157)
for(i in 1:500){
  k <- kmeans(midsuppoints, centers=7)
  kmeans500[i,1:7] <- as.numeric(sort(k$centers))
  kmeans500[i,8:14] <- as.numeric((table(k$cluster)[order(k$centers)]/sum(table(k$cluster)))*7/9)
}

kmeansDesign500 <- vector(mode="list", length=1)
kmeansDesign500$weights <- rep(1/9,9)
kmeansDesign500$supPoints <-  round(c(0, as.numeric(colMeans(kmeans500)[1:7]), 1000))


####################################################################################################

# Develop simultaneous inference design 

# First introduce a 5x5 grid over parameter space with steepness h and EC50

# Initialize dataframe for grid 
# Only use parameters corresponding to genes with present locally D-optimal design 

dfarea <- data.frame(Gene=1:length(relPar_opt$h), h= relPar_opt$h, e=relPar_opt$EC50, 
                       area= numeric(length(relPar_opt$e0)))

# Create areas of grid classification 
dfarea$area[dfarea$h > 0 & dfarea$h <= 2 
                & dfarea$e > 0 & dfarea$e <= 200] <- 1

dfarea$area[dfarea$h > 0 & dfarea$h <= 2 
                & dfarea$e > 200 & dfarea$e <= 400] <- 2

dfarea$area[dfarea$h > 0 & dfarea$h <= 2 
                & dfarea$e > 400 & dfarea$e <= 600] <- 3

dfarea$area[dfarea$h > 0 & dfarea$h <= 2 
                & dfarea$e > 600 & dfarea$e <= 800] <- 4

dfarea$area[dfarea$h > 0 & dfarea$h <= 2 
                & dfarea$e > 800 & dfarea$e <= 1000] <- 5

dfarea$area[dfarea$h > 2 & dfarea$h <= 4 
                & dfarea$e > 0 & dfarea$e <= 200] <- 6

dfarea$area[dfarea$h > 2 & dfarea$h <= 4
                & dfarea$e > 200 & dfarea$e <= 400] <- 7

dfarea$area[dfarea$h > 2 & dfarea$h <= 4
                & dfarea$e > 400 & dfarea$e <= 600] <- 8

dfarea$area[dfarea$h > 2 & dfarea$h <= 4
                & dfarea$e > 600 & dfarea$e <= 800] <- 9

dfarea$area[dfarea$h > 2 & dfarea$h <= 4
                & dfarea$e > 800 & dfarea$e <= 1000] <- 10

dfarea$area[dfarea$h > 4 & dfarea$h <= 6 
                & dfarea$e > 0 & dfarea$e <= 200] <- 11

dfarea$area[dfarea$h > 4 & dfarea$h <= 6 
                & dfarea$e > 200 & dfarea$e <= 400] <- 12

dfarea$area[dfarea$h > 4 & dfarea$h <= 6 
                & dfarea$e > 400 & dfarea$e <= 600] <- 13

dfarea$area[dfarea$h > 4 & dfarea$h <= 6 
                & dfarea$e > 600 & dfarea$e <= 800] <- 14

dfarea$area[dfarea$h > 4 & dfarea$h <= 6 
                & dfarea$e > 800 & dfarea$e <= 1000] <- 15

dfarea$area[dfarea$h > 6 & dfarea$h <= 8 
                & dfarea$e > 0 & dfarea$e <= 200] <- 16

dfarea$area[dfarea$h > 6 & dfarea$h <= 8 
                & dfarea$e > 200 & dfarea$e <= 400] <- 17

dfarea$area[dfarea$h > 6 & dfarea$h <= 8 
                & dfarea$e > 400 & dfarea$e <= 600] <- 18

dfarea$area[dfarea$h > 6 & dfarea$h <= 8 
                & dfarea$e > 600 & dfarea$e <= 800] <- 19

dfarea$area[dfarea$h > 6 & dfarea$h <= 8 
                & dfarea$e > 800 & dfarea$e <= 1000] <- 20

dfarea$area[dfarea$h > 8 & dfarea$h <= 10 
                & dfarea$e > 0 & dfarea$e <= 200] <- 21

dfarea$area[dfarea$h > 8 & dfarea$h <= 10
                & dfarea$e > 200 & dfarea$e <= 400] <- 22

dfarea$area[dfarea$h > 8 & dfarea$h <= 10 
                & dfarea$e > 400 & dfarea$e <= 600] <- 23

dfarea$area[dfarea$h > 8 & dfarea$h <= 10 
                & dfarea$e > 600 & dfarea$e <= 800] <- 24

dfarea$area[dfarea$h > 8 & dfarea$h <= 10 
                & dfarea$e > 800 & dfarea$e <= 1000] <- 25

# Create table of frequencies of all genes present in each area 
tab <- data.frame(table(dfarea$area))
weight <- sum(tab$Freq)
tab$weight <- tab$Freq/weight
t(tab)

# Now draw a representative for each area 
set.seed(197)
Reparea <- data.frame(area=1:25, weight= tab$weight, Rep=numeric(25))

for(i in 1:25){
  Reparea$Rep[i] <- sample(dfarea$Gene[dfarea$area==i],1)
}

Reparea$E0 <- numeric(25)
Reparea$EMax <- numeric(25)
Reparea$EC50 <- numeric(25)
Reparea$h <- numeric(25)

for(i in 1:25){
  Reparea[i,4:7] <- as.numeric(relPar_opt[Reparea$Rep[i],])
}

# Only use areas which are represented higher than 5%, reduces to 7 representatives 
Reparea7 <- Reparea[which(Reparea$weight>0.05),]
weight7 <- Reparea7$weight/sum(Reparea7$weight)
Prior7 <- Reparea7
theta7 <- Reparea7[,4:7]

# Calculate locally D-optimal designs for the 7 representatives
OD7 <- vector(mode="list", length=length(Prior7$E0))
for(i in 1:length(Prior7$E0)){
  # Fix representative parameter values for each area 
  mytheta= as.numeric(theta7[i,])
  
  # Use PSO to calculate each locally D-optimal design 
  opt_des = psoOptDesign_fixed(crit=Dcrit, control = list(numIt = 100, numPart = 300,
                                                    setProgressBar= FALSE),
                         nPoints=4, wFixed=rep(1/4, 4), Lb= Lb0, Ub =Ub0, 
                         gradient= mygrad, theta= mytheta, boundarySup = TRUE)
  
  OD7[[i]] <- opt_des
}

# Calculate D-optimality criterion value
Dcrit7 <- numeric(7)
for(i in 1:7){
  Dcrit7[i] <- Dcrit(w=rep(1/4, 4), x= as.numeric(OD7[[i]]$supPoints),
                      gradient= mygrad, theta= as.numeric(theta7[i,]))
}


# Calculate simultaneous inference design with PSO 
# Use following criterion only considering the 7 representatives 


# Set reciprocal condition number to e^-20 
critnumber <- 0.00000000000000000001

# Only calculate simultaneous inference design with 7 representative genes and corresponding weights 
# for the 7 areas

Sim_Crit7 <- function(w,x,gradient, pillarweight, thetapillar, Dcrit7, ...){
  # Initialize criterion value 
  Sim_Crit <- 0
  for(i in 1:7){
    # Use condition number to avoid singular matrices 
    rcond_var <- sapply(1:7, FUN = function(i){rcond(info(w,x,gradient, 
                                                          theta=as.numeric(thetapillar[i,])))})
    if(min(rcond_var > critnumber )){
      current <- pillarweight[i]*((det(info(w,x,gradient, 
                                            theta=as.numeric(thetapillar[i,]), ...)))^(1/4)/ Dcrit7[i])
      Sim_Crit <- Sim_Crit+current
    }
    else(Sim_Crit <- -Inf)
  }
  return(Sim_Crit)
}


# DO NOT RUN because of high runtime:
# Lb0 = 0
# Ub0 = 1000
# 
# starttime <- Sys.time()
# simDesign = psoOptDesign(crit=Sim_Crit7, control = list(numIt = 2000, numPart = 3000, 
#                                                                setProgressBar= TRUE),
#                               nPoints=9, Lb= Lb0, Ub =Ub0, gradient= mygrad,
#                               thetapillar=theta7, Dcrit7=Dcrit7,
#                               pillarweight= weight7)
# endtime <- Sys.time()
# time <- endtime-starttime
# Increase particles and iterations if optimal design cannot be found 


# Check results with equivalence sentence

simultaneousDesign <- vector(mode="list", length=2)
simultaneousDesign$weights <- c(17.2594,4.96305,11.541,11.9595,11.2259, 14.4871,2.68754, 6.42786,19.4487)/100
simultaneousDesign$supPoints <- c(0.0,145.286,279.822,345.414,456.865,575.21 ,656.139,781.354,1000.0)



# Calculate simultaneous efficiency for each area
simultaneouseff7 <- sapply(1:7, 
                        FUN = function(i){(det(info(w=simultaneousDesign$weights,
                                                    x=simultaneousDesign$supPoints,
                                                    gradient=mygrad,
                                                    theta=as.numeric(theta7[i,]))))^(1/4)/ Dcrit7[i]})
sum(simultaneouseff7 * weight7)# efficiency of 0.833

# Calculate criterion value for every x in design space [0,1000]
xseq =seq(Lb0, Ub0, length= 1001)
yseq = vapply(xseq, FUN= eq_funcsimultaneous7, FUN.VALUE=1, design=simultaneousDesign,
              gradient=mygrad, pillarweight=weight7, 
              simultaneouseff = simultaneouseff7,thetapillar = theta7)

# Plot the equivalence theorem to see if the design is simultaneous D-optimal 
plot(xseq, yseq, type ="l", main="",
     ylab=expression(s(x,xi[Theta[7]],pi)), xlab="concentration")
abline(h= 0, col="grey") #simultaneous D-optimal 


####################################################################################################

# Original design 
origDesign <- kmeansDesign500
origDesign$weights <- c(2/9,rep(1/9,7))
origDesign$supPoints <- c(0,25,150,350,450,550,800,1000)

# Equidistant design 
equiDesign <- kmeansDesign500
equiDesign$supPoints <- c(0,125,250,375,500,625,750,875,1000)

# Log-equidistant design 
logequiDesign <- kmeansDesign500
logequiDesign$supPoints <- round(c(0, lseq(from=1, to=1000, length=8)))



####################################################################################################

### D-efficiency comparison 

# Calculate D-criterion value of locally optimal designs for every gene
Dcrit_all <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all[i] <- Dcrit(w=rep(1/4, 4), x= as.numeric(optDesigns[[i]]$supPoints), 
                        gradient= mygrad, theta= as.numeric(relPar[i,]))
}


# Now for every gene with the different designs 

# Simultaneous design  
Dcrit_all_simultaneous <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all_simultaneous[i] <- Dcrit(w=as.numeric(simultaneousDesign$weights), 
                                     x= as.numeric(simultaneousDesign$supPoints), 
                              gradient= mygrad, theta= as.numeric(relPar[i,]))
}
which(is.na(Dcrit_all_simultaneous))#NA of genes with singular matrix like before 


# Kmeans design 
Dcrit_all_kmeans <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all_kmeans[i] <- Dcrit(w=as.numeric(kmeansDesign500$weights), 
                               x= as.numeric(kmeansDesign500$supPoints), 
                               gradient= mygrad, theta= as.numeric(relPar[i,]))
}
which(is.na(Dcrit_all_kmeans)) %in% which(is.na(Dcrit_all_simultaneous))#same as before


# Original design
Dcrit_all_orig <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all_orig[i] <- Dcrit(w=as.numeric(origDesign$weights), x= as.numeric(origDesign$supPoints), 
                             gradient= mygrad, theta= as.numeric(relPar[i,]))
}
which(is.na(Dcrit_all_orig)) %in% which(is.na(Dcrit_all_simultaneous))#same as before


# Equidistant design
Dcrit_all_equi <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all_equi[i] <- Dcrit(w=as.numeric(equiDesign$weights), x= as.numeric(equiDesign$supPoints), 
                             gradient= mygrad, theta= as.numeric(relPar[i,]))
}
which(is.na(Dcrit_all_equi)) %in% which(is.na(Dcrit_all_simultaneous))#same as before 

# Log-equidistant design
Dcrit_all_logequi <- rep(NA, length(optDesigns))
for(i in opt_new){
  Dcrit_all_logequi[i] <- Dcrit(w=as.numeric(logequiDesign$weights), 
                                x= as.numeric(logequiDesign$supPoints), 
                                gradient= mygrad, theta= as.numeric(relPar[i,]))
}
which(is.na(Dcrit_all_logequi)) %in% which(is.na(Dcrit_all_simultaneous))
# Nearly overall same NAs (except logequidistant design)

simultaneous_DEff <- Dcrit_all_simultaneous/Dcrit_all
kmeans_DEff <- Dcrit_all_kmeans/Dcrit_all
orig_DEff <- Dcrit_all_orig/Dcrit_all
equi_DEff <- Dcrit_all_equi/Dcrit_all
logequi_DEff <- Dcrit_all_logequi/Dcrit_all

length(which(is.na(simultaneous_DEff)))#85 for all except logequi

# Remove the genes which lead to singular matrices 85 and store all calculated efficiency values in 
# D_eff 
remove <- which(is.na(simultaneous_DEff))
D_Eff <- data.frame(Defficiency= c(rep("Simultaneous", length(optDesigns)-length(remove)),
                                        rep("kmeans", length(optDesigns)-length(remove)),
                                        rep("original", length(optDesigns)-length(remove)),
                                        rep("equidistant", length(optDesigns)-length(remove)),
                                        rep("logequidistant", length(optDesigns)-length(remove))),
                         Critvalue = c(simultaneous_DEff[-remove],
                                       kmeans_DEff[-remove],
                                       orig_DEff[-remove],
                                       equi_DEff[-remove],
                                       logequi_DEff[-remove]))



####################################################################################################

# Simulation study 
# Load parameters  
load("relPar.RData")
load("Models.RData")


# Setup simulation parameters for all designs 
simnum <- 500 # Number of simulation steps
space <- 1:1000 # Design space
new <- data.frame(dose = space, iwas=factor(1)) # Initialize sequence for prediction of new curves 
n_opt <- c(18,27,36,45,63,90) # Observation numbers 

# Run simulation study for all genes with a locally D-optimal design 
ind <- opt_new

# Calculate corresponding standard deviation sigma with range of each gene 
sigma <- numeric(length(Models))
for(i in ind){
  sigma[i] <- 0.2*abs(unname(Models[[i]]$coefs[2]))
}

# Calculate original curve of each gene corresponding to the fitted model  
Origin <- vector(mode="list", length=length(Models))
for(i in ind){
  Origin[[i]] <- predict(Models[[i]], newdata=new, predType = "full-model", se.fit = FALSE)
}


# Set maximal dose 
mD <- 1000


# Set seed
set.seed(167)





# Run each simulation separately because of long runtime 

# First for locally D-optimal designs 

# Weights of the support points 
weights <- rep(1/4,4)

# Initialize variable to save results 
results <- vector(mode="list", length = length(relPar))
g_OD <- vector(mode="list", length = 6)
df_results <- data.frame(sim = 1:500, EC50=numeric(500), NRMSE=numeric(500))

for(i in 1:6){
  g_OD[[i]] <- df_results
}

names(g_OD) <- c("opt 18","opt 27","opt 36","opt 45","opt 63","opt 90")

for(i in 1:length(results)){
  results[[i]] <- g_OD
}

# Calculate runtime 
starttime <- Sys.time()

# Start simulation and vary genes in i 
for(i in ind){
  # Fix model 
  mod <- Models[[i]]
  
  # Extract original model course (outsourced because of runtime)
  origin <- Origin[[i]]
  
  # Use model specific sigma  
  sigma2 <- sigma[i]
  
  # Extract support points 
  supp <- as.numeric(sort(optDesigns[[i]]$supPoints))
  
  # Simulation for each observation number varied with j 
  for(j in 1:6){
    # Use rounding procedure of Pukelsheim to fix numbers of observations at each step
    reps <- rndDesign(weights, n_opt[j])
    doses <- rep(supp, reps)
    
    newdat <- data.frame(dose = doses, iwas=factor(1))
    truePred <- predict(mod, newdata=newdat, predType = "full-model", se.fit = FALSE)
    
    # Simulation steps varied in i 
    for(z in 1:simnum){
      # Generate normal errors with corresponding sigma 
      error <- rnorm(length(doses),0,sigma2)
      simdata <- truePred+error
      simdat <- data.frame(dose=doses, resp=simdata)
      
      # Fit model with simulated data points 
      simmod <- fitMod(dose,resp, data=simdat, 
                       model="sigEmax", 
                       bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
      
      # Calculate EC50 
      results[[i]][[j]]$EC50[z] <- unname(simmod$coefs[3])
      
      # Calculate RMSE 
      newPred <- predict(simmod, newdata=new, predType = "full-model", se.fit = FALSE)
      # Standardize with corresponding range (EMax parameter) of each gene to get NRMSE
      results[[i]][[j]]$NRMSE[z] <- rmse(origin, newPred)/as.numeric(abs(mod$coefs[2]))
    }
  }
}
endtime <- Sys.time()
time_opt <- endtime-starttime # Runs approx. 2 days 

opt <- results[ind] 


#######

# Do the same for the other designs 

# Simultaneous D-optimal design: 

# Initialize variable to save results 
results <- vector(mode="list", length = length(relPar))
g_OD <- vector(mode="list", length = 6)
df_results <- data.frame(sim = 1:500, EC50=numeric(500), NRMSE=numeric(500))

for(i in 1:6){
  g_OD[[i]] <- df_results
}

names(g_OD) <- c("simultaneous 18","simultaneous 27","simultaneous 36","simultaneous 45",
                 "simultaneous 63","simultaneous 90")

for(i in 1:length(results)){
  results[[i]] <- g_OD
}


# calculate simultaneous inference design for different numbers of measurements with help of 
# rounding procedure of Pukelsheim
set.seed(345)
ODsimultaneous <- vector(mode="list", length=6)
reps <- c(2,3,4,5,7,10)

weights <- simultaneousDesign$weights
supp <- simultaneousDesign$supPoints

reps2 <- vector(mode="list", length=length(reps))
ODsimultaneous <- vector(mode="list", length=length(reps))

for(i in 1:length(reps)){
  reps2[[i]] <- rndDesign(weights, n_opt[i])
  ODsimultaneous[[i]] <- rep(supp, reps2[[i]])
}



# Calculate runtime 
starttime <- Sys.time()

# Start simulation and vary genes in i (only look at genes in ind)

for(i in ind){
  # Fix model 
  mod <- Models[[i]]
  
  # Extract original model course (outsourced because of runtime)
  origin <- Origin[[i]]
  
  # Use model specific sigma  
  sigma2 <- sigma[i]
  
  # Vary number of measurements  
  for(j in 1:6){
    
    doses <- ODsimultaneous[[j]]
    newdat <- data.frame(dose = doses, iwas=factor(1))
    truePred <- predict(mod, newdata=newdat, predType = "full-model", se.fit = FALSE)
    
    # Simulation steps varied in i 
    for(z in 1:simnum){
      # Generate normal errors with corresponding sigma 
      error <- rnorm(length(doses),0,sigma2)
      simdata <- truePred+error
      simdat <- data.frame(dose=doses, resp=simdata)
      
      # Fit model with simulated data points 
      simmod <- fitMod(dose,resp, data=simdat, 
                       model="sigEmax", 
                       bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
      
      # Calculate EC50 
      results[[i]][[j]]$EC50[z] <- unname(simmod$coefs[3])
      
      # Calculate RMSE 
      newPred <- predict(simmod, newdata=new, predType = "full-model", se.fit = FALSE)
      # Standardize with corresponding range (EMax parameter) of each gene to get NRMSE
      results[[i]][[j]]$NRMSE[z] <- rmse(origin, newPred)/as.numeric(abs(mod$coefs[2]))
    }
  }
}
endtime <- Sys.time()
time_simultaneous <- endtime-starttime # Runs approx. 2 days 

simultaneous <- results[ind]


###

# Log-equidistant 

# Initialize variable to save results 
results <- vector(mode="list", length = length(relPar))
g_OD <- vector(mode="list", length = 6)
df_results <- data.frame(sim = 1:500, EC50=numeric(500), NRMSE=numeric(500))

for(i in 1:6){
  g_OD[[i]] <- df_results
}

names(g_OD) <- c("logequi 18","logequi 27","logequi 36","logequi 45","logequi 63","logequi 90")

for(i in 1:length(results)){
  results[[i]] <- g_OD
}

# First initialize logequidistant designs with different repetitions of measurements 
ODlogequi <- vector(mode="list",length=6)
reps <- c(2,3,4,5,7,10)

for(i in 1:6){
  ODlogequi[[i]] <- rep(logequiDesign$supPoints,each=reps[i])
}

# Set seed again
set.seed(167)

# Calculate runtime 
starttime <- Sys.time()

# Start simulation and vary genes in i (only look at genes in ind)

for(i in ind){
  # Fix model 
  mod <- Models[[i]]
  
  # Extract original model course (outsourced because of runtime)
  origin <- Origin[[i]]
  
  # Use model specific sigma  
  sigma2 <- sigma[i]
  
  # Vary number of measurements  
  for(j in 1:6){
    
    doses <- ODlogequi[[j]]
    newdat <- data.frame(dose = doses, iwas=factor(1))
    truePred <- predict(mod, newdata=newdat, predType = "full-model", se.fit = FALSE)
    
    # Simulation steps varied in i 
    for(z in 1:simnum){
      # Generate normal errors with corresponding sigma 
      error <- rnorm(length(doses),0,sigma2)
      simdata <- truePred+error
      simdat <- data.frame(dose=doses, resp=simdata)
      
      # Fit model with simulated data points 
      simmod <- fitMod(dose,resp, data=simdat, 
                       model="sigEmax", 
                       bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
      
      # Calculate EC50 
      results[[i]][[j]]$EC50[z] <- unname(simmod$coefs[3])
      
      # Calculate RMSE 
      newPred <- predict(simmod, newdata=new, predType = "full-model", se.fit = FALSE)
      # Standardize with corresponding range (EMax parameter) of each gene to get NRMSE
      results[[i]][[j]]$NRMSE[z] <- rmse(origin, newPred)/as.numeric(abs(mod$coefs[2]))
    }
  }
}
endtime <- Sys.time()
time_logequi <- endtime-starttime # Runs approx. 2 days 

logequi <- results[ind]



######

# Equidistant design 

# Initialize variable to save results 
results <- vector(mode="list", length = length(relPar))
g_OD <- vector(mode="list", length = 6)
df_results <- data.frame(sim = 1:500, EC50=numeric(500), NRMSE=numeric(500))

for(i in 1:6){
  g_OD[[i]] <- df_results
}

names(g_OD) <- c("equi 18","equi 27","equi 36","equi 45","equi 63","equi 90")

for(i in 1:length(results)){
  results[[i]] <- g_OD
}

# First initialize equidistant designs with different repetitions of measurements 
ODequi <- vector(mode="list",length=6)
reps <- c(2,3,4,5,7,10)

for(i in 1:6){
  ODequi[[i]] <- rep(equiDesign$supPoints,each=reps[i])
}

# Set seed again
set.seed(167)

# Calculate runtime 
starttime <- Sys.time()

# Start simulation and vary genes in i (only look at genes in ind)

for(i in ind){
  # Fix model 
  mod <- Models[[i]]
  
  # Extract original model course (outsourced because of runtime)
  origin <- Origin[[i]]
  
  # Use model specific sigma  
  sigma2 <- sigma[i]
  
  # Vary number of measurements  
  for(j in 1:6){
    
    doses <- ODequi[[j]]
    newdat <- data.frame(dose = doses, iwas=factor(1))
    truePred <- predict(mod, newdata=newdat, predType = "full-model", se.fit = FALSE)
    
    # Simulation steps varied in i 
    for(z in 1:simnum){
      # Generate normal errors with corresponding sigma 
      error <- rnorm(length(doses),0,sigma2)
      simdata <- truePred+error
      simdat <- data.frame(dose=doses, resp=simdata)
      
      # Fit model with simulated data points 
      simmod <- fitMod(dose,resp, data=simdat, 
                       model="sigEmax", 
                       bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
      
      # Calculate EC50 
      results[[i]][[j]]$EC50[z] <- unname(simmod$coefs[3])
      
      # Calculate RMSE 
      newPred <- predict(simmod, newdata=new, predType = "full-model", se.fit = FALSE)
      # Standardize with corresponding range (EMax parameter) of each gene to get NRMSE
      results[[i]][[j]]$NRMSE[z] <- rmse(origin, newPred)/as.numeric(abs(mod$coefs[2]))
    }
  }
}
endtime <- Sys.time()
time_equi <- endtime-starttime # runs approx. 2 days 

equi <- results[ind]


######

# Original design 

# Initialize variable to save results 
results <- vector(mode="list", length = length(relPar))
g_OD <- vector(mode="list", length = 6)
df_results <- data.frame(sim = 1:500, EC50=numeric(500), NRMSE=numeric(500))

for(i in 1:6){
  g_OD[[i]] <- df_results
}

names(g_OD) <- c("original 18","original 27","original 36","original 45",
                 "original 63","original 90")

for(i in 1:length(results)){
  results[[i]] <- g_OD
}

# First initialize original designs with different repetitions of measurements 
ODoriginal <- vector(mode="list",length=6)
reps <- c(2,3,4,5,7,10)

# Initialize original designs with different numbers of measurements 
for(i in 1:6){
  ODoriginal[[i]] <- rep(c(0,origDesign$supPoints),each=reps[i])
}

# Set seed again
set.seed(167)

# Calculate runtime 
starttime <- Sys.time()

# Start simulation and vary genes in i (only look at genes in ind)

for(i in ind){
  # Fix model 
  mod <- Models[[i]]
  
  # Extract original model course (outsourced because of runtime)
  origin <- Origin[[i]]
  
  # Use model specific sigma  
  sigma2 <- sigma[i]
  
  # Vary number of measurements  
  for(j in 1:6){
    
    doses <- ODoriginal[[j]]
    newdat <- data.frame(dose = doses, iwas=factor(1))
    truePred <- predict(mod, newdata=newdat, predType = "full-model", se.fit = FALSE)
    
    # Simulation steps varied in i 
    for(z in 1:simnum){
      # Generate normal errors with corresponding sigma 
      error <- rnorm(length(doses),0,sigma2)
      simdata <- truePred+error
      simdat <- data.frame(dose=doses, resp=simdata)
      
      # Fit model with simulated data points 
      simmod <- fitMod(dose,resp, data=simdat, 
                       model="sigEmax", 
                       bnds = matrix(c(0.001 * mD, 0.05, 1.5 * mD, 10), 2))
      
      # Calculate EC50 
      results[[i]][[j]]$EC50[z] <- unname(simmod$coefs[3])
      
      # Calculate RMSE 
      newPred <- predict(simmod, newdata=new, predType = "full-model", se.fit = FALSE)
      # Standardize with corresponding range (EMax parameter) of each gene to get NRMSE
      results[[i]][[j]]$NRMSE[z] <- rmse(origin, newPred)/as.numeric(abs(mod$coefs[2]))
    }
  }
}
endtime <- Sys.time()
time_original <- endtime-starttime # Runs approx. 2 days 

original <- results[ind]

