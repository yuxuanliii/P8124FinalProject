################
# Load Libraries
################
library(micd)
library(pcalg)
source("/Users/lyx/Desktop/P8124/lv-ida-master/lvida.R")
source("/Users/lyx/Desktop/P8124/lv-ida-master/iscyclic.R")

#############
# Clean Data #
#############

load("/Users/lyx/Desktop/P8124/public_dataset/hsls_17_student_pets_sr_v1_0.rdata")

# Select variables that happens before college enrollment
vars <- c("X1SEX", "X1RACE", "X1DUALLANG", "X1TXMTH",
          "X1PAREDU", "X1PARPATTERN", "X1SES", "X1MTHID", "X1SCIID",
          "X1SCHOOLBEL", "X1SCHOOLENG", "X1STUEDEXPCT", "X1PAREDEXPCT",
          "X1TMCOMM", "X1TMEFF", "X1TMEXP", "X1TMPRINC", "X1TMRESP",
          "X1TSCOMM", "X1TSEFF", "X1TSEXP", "X1TSPRINC", "X1TSRESP",
          "X1SCHOOLCLI", "X1COUPERTEA", "X1COUPERCOU", "X1COUPERPRI",
          "X2ENROLSTAT", "X2EVERDROP", "X2TXMTH", "X2SES", "X2MTHID", 
          "X2SCIID", "X2STUEDEXPCT", "X2PAREDEXPCT", "X2REQLEVEL", 
          "X2SCHOOLCLI", "X3TGPAENG", "X3TCREDENG", "X3TCREDAPENG",
          "X3TCREDAPMTH", "X3TCREDMAT", "X3TGPAMAT", "X3TCREDAPSCI",
          "X3TCREDSCI", "X3TGPASCI", "X3TCREDSOCST", "X3TCREDAPSS", 
          "X3TGPASOCST", "X3TCREDART", "X3TCREDLANG", "X3TCREDAPLNG",
          "X3TGPALANG", "X4PSENRSTLV")
data <- hsls_17_student_pets_sr_v1_0[,c(vars, colnames(hsls_17_student_pets_sr_v1_0)[406:483])]

# Code missing as NA
data[data == -8 | data == -9 | data == -1] <- NA
# Remove observations with missing outcome: X4PSENRSTLV
data <- data[data[, "X4PSENRSTLV"] %in% c("Not enrolled", "Enrolled in a 4-yr institution", 
                                          "Enrolled in a 2-yr institution", 
                                          "Enrolled in a less-than-2-yr institution",
                                          "Enrolled, institutional level unknown"),]
# Transform categorical variables to factors
categoricalCols <- c(
  "X1SEX", "X1RACE", "X1DUALLANG", "X1PAREDU", "X1PARPATTERN",
  "X1STUEDEXPCT", "X1PAREDEXPCT",
  "X2ENROLSTAT", "X2EVERDROP", "X2STUEDEXPCT",
  "X2PAREDEXPCT", "X2REQLEVEL", "X4PSENRSTLV",
  colnames(data)[128:132]
)
data[, categoricalCols] <- lapply(data[, categoricalCols], as.factor)

# Remove columns with constant values
data <- data[, !colnames(data) %in% c("X3TGPAIB", "X3TGPAAP", "X3TCREDAP", "X3TCREDIB")]


###################################
# Select variables due to runtime #
###################################
varSelected <- c("X1SEX", "X1SES", "X1SCIID", "X3TGPA12TH", "X3TGPAAPIB", 
                 "X3TGPAENG", "X3TGPASCI", "X3TGPAMAT", "X4PSENRSTLV")

data_est <- data[,varSelected]
data_est <- data_est[data$X1SEX %in% c("Male", "Female"),]
data_est$X1SEX <- as.factor(ifelse(data_est$X1SEX == "Female", 1, 0))
Y <- data_est[complete.cases(data_est), "X4PSENRSTLV"]
Y_binary <- rep(1, length(Y))
Y_binary[Y  == "Not enrolled"] <- 0
data_est <- data_est[complete.cases(data_est),]
data_est <- cbind(as.factor(Y_binary), data_est[,1:(length(varSelected)-1)])
data_org <- data_est

#################################################################
# Run PC and FCI for 10 times for alpha = 0.01, 0.05, 0.1, 0.15 #
#################################################################

# alpha = 0.1

effects_fci <- list()
effects_pc <- list()

B = 10

for (b in 1:B) {
  half_idx <- sample(1:nrow(data_org), size = floor(nrow(data_org)/2), replace = FALSE)
  data_est <- data_org[half_idx, ]
  ## use mixCItest within pcalg::fci
  fci.fit <- fci(suffStat = data_est, indepTest = mixCItest, alpha = 0.1, p=9)
  pc.fit <- pc(suffStat = data_est, indepTest = mixCItest, alpha = 0.1, p=9)
  
  if(is.cyclic(fci.fit@amat)){
    cat("#### FOUND CYCLIC GRAPH #### \n LV-IDA won't work here! \n try again! \n")
  }
  data_est <- apply(data_est, 2, as.numeric)
  effects_fci[[b]] <- sapply(2:9, function(p) min(lv.ida(p,1,cov(data_est),fci.fit@amat,method="local"), na.rm=T))
  effects_pc[[b]] <- sapply(2:9, function(p) min(ida(p,1,cov(data_est), pc.fit@graph), na.rm=T))
}

# alpha = 0.15

effects_fci_015 <- list()
effects_pc_015 <- list()

for (b in 1:B) {
  half_idx <- sample(1:nrow(data_org), size = floor(nrow(data_org)/2), replace = FALSE)
  data_est <- data_org[half_idx, ]
  ## use mixCItest within pcalg::fci
  fci.fit <- fci(suffStat = data_est, indepTest = mixCItest, alpha = 0.15, p=9)
  pc.fit <- pc(suffStat = data_est, indepTest = mixCItest, alpha = 0.15, p=9)
  
  if(is.cyclic(fci.fit@amat)){
    cat("#### FOUND CYCLIC GRAPH #### \n LV-IDA won't work here! \n try again! \n")
  }
  data_est <- apply(data_est, 2, as.numeric)
  effects_fci_015[[b]] <- sapply(2:9, function(p) min(lv.ida(p,1,cov(data_est),fci.fit@amat,method="local"), na.rm=T))
  effects_pc_015[[b]] <- sapply(2:9, function(p) min(ida(p,1,cov(data_est), pc.fit@graph), na.rm=T))
  
}

# alpha = 0.05

effects_fci_005 <- list()
effects_pc_005 <- list()

for (b in 1:B) {
  half_idx <- sample(1:nrow(data_org), size = floor(nrow(data_org)/2), replace = FALSE)
  data_est <- data_org[half_idx, ]
  ## use mixCItest within pcalg::fci
  fci.fit <- fci(suffStat = data_est, indepTest = mixCItest, alpha = 0.05, p=9)
  pc.fit <- pc(suffStat = data_est, indepTest = mixCItest, alpha = 0.05, p=9)
  
  if(is.cyclic(fci.fit@amat)){
    cat("#### FOUND CYCLIC GRAPH #### \n LV-IDA won't work here! \n try again! \n")
  }
  data_est <- apply(data_est, 2, as.numeric)
  effects_fci_005[[b]] <- sapply(2:9, function(p) min(lv.ida(p,1,cov(data_est),fci.fit@amat,method="local"), na.rm=T))
  effects_pc_005[[b]] <- sapply(2:9, function(p) min(ida(p,1,cov(data_est), pc.fit@graph), na.rm = T))
  
}

# alpha = 0.01

effects_fci_001 <- list()
effects_pc_001 <- list()

for (b in 1:B) {
  half_idx <- sample(1:nrow(data_org), size = floor(nrow(data_org)/2), replace = FALSE)
  data_est <- data_org[half_idx, ]
  ## use mixCItest within pcalg::fci
  fci.fit <- fci(suffStat = data_est, indepTest = mixCItest, alpha = 0.01, p=9)
  pc.fit <- pc(suffStat = data_est, indepTest = mixCItest, alpha = 0.01, p=9)
  
  if(is.cyclic(fci.fit@amat)){
    cat("#### FOUND CYCLIC GRAPH #### \n LV-IDA won't work here! \n try again! \n")
  }
  data_est <- apply(data_est, 2, as.numeric)
  effects_fci_001[[b]] <- sapply(2:9, function(p) min(lv.ida(p,1,cov(data_est),fci.fit@amat,method="local"),na.rm = T))
  effects_pc_001[[b]] <- sapply(2:9, function(p) min(ida(p,1,cov(data_est), pc.fit@graph),na.rm = T))
  
}


# Results and Evaluation

top3_idx <- lapply(effects_fci_015, function(x) {
  order(x, decreasing = TRUE)[1:3]
})
all_idx <- unlist(top3_idx)

freq <- sort(table(all_idx), decreasing = TRUE)
top3_most_frequent <- head(freq, 3)





