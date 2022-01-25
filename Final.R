

rm(list = ls())
library(quadprog)
library(xts)
library(zoo)
library(quantmod)
library(lubridate)
library(PerformanceAnalytics)
library(rpart)
library(writexl)
library(xlsx)
library(Matrix)

setwd("/Users/metuhead/Desktop/Graduate/Senior/FE630/Final/Data/Data_Cleaned")

################################################################################
## 1.Data Handling

#1) DownLoad Data

tickers <- c("FXE","EWJ","GLD","QQQ","SPY","SHV","DBA","USO","XBI","ILF","EPP","FEZ")

getSymbols(tickers,from="2007-02-28", to="2021-10-31")

# Download the Vix index for additional project(Endogenic Target Beta and Weight)
getSymbols("^VIX",from="2007-02-28", to="2021-10-31")
V <- VIX$VIX.Close

P <- merge(FXE$FXE.Close,EWJ$EWJ.Close,GLD$GLD.Close,QQQ$QQQ.Close,SPY$SPY.Close,
           SHV$SHV.Close,DBA$DBA.Close,USO$USO.Close,XBI$XBI.Close,ILF$ILF.Close,
           EPP$EPP.Close,FEZ$FEZ.Close)

R <- merge(dailyReturn(FXE$FXE.Close),dailyReturn(EWJ$EWJ.Close),dailyReturn(GLD$GLD.Close),
           dailyReturn(QQQ$QQQ.Close),dailyReturn(SPY$SPY.Close),dailyReturn(SHV$SHV.Close),
           dailyReturn(DBA$DBA.Close),dailyReturn(USO$USO.Close),dailyReturn(XBI$XBI.Close),
           dailyReturn(ILF$ILF.Close),dailyReturn(EPP$EPP.Close),dailyReturn(FEZ$FEZ.Close))

rm(FXE,EWJ,GLD,QQQ,SPY,SHV,DBA,USO,XBI,ILF,EPP,FEZ,VIX)

colnames(P) <- tickers
colnames(R) <- tickers
R <- R[-1,]
Date <- index(R)
Fa <- read.csv("Fa.csv")
Fa <- Fa[,-1]
rownames(Fa) <- Date
colnames(Fa) <- c("RMRF","SMB","HML","RF")
Fa <- as.data.frame(Fa)

#2) Check the yearly Return and covariance(3 Period)

R_year <- as.matrix(apply.yearly(R,mean) *252)
rownames(R_year) <- c("07","08","09","10","11","12","13","14","15","16","17","18","19","20","21")
R_year

write.xlsx(R_year,"R_year.xlsx")

BC <- R["2007-03-01/2008-08-28",]
DC <- R["2008-09-01/2010-09-01",]
AC <- R["2010-09-01/2021-10-31",]


Q_BC <- cor(BC)
Q_DC <- cor(DC)
Q_AC <- cor(AC)


write.xlsx(Q_BC,"Q_BC.xlsx")
write.xlsx(Q_DC,"Q_DC.xlsx")
write.xlsx(Q_AC,"Q_AC.xlsx")

#3) Make some matrix(Beta, Mu, RiRf, Coef)

FD <- Fa/252 # FD is daily of Fama factors by dividing 252
RiRf <- R-FD[,"RF"]
Coef <- as.data.frame(matrix(0,nrow=12,ncol=4))
Beta <- as.data.frame(matrix(0,nrow=12,ncol=1))

colnames(Coef) <- c("Alpha","RMRF","SMB","HML")
colnames(Beta) <- c("Beta")
rownames(Coef) <- tickers
rownames(Beta) <- tickers
FD1 <- cbind(rep(1,length(FD[,1])),FD[,1:3])
colnames(FD1) <- c("Alpha","RMRF","SMB","HML")

Mu <- as.data.frame(matrix(0,nrow=12,ncol=1))
colnames(Mu) <- c("Mean") 


################################################################################
## 2.ProcessData function

processdata <- function(LB_Mu,LB_Q,startday){

  # 1. Set the data for getting Mu and Q and Beta

  point <- which(Date == startday)
  LBst_M <- point - LB_Mu
  LBed_M <- point - 1
  LBst_Q <- point - LB_Q
  LBed_Q <- point - 1
  FD_M <- FD[LBst_M:LBed_M,]
  FD1_M <- FD1[LBst_M:LBed_M,]
  R_M <- R[LBst_M:LBed_M,]
  R_Q <- R[LBst_Q:LBed_Q,]
  RiRf_M <- RiRf[LBst_M:LBed_M,]


  #2 Regression of Fama 3 factor model & CAPM(Getting Beta)

  for (i in 1:12) {
  
    Coef[i,] <- coefficients(lm(RiRf_M[,i] ~ FD_M[,"RMRF"]+FD_M[,"SMB"]+FD_M[,"HML"]))
    Beta[i,] <- coefficients(lm(RiRf_M[,i] ~ RiRf_M[,"SPY"]))[2]
  
    }

  #3 Calculate the Mu, Q
  
  Rho <- as.data.frame(t(as.matrix(Coef) %*% t(as.matrix(FD1_M))))
  Mu <- as.data.frame(sapply(Rho,mean))
  colnames(Mu) <- c("Mean")
  
  Q <- as.data.frame(cov(Rho))
  
  save(Mu,Q,Beta,file="inputs.rData")
  }

processdata(60,60,"2008-09-15")

load("inputs.rData")
Mu
Q
Beta

################################################################################
## 3. Optimizer function

Optimizer <- function(lambda,T_Beta,wp) {

  load("Inputs.rData")
  
  if (det(as.matrix(2*lambda*Q)) < 0.001) {
    
    Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
  
  } else {
    
    Dmat <- as.matrix(2*lambda*Q)
  } 
  
  dvec <- as.matrix(Mu) + as.matrix(2* lambda * Q) %*% as.matrix(wp)
  Amat <- as.matrix(cbind(rep(1,12),Beta,diag(12),-diag(12)))
  bvec <- c(1,T_Beta,rep(-2,24))
  
  w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
  rownames(w) <- tickers
  
  return(w)
  
} 

wp <- matrix(c(0.1,0.2,-0.3,0.1,0.2,-0.1,0.2,0.4,-0.2,0.5,-0.1,0),ncol=1)
Optimizer(0.5,2,wp)
sum(wp)

###############################################################################
## 4. Realized Portfolio Return

Port_Return <- function(startday,endday,LB_M,LB_Q,T_Beta,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex    # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex
  
  Nrebal <- PlayTime/5
  
  Date_RR <- index(R[startindex:(endindex-1),])
  RR <- as.matrix(R[startindex:(endindex-1),]) # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/12,12),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  
  ret <- NULL # ret is the daily return for 5days
  Rret <- ret # Rret stored the realized return (Daily)
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  for (i in 1:Nrebal) {
    
  
    processdata(LB_M,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    wp <- Optimizer(lambda,T_Beta,wp,startday) # Get new portfolio wp from the Beta, Mu, Q
    W <- rbind(W,t(wp)) # Record the new portfolio
    
    ret <- RR[(5*(i-1)+1):(5*(i-1)+5),] %*% wp     # Make the Realized Return x weights
    Rret <- rbind(Rret,ret)    # Record the Realized Return reflecting the weights
    
    startindex <- startindex + 5
    startday <- as.Date(Date[startindex])
  }
  
  Rret <- xts(Rret,order.by = Date_RR)
  colnames(Rret) <- c("Return")

  save(W, Rret, file="Result.rData")
  
  }


############# Check everything is OK ###########################################


## Before the excution of the function. 
## Please, check your startday and endday if the market is open
## if startindex or endindex do not have a certain number, the market is closed

startday <- "2008-01-02"
startindex <- which(Date == startday)
startindex

endday <- "2009-01-02"
endindex <- which(Date == endday)
endindex

N <- (endindex - ((endindex - startindex) %% 5) - startindex) / 5
N

Port_Return(startday,endday,60,60,1,0.5)
load("Result.rData")


plot(W[,"GLD"],type="o",ylab="GLD&QQQ",col="orange",main="Weights on GLD and QQQ")
points(W[,"QQQ"],type="o",col="Blue")

sapply(as.data.frame(W),mean)
apply.monthly(Rret,mean)

Rret[232:250]

chart.CumReturns(Rret)

R["2008-12-10"]

###############################################################################
## 5. Analysis

# 1. Divide the period
# We assume the 2008.9.1 is the start of the Crisis when Leman Brothers defaulted
# and it ends 2010.9.1.


BC <- R["2007-03-01/2008-08-28",]
DC <- R["2008-09-01/2010-09-01",]
AC <- R["2010-09-01/2021-10-31",]

Mean_ETF <- NULL
Mean_ETF <- rbind(Mean_ETF,Mean_BC,Mean_DC,Mean_AC)
Mean_BC <-sapply(BC,mean) * 252
Mean_DC <- sapply(DC,mean) * 252
Mean_AC <- sapply(AC,mean) * 252

rownames(Mean_ETF) <- c("Before","During","After")

# BC is Before Crisis, DC is during Crisis, AC is After Crisis

# The Period before Crisis, the investment starts at 07-07-16 ends at 08-08-28

chart.CumReturns(BC[,"SPY"])
chart.CumReturns(DC[,"SPY"])
chart.CumReturns(AC[,"SPY"])

####################### Before Crisis ##########################################

startday <- "2007-07-16"
startindex <- which(Date == startday)
startindex

endday <- "2008-01-31"
endindex <- which(Date == endday)
endindex

startdayBC <- "2007-07-16"
enddayBC <- "2008-01-31"

BC_30.90.0.5 <- Port_Return(startdayBC,enddayBC,30,90,0.5,1)
BC_60.90.0.5 <- Port_Return(startdayBC,enddayBC,60,90,0.5,1)
BC_90.90.0.5 <- Port_Return(startdayBC,enddayBC,90,90,0.5,1)

BC_SPY <- R["2007-07-16/2008-01-25"][,"SPY"]
BC_R <- R["2007-03-01/2008-08-28"]

sum(BC_SPY)

sapply(BC_R,mean) * 252

R_BC <- cbind(BC_SPY,BC_30.90.0.5,BC_60.90.0.5,BC_90.90.0.5)
colnames(R_BC) <- c("SPY","30/90","60/90","90/90") 

Mean_BC <- sapply(R_BC,mean)*252 # Annualized Mean Return Before Crisis
Vol_BC <- sapply(R_BC,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_BC <- sapply(R_BC,skewness)
Kur_BC <- sapply(R_BC,kurtosis)
SR_BC <- sapply(R_BC,SharpeRatio.annualized)

Table_BC <- rbind(Mean_BC,Vol_BC,Ske_BC,Kur_BC,SR_BC)
rownames(Table_BC) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_BC

chart.CumReturns(R_BC,wealth.index = TRUE,colorset=c("Blue","orange","green","red"),legend.loc = "topleft")




BC_60.90.m1 <- Port_Return(startdayBC,enddayBC,60,90,-1,1)
BC_60.90.m0.5 <- Port_Return(startdayBC,enddayBC,60,90,-0.5,1)
BC_60.90.0 <- Port_Return(startdayBC,enddayBC,60,90,0,1)
BC_60.90.0.5 <- Port_Return(startdayBC,enddayBC,60,90,0.5,1)
BC_60.90.1 <- Port_Return(startdayBC,enddayBC,60,90,1,1)

R_BC2 <- cbind(BC_SPY,BC_60.90.m1,BC_60.90.m0.5,BC_60.90.0,BC_60.90.0.5,BC_60.90.1)
colnames(R_BC2) <- c("SPY","beta -1","beta -0.5","beta 0","beta 0.5", "beta 1") 

Mean_BC2 <- sapply(R_BC2,mean)*252 # Annualized Mean Return Before Crisis
Vol_BC2 <- sapply(R_BC2,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_BC2 <- sapply(R_BC2,skewness)
Kur_BC2 <- sapply(R_BC2,kurtosis)
SR_BC2 <- sapply(R_BC2,SharpeRatio.annualized)

Table_BC2 <- rbind(Mean_BC2,Vol_BC2,Ske_BC2,Kur_BC2,SR_BC2)
rownames(Table_BC2) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_BC2

chart.CumReturns(R_BC2,colorset=c("Blue","orange","green","red","purple","black"),legend.loc = "topleft")

hist(BC_30.90.0.5,main="Before Crisis(30/90 histogram)",xlab="Beta = 0.5")
hist(BC_60.90.0,main="Before Crisis(60/90 histogram)",xlab="Beta = 0")



####################### During Crisis ##########################################

startday <- "2009-07-16"
startindex <- which(Date == startday)
startindex

endday <- "2010-07-15"
endindex <- which(Date == endday)
endindex

startdayDC <- "2009-07-16"
enddayDC <- "2010-07-15"

DC_30.90.0.5 <- Port_Return(startdayDC,enddayDC,30,90,0.5,1)
DC_60.90.0.5 <- Port_Return(startdayDC,enddayDC,60,90,0.5,1)
DC_90.90.0.5 <- Port_Return(startdayDC,enddayDC,90,90,0.5,1)

DC_SPY <- R["2009-07-16/2010-07-13"][,"SPY"]

R_DC <- cbind(DC_SPY,DC_30.90.0.5,DC_60.90.0.5,DC_90.90.0.5)
colnames(R_DC) <- c("SPY","30/90","60/90","90/90") 

Mean_DC <- sapply(R_DC,mean)*252 # Annualized Mean Return Before Crisis
Vol_DC <- sapply(R_DC,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_DC <- sapply(R_DC,skewness)
Kur_DC <- sapply(R_DC,kurtosis)
SR_DC <- sapply(R_DC,SharpeRatio.annualized)

Table_DC <- rbind(Mean_DC,Vol_DC,Ske_DC,Kur_DC,SR_DC)
rownames(Table_DC) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_DC

chart.CumReturns(R_DC,wealth.index = TRUE,colorset=c("Blue","orange","green","red"),legend.loc = "topleft")


DC_60.90.m1 <- Port_Return(startdayDC,enddayDC,60,90,-1,1)
DC_60.90.m0.5 <- Port_Return(startdayDC,enddayDC,60,90,-0.5,1)
DC_60.90.0 <- Port_Return(startdayDC,enddayDC,60,90,0,1)
DC_60.90.0.5 <- Port_Return(startdayDC,enddayDC,60,90,0.5,1)
DC_60.90.1 <- Port_Return(startdayDC,enddayDC,60,90,1,1)

R_DC2 <- cbind(DC_SPY,DC_60.90.m1,DC_60.90.m0.5,DC_60.90.0,DC_60.90.0.5,DC_60.90.1)
colnames(R_DC2) <- c("SPY","beta -1","beta -0.5","beta 0","beta 0.5", "beta 1") 

Mean_DC2 <- sapply(R_DC2,mean)*252 # Annualized Mean Return Before Crisis
Vol_DC2 <- sapply(R_DC2,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_DC2 <- sapply(R_DC2,skewness)
Kur_DC2 <- sapply(R_DC2,kurtosis)
SR_DC2 <- sapply(R_DC2,SharpeRatio.annualized)

Table_DC2 <- rbind(Mean_DC2,Vol_DC2,Ske_DC2,Kur_DC2,SR_DC2)
rownames(Table_DC2) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_DC2

chart.CumReturns(R_DC2,colorset=c("Blue","orange","green","red","purple","black"),legend.loc = "topleft")

hist(DC_30.90.0.5,main="During Crisis(30/90 histogram)",xlab="Beta = 0.5")
hist(DC_60.90.0.5,main="During Crisis(60/90 histogram)",xlab="Beta = 0.5")



####################### After Crisis ##########################################

startday <- "2012-07-16"
startindex <- which(Date == startday)
startindex

endday <- "2013-07-15"
endindex <- which(Date == endday)
endindex

startdayAC <- "2012-07-16"
enddayAC <- "2013-07-15"

AC_30.90.0.5 <- Port_Return(startdayAC,enddayAC,30,90,0.5,1)
AC_60.90.0.5 <- Port_Return(startdayAC,enddayAC,60,90,0.5,1)
AC_90.90.0.5 <- Port_Return(startdayAC,enddayAC,90,90,0.5,1)

AC_SPY <- R["2012-07-16/2013-07-08"][,"SPY"]

R_AC <- cbind(AC_SPY,AC_30.90.0.5,AC_60.90.0.5,AC_90.90.0.5)
colnames(R_AC) <- c("SPY","30/90","60/90","90/90") 

Mean_AC <- sapply(R_AC,mean)*252 # Annualized Mean Return Before Crisis
Vol_AC <- sapply(R_AC,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_AC <- sapply(R_AC,skewness)
Kur_AC <- sapply(R_AC,kurtosis)
SR_AC <- sapply(R_AC,SharpeRatio.annualized)

Table_AC <- rbind(Mean_AC,Vol_AC,Ske_AC,Kur_AC,SR_AC)
rownames(Table_AC) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_AC

chart.CumReturns(R_AC,wealth.index = TRUE,colorset=c("Blue","orange","green","red"),legend.loc = "topleft")


AC_60.90.m1 <- Port_Return(startdayAC,enddayAC,60,90,-1,1)
AC_60.90.m0.5 <- Port_Return(startdayAC,enddayAC,60,90,-0.5,1)
AC_60.90.0 <- Port_Return(startdayAC,enddayAC,60,90,0,1)
AC_60.90.0.5 <- Port_Return(startdayAC,enddayAC,60,90,0.5,1)
AC_60.90.1 <- Port_Return(startdayAC,enddayAC,60,90,1,1)

R_AC2 <- cbind(AC_SPY,AC_60.90.m1,AC_60.90.m0.5,AC_60.90.0,AC_60.90.0.5,AC_60.90.1)
colnames(R_AC2) <- c("SPY","beta -1","beta -0.5","beta 0","beta 0.5", "beta 1") 

Mean_AC2 <- sapply(R_AC2,mean)*252 # Annualized Mean Return Before Crisis
Vol_AC2 <- sapply(R_AC2,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_AC2 <- sapply(R_AC2,skewness)
Kur_AC2 <- sapply(R_AC2,kurtosis)
SR_AC2 <- sapply(R_AC2,SharpeRatio.annualized)

Table_AC2 <- rbind(Mean_DC2,Vol_DC2,Ske_DC2,Kur_DC2,SR_DC2)
rownames(Table_DC2) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_AC2

chart.CumReturns(R_AC2,colorset=c("Blue","orange","green","red","purple","black"),legend.loc = "topleft")

hist(AC_30.90.0.5,main="After Crisis(30/90 histogram)",xlab="Beta = 0.5")
hist(AC_60.90.0.5,main="After Crisis(60/90 histogram)",xlab="Beta = 0.5")


#################################################################################
## Additional Project : K-Surfing Vol Wave Strategy (Endogenic Target Beta and Weights)

V_week <- apply.weekly(V,mean)

hist(V_week,main="Histogram of Weekly Vix index",xlab="Vix")

V_Table <- matrix(rep(0,3),nrow=1,ncol=3)
V_Table[1] <- quantile(V_week,0.25)
V_Table[2] <- quantile(V_week,0.5)
V_Table[3] <- quantile(V_week,0.75)
colnames(V_Table) <- c("1Q","Median","3Q")
V_Table

K_Surf_Wave <- function(startday){
  
  point <- which(Date == startday)
  LBst_V <- point - 5
  LBed_V <- point - 1
  
  V_Mean <- mean(V[LBst_V:LBed_V])
  
  if (V_Mean < 13.78) {
    
    T_Beta <- 3
    
  } else if (V_Mean < 17.5) {
    
    T_Beta <- 1.5
    
  } else if (V_Mean < 22.5) {
  
    T_Beta <- 1
    
  } else {
    
    T_Beta <- -5
    
  }

  return(T_Beta)  
}


K_Optimizer <- function(lambda,T_Beta,wp) {
  
  load("Inputs.rData")
  
  if (det(as.matrix(2*lambda*Q)) < 0.001) {
    
    Dmat <- nearPD(as.matrix(2*lambda*Q))$mat
    
  } else {
    
    Dmat <- as.matrix(2*lambda*Q)
  } 
  
  dvec <- as.matrix(Mu) + as.matrix(2* lambda * Q) %*% as.matrix(wp)
  Amat <- as.matrix(cbind(rep(1,12),Beta,diag(12),-diag(12)))
  bvec <- c(1,T_Beta,rep(-2,24))
  
  w <- matrix(solve.QP(Dmat,dvec,Amat,bvec,meq=2)$sol,ncol=1)
  rownames(w) <- tickers
  
  if (T_Beta < -3) {
    
    w <- matrix(rep(0,12),ncol=1)
  }
  
  return(w)
  
} 


Port_Return_K <- function(startday,endday,LB_M,LB_Q,lambda){
  
  startindex <- which(Date == startday)
  endindex <- which(Date == endday)
  PlayTime <- endindex - startindex    # PlayTime is the period we are investing
  
  endindex <- endindex - (PlayTime %% 5)
  PlayTime <- endindex - startindex
  
  Nrebal <- PlayTime/5
  
  Date_RR <- index(R[startindex:(endindex-1),])
  RR <- as.matrix(R[startindex:(endindex-1),]) # RR is dailyReturn in PlayTime
  
  wp <- matrix(rep(1/12,12),ncol=1) # wp is the initial portfolio
  rownames(wp) <- tickers; colnames(wp) <- c("Weight")
  
  
  ret <- NULL # ret is the daily return for 5days
  Rret <- ret # Rret stored the realized return (Daily)
  W <- NULL  # W stores the weights set everyweek
  
  # Loop : Rebalance everyweek and find dailyreturn of the Portfolio
  
  for (i in 1:Nrebal) {
    
    
    processdata(LB_M,LB_Q,startday) # Get Beta, Mu, Q from the lookback period
    T_Beta <- K_Surf_Wave(startday)
    wp <- K_Optimizer(lambda,T_Beta,wp) # Get new portfolio wp from the Beta, Mu, Q
    W <- rbind(W,t(wp)) # Record the new portfolio
    
    ret <- RR[(5*(i-1)+1):(5*(i-1)+5),] %*% wp     # Make the Realized Return x weights
    Rret <- rbind(Rret,ret)    # Record the Realized Return reflecting the weights
    
    startindex <- startindex + 5
    startday <- Date[startindex]
  }
  
  Rret <- xts(Rret,order.by = Date_RR)
  colnames(Rret) <- c("DailyReturn")
  
  return(Rret)
  
}


a <- R["2015-01-01/2016-06-30"]
chart.CumReturns(a)
b <- V["2015-01-01/2016-06-30"]
plot(b)

startday <- "2007-07-13"
startindex <- which(Date == startday)
startindex

endday <- "2021-10-29"
endindex <- which(Date == endday)
endindex

A_6090B1 <- Port_Return(startday,endday,60,90,1,1)
A_6090B0.5 <- Port_Return(startday,endday,60,90,0.5,1)
A_6090B3 <- Port_Return(startday,endday,60,90,3,1)
A_6090B0 <- Port_Return(startday,endday,60,90,0,1)
A_6090_K <- Port_Return_K(startday,endday,60,90,1)


Return.cumulative(A_6090_K)
Return.cumulative(A_6090B0)
Return.cumulative(A_6090B0.5)
Return.cumulative(A_6090B1)
Return.cumulative(A_6090B3)


Kdata <- cbind(A_6090_K,A_6090B0,A_6090B0.5,A_6090B1,A_6090B3)

colnames(Kdata) <- c("K_Surf","beta=0","beta=0.5","beta=1","beta=3") 

Mean_K <- sapply(Kdata,mean)*252 # Annualized Mean Return Before Crisis
Vol_K <- sapply(Kdata,sd) * sqrt(252) # Annualized Volatility Before Crisis
Ske_K <- sapply(Kdata,skewness)
Kur_K <- sapply(Kdata,kurtosis)
SR_K <- sapply(Kdata,SharpeRatio.annualized)

Table_K <- rbind(Mean_K,Vol_K,Ske_K,Kur_K,SR_K)
rownames(Table_K) <- c("Mean","Volatility","Skewness","Kurtosis","Sharpe Ratio")
Table_K

chart.CumReturns(Kdata,wealth.index = TRUE,colorset=c("Blue","orange","green","red","purple"),legend.loc = "topleft")
























