library(fGarch)
library(Ecdat)
library(sn)
library(MASS)
library(LambertW)
library(ks)
library(htmlTable)
library(copula)
library(ggplot2)
library(patchwork)
library(pastecs)
library(rugarch)
library(forecast)
library(quantmod)
library(tseries)
library(modelsummary)

## download historical data
getSymbols("SPXL;TFC",
           from = "2012/11/01",
           to = "2022/11/01",
           periodicity = "daily")

## log return conversion
SPXLr<-diff(log(SPXL$SPXL.Adjusted))%>%
  na.omit()
names(SPXLr)<-c("logreturns")

TFCr<-diff(log(TFC$TFC.Adjusted))%>%
  na.omit()
names(TFCr)<-c("logreturns")

mean(SPXLr)
sd(SPXLr)

plot(density(SPXLr$logreturns),lwd = 2,"Density Plot of SPXL Log Returns") + polygon(density(SPXLr$logreturns))

## historical stock price and log returns
par(mfrow = c(2, 1))
plot(SPXL$SPXL.Adjusted, main = "SPXL Adjusted Stock Prices")
plot(SPXLr, main = "SPXL Daily Log Returns")

plot(SPXLr^2, ylim=c(0,.08),main = "SPXL Return Volatility")

par(mfrow = c(1, 2))

acf(SPXLr, main= "ACF of SPXL Log Returns")
acf(SPXLr^2, main= "ACF of SPXL Squared Log Returns")
Box.test(SPXLr, lag=10, type = 'Ljung')
adf.test(SPXLr)
Box.test(SPXLr^2, type = 'Ljung')

mean(TFCr)
sd(TFCr)

plot(density(TFCr$logreturns),lwd = 2,"Density Plot of TFC Log Returns") + polygon(density(TFCr$logreturns))
## historical stock price and log returns
par(mfrow = c(2, 1))
plot(TFC$TFC.Adjusted, main = "TFC Adjusted Stock Prices")
plot(TFCr, main = "TFC Daily Log Returns")

plot(TFCr^2, ylim=c(0,.05),main = "TFC Return Volatility")

par(mfrow = c(1, 2))

acf(TFCr, main= "ACF of TFC Log Returns")
acf(TFCr^2, main= "ACF of TFC Squared Log Returns")

Box.test(SPXLr, lag=10, type = 'Ljung')
adf.test(SPXLr)
Box.test(SPXLr^2, type = 'Ljung')

arma.garch = ugarchspec(mean.model=list(armaOrder=c(1,0)),
                        variance.model=list(garchOrder=c(1,1)),
                        distribution.model = "std")

spxl.ar.garch = ugarchfit(data=SPXLr, spec=arma.garch)
show(spxl.ar.garch)


spxle<-residuals(spxl.ar.garch, standardize=TRUE)
spxlfit_t<-fitdistr(spxle,"t")

spxl_para<-spxlfit_t$estimate

u1<-pstd(spxle,mean = spxl_para[1],sd = spxl_para[2],nu = spxl_para[3])%>%
  as.vector()


tfc.ar.garch = ugarchfit(data=TFCr, spec=arma.garch)
show(tfc.ar.garch)


tfce<-residuals(tfc.ar.garch, standardize=TRUE)

tfcfit_t<-fitdistr(tfce,"t")
tfc_para<-tfcfit_t$estimate
u2<-pstd(tfce,mean = tfc_para[1],sd = tfc_para[2],nu = tfc_para[3])%>%
  as.vector()

par(mfrow = c(1, 2))
plot(spxl.ar.garch, which = 10)
plot(spxl.ar.garch, which = 11)
mtext("SPXL",side = 3, line = - 2, outer = TRUE, font = 2)
par(mfrow = c(1, 1))
plot(spxl.ar.garch, which = 8)

Box.test(spxle, type = 'Ljung')
Box.test(spxle^2, type = 'Ljung')

par(mfrow = c(1, 2))
plot(tfc.ar.garch, which = 10)
plot(tfc.ar.garch, which = 11)
mtext("TFC",side = 3, line = - 2, outer = TRUE, font = 2)
par(mfrow = c(1, 1))
plot(tfc.ar.garch, which = 8)
Box.test(tfce, type = 'Ljung')
Box.test(tfce^2, type = 'Ljung')


U.hat<-cbind(u1,u2)
### Non-Parametric Density estimation plot
fhatU=kde(x=U.hat,H=Hscv(x=U.hat))
plot(fhatU,cont=seq(10,80,10), main = "Contour plot",xlab = "SPXL Standardized Residuals", ylab = "TFC Standardized Residuals")


tau=cor.test(as.numeric(u1),as.numeric(u2),method="kendall")$estimate
r=cor.test(as.numeric(u1),as.numeric(u2),method="pearson")$estimate
## estimator for spearman's rho used in copulas
omega=sin((tau*pi)/2)
## Guass
Cgauss <- fitCopula(copula=normalCopula(dim=2),data=U.hat,method="ml",start=c(omega))
guass_para<-coef(Cgauss)

guass_loglike <- loglikCopula(param=guass_para,u=U.hat,copula=normalCopula(dim=2))
guass_AIC <- -2*guass_loglike + 2*length(guass_para)

## T
Ct<-fitCopula(copula=tCopula(dim=2),data=U.hat,method="ml",start=c(omega,10))
t_para<-coef(Ct)

t_loglik<-loglikCopula(param=t_para,u=U.hat,copula=tCopula(dim=2))
t_AIC<- -2*t_loglik + 2*length(t_para)

## Clayton
Cclay <- fitCopula(copula=claytonCopula(dim=2),data=U.hat,method="ml",start=1) 
clay_para <- coef(Cclay)
clay_loglik <- loglikCopula(param=clay_para,u=(U.hat),copula=claytonCopula(dim=2))

clay_AIC<- -2*clay_loglik + 2*length(clay_para)


## Gumbel
Cgum <- fitCopula(copula=gumbelCopula(dim=2),data=U.hat,method="ml",start=1) 
gum_para<-coef(Cgum)

gum_loglik <- loglikCopula(param=gum_para,u=as.matrix(U.hat),copula=gumbelCopula(dim=2)) 
gum_AIC <- -2*gum_loglik + 2*length(gum_para)

Copula<-matrix(nrow = 4, ncol = 1)
Copula[1,]=c(guass_AIC)
Copula[2,]=c(t_AIC) 
Copula[3,]=c(clay_AIC) 
Copula[4,]=c(gum_AIC) 

Copula<-data.frame(Copula)
names<-c("Guassian", "T", "Clayton", "Gumbel")
row.names(Copula) <- names
names(Copula) <- c("AIC")
Copula<-round(Copula,2)

htmlTable(Copula, caption = "Copula Results")
set.seed(1234)

#generate a random sample with size n from the fitted t copula
Sim_prob <- rCopula(copula = tCopula(t_para[1], df = t_para[2], dim = 2), n = 10000)
Sim_prob <-data.frame(Sim_prob)
names(Sim_prob)<-c("SPXL","TFC")

#transform marginals into the estimated t
Sim_SPXLe <- qstd(Sim_prob[,1],mean = spxl_para[1],sd = spxl_para[2],nu = spxl_para[3])

Sim_TFCe <- qstd(Sim_prob[,2],mean = tfc_para[1],sd = tfc_para[2],nu = tfc_para[3])

## sigma for next period
spxl_model<-coef(spxl.ar.garch)
tfc_model<-coef(tfc.ar.garch)

sigma1spxl<- sqrt(spxl_model[3]+(spxl_model[4]*(tail(spxle,1)^2))+(spxl_model[5]*(tail(sigma(spxl.ar.garch),1)^2
)))%>%as.numeric()

sigma1tfc<- sqrt(tfc_model[3]+(tfc_model[4]*(tail(tfce,1)^2))+(tfc_model[5]*(tail(sigma(tfc.ar.garch),1)^2
)))%>%as.numeric()

## log return
spxlr0 <- as.numeric(tail(SPXLr,1))
spxl1<-spxl_model[1]+(spxl_model[2]*(spxlr0))+(sigma1spxl*Sim_SPXLe)

tfcr0 <- as.numeric(tail(TFCr,1))
tfc1<-tfc_model[1]+(tfc_model[2]*(tfcr0))+(sigma1tfc*Sim_TFCe)


## store Var in table 
rho<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
n=length(rho)


vartable<-matrix(data=NA,nrow=9,ncol=2)
vartable[,1]<-rho

fill= rep(NA,9)
for(i in 1:n){
  fill[i]<- (rho[i]*(spxl1)+((1-rho[i])*(tfc1)))%>%
    quantile(0.01)*-1
  vartable[,2]<-fill
}
vartable<-data.frame(vartable)
names(vartable)<-c('SPXL Weight','VaR')
vartable<-round(vartable,3)

htmlTable(vartable, caption = "VaR Table")
