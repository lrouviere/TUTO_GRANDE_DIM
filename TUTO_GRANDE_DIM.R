## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message = FALSE, size = "script",fig.height=3.5, fig.width=5,prompt=TRUE,collapse=TRUE,fig.align='center',cache=TRUE)
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message = FALSE, size = "script",prompt=FALSE,collapse=TRUE,fig.align='center',cache=FALSE,comment=NA)
knitr::knit_hooks$set(purl = knitr::hook_purl)

## ---- include=FALSE---------------------------------------
options(width = 60)
local({
  hook_output <- knitr::knit_hooks$get('output')
  knitr::knit_hooks$set(output = function(x, options) {
    if (!is.null(options$max.height)) options$attr.output <- c(
      options$attr.output,
      sprintf('style="max-height: %s;"', options$max.height)
    )
    hook_output(x, options)
  })
})

## ----message=FALSE, warning=FALSE, echo=FALSE-------------
library(tidyverse)
theme_set(theme_classic())

## ----include=FALSE,eval=FALSE-----------------------------
#  # automatically create a bib database for R packages
#  knitr::write_bib(c(
#    .packages(), 'bookdown', 'knitr', 'rmarkdown'
#  ), 'packages.bib')

## ---------------------------------------------------------
simu <- function(napp=300,ntest=500,p=3,graine=1234){
  set.seed(graine)
  n <- napp+ntest
  X <- matrix(runif(n*p),ncol=p)
  Y <- apply(X^2,1,sum)+rnorm(n,sd=0.5)
  Yapp <- Y[1:napp]
  Ytest <- Y[-(1:napp)]
  Xapp <- data.frame(X[1:napp,])
  Xtest <- data.frame(X[-(1:napp),])
  return(list(Xapp=Xapp,Yapp=Yapp,Xtest=Xtest,Ytest=Ytest))
}
df <- simu(napp=300,ntest=500,p=3,graine=1234)

## ---------------------------------------------------------
library(FNN)
mod3ppv <- knn.reg(train=df$Xapp,y=df$Yapp,k=3)

## ---------------------------------------------------------
mod3ppv$PRESS

## ---------------------------------------------------------
mod3ppv$PRESS/max(c(nrow(df$Xapp),1))

## ----teacher=correct--------------------------------------
sel.k <- function(K_cand=seq(1,50,by=5),Xapp,Yapp){
  ind <- 1
  err <- rep(0,length(K_cand))
  for (k in K_cand){
modkppv <- knn.reg(train=Xapp,y=Yapp,k=k)
err[ind] <- modkppv$PRESS/max(c(nrow(Xapp),1))
ind <- ind+1
  }
  return(K_cand[which.min(err)])
}

## ----eval=cor,echo=TRUE-----------------------------------
k.opt <- sel.k(seq(1,50,by=5),df$Xapp,df$Yapp)
prev <- knn.reg(train=df$Xapp,y=df$Yapp,test=df$Xtest,k=k.opt)$pred
mean((prev-df$Ytest)^2)

## ---------------------------------------------------------
DIM <- c(1,5,10,50)
K_cand <- seq(1,50,by=5)

## ----teacher=correct--------------------------------------
df <- data.frame(mat.err)
nom.dim <- paste("D",DIM,sep="")
names(df) <- nom.dim

## ----teacher=correct--------------------------------------
df %>% summarise_all(mean)
df %>% summarise_all(var)

## ----teacher=correct--------------------------------------
df1 <- pivot_longer(df,cols=everything(),names_to="dim",values_to="erreur")
df1 <- df1 %>% mutate(dim=fct_relevel(dim,nom.dim))
ggplot(df1)+aes(x=dim,y=erreur)+geom_boxplot()

## ---------------------------------------------------------
DIM <- c(0,50,100,200)

## ---------------------------------------------------------
n <- 250
p <- 1000
X <- matrix(runif(n*p),ncol=p)
simu.lin <- function(X,graine){
  set.seed(graine)
  Y <- X[,1]+rnorm(nrow(X),sd=0.5)
  df <- data.frame(Y,X)
  return(df)
}

## ----teacher=correct--------------------------------------
df <- data.frame(matbeta1)
nom.dim <- paste("D",DIM,sep="")
names(df) <- nom.dim

## ----teacher=correct--------------------------------------
df %>% summarise_all(mean)
df %>% summarise_all(var)

## ----teacher=correct--------------------------------------
df1 <- gather(df,key="dim",value="erreur")
df1 <- df1 %>% mutate(dim=fct_relevel(dim,nom.dim))
ggplot(df1)+aes(x=dim,y=erreur)+geom_boxplot()+theme_classic()

## ---------------------------------------------------------
ozone <- read.table("data/ozone.txt")
head(ozone)

## ----teacher=correct--------------------------------------
lin.complet <- lm(maxO3~.,data=ozone)
summary(lin.complet)
anova(lin.complet)

## ---------------------------------------------------------
library(leaps)
mod.sel <- regsubsets(maxO3~.,data=ozone,nvmax=14)
summary(mod.sel)

## ----teacher=correct--------------------------------------
plot(mod.sel,scale="r2")

## ----teacher=correct--------------------------------------
plot(mod.sel,scale="bic")
plot(mod.sel,scale="Cp")

## ----teacher=correct,indent='        '--------------------
mod.step <- step(lin.complet,direction="backward",trace=0)
mod.step

## ---------------------------------------------------------
library(ISLR)
Hitters <- na.omit(Hitters)

## ----teacher=correct--------------------------------------
X <- model.matrix(Salary~.,data=Hitters)[,-1]

## ----teacher=correct--------------------------------------
Xcr <- scale(X)
Xbar  <- apply(X,2,mean)
stdX <- apply(X,2,sd)

## ----teacher=correct--------------------------------------
library(FactoMineR)
acp.hit <- PCA(Xcr,scale.unit=FALSE,graph=TRUE)

## ----teacher=correct--------------------------------------
CC <- acp.hit$ind$coord

## ----teacher=correct--------------------------------------
donnees <- cbind.data.frame(CC,Salary=Hitters$Salary)
mod <- lm(Salary~.,data=donnees)
theta <- coef(mod)
theta

## ----teacher=correct--------------------------------------
acp.main <- eigen(t(Xcr)%*%Xcr)
U <- acp.main$vectors
CC <- Xcr%*%(-U[,1:5])
D <- cbind.data.frame(CC,Salary=Hitters$Salary)
modS <- lm(Salary~.,data=D)
coefS <- modS$coefficients
coef(modS)

## ----teacher=correct--------------------------------------
U <- acp.hit$svd$V
V <- t(U)
beta0.cr <- mean(Hitters$Salary)
beta.cr <- as.vector(theta[2:6])%*%V
beta.cr

## ----teacher=correct--------------------------------------
beta0 <- beta0.cr-sum(beta.cr*Xbar/stdX)
beta <- beta.cr/stdX
beta0
beta

## ----teacher=correct--------------------------------------
library(pls)
pcr.fit <- pcr(Salary~.,data=Hitters,scale=TRUE,ncomp=19)
coefficients(pcr.fit,ncomp=5)

## ---------------------------------------------------------
df.new <- Hitters[c(1,100,80),]

## ----teacher=correct--------------------------------------
predict(pcr.fit,newdata=df.new,ncomp=5)

## ----echo=cor,eval=cor------------------------------------
t(as.matrix(coefficients(pcr.fit,ncomp=5))) %*% 
  t(as.matrix(Xcr[c(1,100,80),]))+mean(Hitters$Salary)
#ou
beta0.cr+beta.cr%*%t(as.matrix(Xcr[c(1,100,80),]))

## ----teacher=correct--------------------------------------
beta0+beta %*% t(as.matrix(X[c(1,100,80),]))

## ----echo=FALSE,eval=FALSE--------------------------------
#  mod.lin <- lm(Salary~.,data=Hitters)
#  predict(mod.lin,newdata = df.new)
#  predict(pcr.fit,newdata=df.new,ncomp=19)
#  

## ----teacher=correct--------------------------------------
Y <- as.vector(Hitters$Salary)
w1 <- t(Xcr)%*%Y
w1
Z1 <- Xcr%*%w1

## ----teacher=correct--------------------------------------
df <- data.frame(Z1,Y)
mod1 <- lm(Y~Z1-1,data=df)
alpha1 <- coef(mod1)
alpha1

## ----teacher=correct--------------------------------------
alpha1*w1

## ----teacher=correct--------------------------------------
pls.fit <- plsr(Salary~.,data=Hitters,scale=TRUE)
coefficients(pls.fit,ncomp = 1)

## ----teacher=correct--------------------------------------
set.seed(1234)
perm <- sample(nrow(Hitters))
dapp <- Hitters[perm[1:200],]
dtest <- Hitters[perm[201:nrow(Hitters)],]

## ----teacher=correct--------------------------------------
choix.pcr <- pcr(Salary~.,data=dapp,validation="CV")
ncomp.pcr <- which.min(choix.pcr$validation$PRESS)
ncomp.pcr

## ----teacher=correct--------------------------------------
choix.pls <- plsr(Salary~.,data=dapp,validation="CV")
ncomp.pls <- which.min(choix.pls$validation$PRESS)
ncomp.pls

## ----teacher=correct--------------------------------------
mod.lin <- lm(Salary~.,data=dapp)

## ----teacher=correct--------------------------------------
prev <- data.frame(
  lin=predict(mod.lin,newdata=dtest),
  pcr=as.vector(predict(choix.pcr,newdata = dtest,ncomp=ncomp.pcr)),
  pls=as.vector(predict(choix.pls,newdata = dtest,ncomp=ncomp.pls)),
  obs=dtest$Salary
)

## ----teacher=correct--------------------------------------
prev %>% summarize_at(1:3,~(mean((.-obs)^2))) %>% sqrt()

## ----teacher=correct--------------------------------------
set.seed(1234)
bloc <- sample(1:10,nrow(Hitters),replace=TRUE)
table(bloc)

## ----teacher=correct--------------------------------------
set.seed(4321)
prev <- data.frame(matrix(0,nrow=nrow(Hitters),ncol=3))
names(prev) <- c("lin","PCR","PLS")
for (k in 1:10){
#  print(k)
  ind.test <- bloc==k
  dapp <- Hitters[!ind.test,]
  dtest <- Hitters[ind.test,]
  choix.pcr <- pcr(Salary~.,data=dapp,validation="CV")
  ncomp.pcr <- which.min(choix.pcr$validation$PRESS)
  choix.pls <- plsr(Salary~.,data=dapp,validation="CV")
  ncomp.pls <- which.min(choix.pls$validation$PRESS)
  mod.lin <- lm(Salary~.,data=dapp)
  prev[ind.test,] <- data.frame(
lin=predict(mod.lin,newdata=dtest),
PCR=as.vector(predict(choix.pcr,newdata = dtest,ncomp=ncomp.pcr)),
PLS=as.vector(predict(choix.pls,newdata = dtest,ncomp=ncomp.pls)))
}

## ----teacher=correct--------------------------------------
prev %>% mutate(obs=Hitters$Salary) %>% summarize_at(1:3,~(mean((.-obs)^2))) %>% sqrt()

## ----teacher=correct--------------------------------------
var(Hitters$Salary) %>% sqrt()

## ----teacher=correct--------------------------------------
prev1 %>% mutate(obs=Hitters$Salary) %>% summarize_at(1:3,~(mean((.-obs)^2))) %>% sqrt()

## ---------------------------------------------------------
ozone <- read.table("data/ozone.txt")
head(ozone)

## ---------------------------------------------------------
ozone.X <- model.matrix(maxO3~.,data=ozone)[,-1]
ozone.Y <- ozone$maxO3

## ----teacher=correct--------------------------------------
library(glmnet)
mod.R <- glmnet(ozone.X,ozone.Y,alpha=0)
mod.L <- glmnet(ozone.X,ozone.Y,alpha=1)

## ----teacher=correct--------------------------------------
mod.R$lambda %>% head()

## ----echo=cor,eval=cor------------------------------------
mod.R$beta[,1]

## ----teacher=correct--------------------------------------
plot(mod.R,label=TRUE)
plot(mod.L,label=TRUE)
plot(mod.R,xvar="lambda",label=TRUE)
plot(mod.L,xvar="lambda",label=TRUE)

## ----teacher=correct--------------------------------------
ridgeCV <- cv.glmnet(ozone.X,ozone.Y,alpha=0)
plot(ridgeCV)

## ----teacher=correct--------------------------------------
ridgeCV$lambda.min
ridgeCV$lambda.1se

## ----teacher=correct--------------------------------------
lassoCV <- cv.glmnet(ozone.X,ozone.Y,alpha=1)
plot(lassoCV)

## ----teacher=correct--------------------------------------
predict(ridgeCV,newx = ozone.X[50:51,],s="lambda.min")
predict(ridgeCV,newx = ozone.X[50:51,],s="lambda.1se")

## ----teacher=correct--------------------------------------
predict(lassoCV,newx = ozone.X[50:51,],s="lambda.min")
predict(lassoCV,newx = ozone.X[50:51,],s="lambda.1se")

## ---------------------------------------------------------
ozone1 <- read.table("data/ozone_complet.txt",sep=";") %>% na.omit()
ozone1.X <- model.matrix(maxO3~.,data=ozone1)[,-1]
ozone1.Y <- ozone1$maxO3

## ----teacher=correct--------------------------------------
cv.ridge.lasso <- function(data,form){
  set.seed(1234)
  data.X <- model.matrix(form,data=data)[,-1]
  data.Y <- data$maxO3
  blocs <- caret::createFolds(1:nrow(data),k=10)
  prev <- matrix(0,ncol=3,nrow=nrow(data)) %>% as.data.frame()
  names(prev) <- c("lin","ridge","lasso")
  for (k in 1:10){
app <- data[-blocs[[k]],]
test <- data[blocs[[k]],]
app.X <- data.X[-blocs[[k]],]
app.Y <- data.Y[-blocs[[k]]]
test.X <- data.X[blocs[[k]],]
test.Y <- data.Y[blocs[[k]]]
ridge <- cv.glmnet(app.X,app.Y,alpha=0)
lasso <- cv.glmnet(app.X,app.Y,alpha=1)
lin <- lm(form,data=app)
prev[blocs[[k]],] <- tibble(lin=predict(lin,newdata=test),
           ridge=as.vector(predict(ridge,newx=test.X)),
           lasso=as.vector(predict(lasso,newx=test.X)))
  }
  err <- prev %>% mutate(obs=data$maxO3) %>% summarise_at(1:3,~mean((obs-.)^2))
  return(err)
}

## ----echo=cor,eval=cor------------------------------------
ozone2.X <- model.matrix(maxO3~.^2,data=ozone1)[,-1]
dim(ozone2.X)

## ---------------------------------------------------------
signal <- read_csv("data/signal.csv")
ggplot(signal)+aes(x=x,y=y)+geom_line()

## ---------------------------------------------------------
donnees <- read_csv("data/ech_signal.csv")
ggplot(signal)+aes(x=x,y=y)+geom_line()+
  geom_point(data=donnees,aes(x=X,y=Y))

## ----teacher=correct--------------------------------------
mat.dict <- function(K,x){
	res <- matrix(0,nrow=length(x),ncol=2*K) %>% as_tibble()
	for (j in 1:K){
	  res[,2*j-1] <- cos(2*j*pi*x)
		res[,2*j] <- sin(2*j*pi*x)
	}
	return(res)
}

## ----teacher=correct--------------------------------------
D25 <- mat.dict(25,donnees$X) %>% mutate(Y=donnees$Y)
mod.lin <- lm(Y~.,data=D25)

## ----teacher=correct--------------------------------------
S25 <- mat.dict(25,signal$x)
prev.MCO <- predict(mod.lin,newdata = S25)
signal1 <- signal %>% mutate(MCO=prev.MCO) %>% rename(signal=y)
signal2 <- signal1 %>% pivot_longer(-x,names_to="meth",values_to="y")
ggplot(signal2)+aes(x=x,y=y)+geom_line(aes(color=meth))+
  scale_y_continuous(limits = c(-2,2))+geom_point(data=donnees,aes(x=X,y=Y))

## ----teacher=correct--------------------------------------
X.25 <- model.matrix(Y~.,data=D25)[,-1]
lasso1 <- glmnet(X.25,D25$Y,alpha=1)
plot(lasso1)

## ----teacher=correct--------------------------------------
lasso.cv <- cv.glmnet(X.25,D25$Y,alpha=1)
plot(lasso.cv)

## ----teacher=correct--------------------------------------
prev.lasso <- as.vector(predict(lasso.cv,newx=as.matrix(S25)))
signal1$lasso <- prev.lasso
signal2 <- signal1 %>% pivot_longer(-x,names_to="meth",values_to="y")
ggplot(signal2)+aes(x=x,y=y)+geom_line(aes(color=meth))+
  scale_y_continuous(limits = c(-2,2))+geom_point(data=donnees,aes(x=X,y=Y))

## ----teacher=correct--------------------------------------
v.sel <- which(coef(lasso.cv)!=0)
v.sel

## ----teacher=correct--------------------------------------
pcr.fit <- pcr(Y~.,data=D25,validation="CV")
ncomp.pcr <- which.min(pcr.fit$validation$PRESS)
ncomp.pcr
prev.pcr <- predict(pcr.fit,newdata=S25,ncomp=ncomp.pcr)

## ----teacher=correct--------------------------------------
pls.fit <- plsr(Y~.,data=D25,validation="CV")
ncomp.pls <- which.min(pls.fit$validation$PRESS)
ncomp.pls
prev.pls <- predict(pls.fit,newdata=S25,ncomp=ncomp.pls)

## ----teacher=correct--------------------------------------
signal1$pcr <- prev.pcr
signal1$pls <- prev.pls
signal2 <- signal1 %>% pivot_longer(-x,names_to="meth",values_to="y")
ggplot(signal2)+aes(x=x,y=y)+geom_line(aes(color=meth))+
  scale_y_continuous(limits = c(-2,2))+geom_point(data=donnees,aes(x=X,y=Y))


## ----teacher=correct--------------------------------------
signal1 %>% summarise_at(-(1:2),~mean((.-signal)^2)) %>%
  sort() %>% round(3)

## ---------------------------------------------------------
ad.data <- read.table("data/ad_data.txt",header=FALSE,sep=",",dec=".",na.strings = "?",strip.white = TRUE)
names(ad.data)[ncol(ad.data)] <- "Y"
ad.data$Y <- as.factor(ad.data$Y)

## ---------------------------------------------------------
summary(ad.data$Y)

## ----teacher=correct--------------------------------------
sum(is.na(ad.data))
var.na <- apply(is.na(ad.data),2,any)
names(ad.data)[var.na]
ind.na <- apply(is.na(ad.data),1,any)
sum(ind.na)

## ----teacher=correct--------------------------------------
ad.data1 <- ad.data[,var.na==FALSE]
dim(ad.data1)
sum(is.na(ad.data1))

## ----teacher=correct--------------------------------------
X.ad <- model.matrix(Y~.,data=ad.data1)[,-1]
Y.ad <- ad.data1$Y

## ----cv-ad,echo=correct,eval=FALSE------------------------
#  set.seed(5678)
#  blocs <- caret::createFolds(1:nrow(ad.data1),k=10)
#  score <- matrix(0,ncol=3,nrow=nrow(ad.data1)) %>% as.data.frame()
#  names(score) <- c("MV","ridge","lasso")
#  for (k in 1:10){
#    print(k)
#    app <- ad.data1[-blocs[[k]],]
#    test <- ad.data1[blocs[[k]],]
#    app.X <- X.ad[-blocs[[k]],]
#    app.Y <- Y.ad[-blocs[[k]]]
#    test.X <- X.ad[blocs[[k]],]
#    test.Y <- Y.ad[blocs[[k]]]
#    ridge <- cv.glmnet(app.X,app.Y,family="binomial",alpha=0)
#    lasso <- cv.glmnet(app.X,app.Y,family="binomial",alpha=1)
#    MV <- glm(Y~.,data=app,family="binomial")
#    score[blocs[[k]],] <- tibble(MV=predict(MV,newdata=test,type="response"),
#               ridge=as.vector(predict(ridge,newx=test.X,type="response")),
#               lasso=as.vector(predict(lasso,newx=test.X,type="response")))
#  }

## ----echo=FALSE,eval=FALSE--------------------------------
#  #Pour compiler plus vite
#  write_csv(score,path="score_cv_ad.csv")

## ----echo=FALSE,eval=correct------------------------------
score <- read_csv("score_cv_ad.csv")

## ----teacher=correct--------------------------------------
score1 <- score %>% 
  mutate(obs=fct_recode(ad.data1$Y,"0"="ad.","1"="nonad.")) %>%
  pivot_longer(-obs,names_to="Methode",values_to="score")
ggplot(score1)+aes(m=score,d=as.numeric(obs),color=Methode)+plotROC::geom_roc()

## ----teacher=correct--------------------------------------
score1 %>% group_by(Methode) %>% 
      summarize(AUC=round(as.numeric(pROC::auc(obs,score)),3)) %>% 
      arrange(desc(AUC)) 


## ----teacher=correct--------------------------------------
score1 %>% mutate(prev=round(score),err=prev!=obs) %>% 
  group_by(Methode) %>% summarize(Err_classif=round(mean(err),3)) %>%
  arrange(Err_classif) 

## ----teacher=correct--------------------------------------
set.seed(1234)
n <- 300
X1<-runif(n)
X2<-runif(n)
bruit<-rnorm(n)
Y<-1+3*X1+5*X2+bruit
donnees<-data.frame(Y,X1,X2)

## ----teacher=correct--------------------------------------
pseudo_back <- function(df,eps=0.00001){
  mat.X <- model.matrix(Y~.,data=df)
  beta_i <- rep(0,ncol(mat.X))
  beta <- rep(1,ncol(mat.X))
  while (min(abs(beta_i-beta))>eps){
beta_i <- beta
for (k in 1:ncol(mat.X)){
  Yk <- Y-mat.X[,-k]%*%(beta[-k])
  dfk <- data.frame(Yk=Yk,Xk=mat.X[,k])
  beta[k]<-coef(lm(Yk~Xk-1,data=dfk))
}
  }
  return(beta)
}

## ----teacher=correct--------------------------------------
pseudo_back(donnees)

## ----teacher=correct--------------------------------------
lm(Y~.,data=donnees)

## ---------------------------------------------------------
n <- 1000
set.seed(1465)
X1 <- 2*runif(n)
X2 <- 2*runif(n)
bruit <- rnorm(n)
Y <- 2*X1+sin(8*pi*X2)+bruit
donnees<-data.frame(Y,X1,X2)

## ----teacher=correct--------------------------------------
library(gam)
model1 <- gam(Y~s(X1,df=1)+s(X2,df=24.579)-1,data=donnees)
plot(model1)

## ----teacher=correct--------------------------------------
model2 <- gam(Y~s(X1,df=0.001)+s(X2,df=0.001)-1,data=donnees)
plot(model2)

## ----teacher=correct--------------------------------------
model2 <- gam(Y~s(X1,df=100)+s(X2,df=100)-1,data=donnees)
plot(model2)

## ----teacher=correct--------------------------------------
model4 <- gam(Y~lo(X1,span=3)+lo(X2,span=0.15,degree=2)-1,data=donnees)
plot(model4)

## ----teacher=correct--------------------------------------
model5 <- gam(Y~lo(X1,span=5)+lo(X2,span=5,degree=2)-1,data=donnees)
plot(model5)

## ----teacher=correct--------------------------------------
model6 <- gam(Y~lo(X1,span=0.01)+lo(X2,span=0.01,degree=2)-1,data=donnees)
plot(model6)

## ----teacher=correct--------------------------------------
mod.mgcv <- mgcv::gam(Y~s(X1)+s(X2),data=donnees)
plot(mod.mgcv)

## ----teacher=correct--------------------------------------
panne <- read.table("data/panne.txt",header=TRUE)
mod1 <- glm(etat~age,data=panne,family=binomial)
summary(mod1)

## ----teacher=correct--------------------------------------
mod.panne <- mgcv::gam(etat~s(age),data=panne)
plot(mod.panne)

## ----teacher=correct--------------------------------------
mod2 <- glm(etat~age+I(age^2),data=panne,family=binomial)
summary(mod2)

