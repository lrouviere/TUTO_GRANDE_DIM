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

## ----eval=cor,echo=TRUE-----------------------------------
#  k.opt <- sel.k(seq(1,50,by=5),df$Xapp,df$Yapp)
#  prev <- knn.reg(train=df$Xapp,y=df$Yapp,test=df$Xtest,k=k.opt)$pred
#  mean((prev-df$Ytest)^2)

## ---------------------------------------------------------
DIM <- c(1,5,10,50)
K_cand <- seq(1,50,by=5)

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

## ---------------------------------------------------------
ozone <- read.table("data/ozone.txt")
head(ozone)

## ---------------------------------------------------------
library(leaps)
mod.sel <- regsubsets(maxO3~.,data=ozone,nvmax=14)
summary(mod.sel)

## ----chunk-bestglm,echo=cor,eval=cor,cache=TRUE,indent='        '----
#  ozone1 <- ozone %>% mutate(vent=as.factor(vent),pluie=as.factor(pluie)) %>%
#    select(-maxO3,everything())
#  library(bestglm)
#  model.bglm <- bestglm(ozone1,IC="BIC")
#  model.bglm$BestModel %>% summary()

## ---------------------------------------------------------
library(ISLR)
Hitters <- na.omit(Hitters)

## ---------------------------------------------------------
df.new <- Hitters[c(1,100,80),]

## ----echo=cor,eval=cor------------------------------------
#  t(as.matrix(coefficients(pcr.fit,ncomp=5))) %*%
#    t(as.matrix(Xcr[c(1,100,80),]))+mean(Hitters$Salary)
#  #ou
#  beta0.cr+beta.cr%*%t(as.matrix(Xcr[c(1,100,80),]))

## ----echo=FALSE,eval=FALSE--------------------------------
#  mod.lin <- lm(Salary~.,data=Hitters)
#  predict(mod.lin,newdata = df.new)
#  predict(pcr.fit,newdata=df.new,ncomp=19)
#  

## ---------------------------------------------------------
ozone <- read.table("data/ozone.txt")
head(ozone)

## ---------------------------------------------------------
ozone.X <- model.matrix(maxO3~.,data=ozone)[,-1]
ozone.Y <- ozone$maxO3

## ----echo=cor,eval=cor------------------------------------
#  mod.R$beta[,1]

## ---------------------------------------------------------
ozone1 <- read.table("data/ozone_complet.txt",sep=";") %>% na.omit()
ozone1.X <- model.matrix(maxO3~.,data=ozone1)[,-1]
ozone1.Y <- ozone1$maxO3

## ----echo=cor,eval=cor------------------------------------
#  ozone2.X <- model.matrix(maxO3~.^2,data=ozone1)[,-1]
#  dim(ozone2.X)

## ---------------------------------------------------------
signal <- read_csv("data/signal.csv")
ggplot(signal)+aes(x=x,y=y)+geom_line()

## ---------------------------------------------------------
donnees <- read_csv("data/ech_signal.csv")
ggplot(signal)+aes(x=x,y=y)+geom_line()+
  geom_point(data=donnees,aes(x=X,y=Y))

## ---------------------------------------------------------
ad.data <- read.table("data/ad_data.txt",header=FALSE,sep=",",dec=".",na.strings = "?",strip.white = TRUE)
names(ad.data)[ncol(ad.data)] <- "Y"
ad.data$Y <- as.factor(ad.data$Y)

## ---------------------------------------------------------
summary(ad.data$Y)

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
#  score <- read_csv("score_cv_ad.csv")

## ---------------------------------------------------------
n <- 1000
set.seed(1465)
X1 <- 2*runif(n)
X2 <- 2*runif(n)
bruit <- rnorm(n)
Y <- 2*X1+sin(8*pi*X2)+bruit
donnees<-data.frame(Y,X1,X2)

