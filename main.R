
# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
library(compositions); library(DirichletReg); library(glmmTMB)
library(tikzDevice); library(xtable); library(kableExtra)

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


# Data --------------------------------------------------------------------
load("data_casestudy.RData")
str(datax)
summary(datax)

## Split dataset into two parts: 
# (A) IRTree-based quantification, (B) Compositional data analysis
set.seed(080224)
iid = sample(x = 1:NROW(datax),size = NROW(datax)*0.5,replace = FALSE) 
iidC = setdiff(1:NROW(datax),iid)
datax_A = datax[iid,]; datax_B = datax[iidC,]


# IRTree-based quantification ---------------------------------------------
# For info about the 'irtrees' library, see: https://cran.r-project.org/web/packages/irtrees/vignettes/pretree_vignette.html

## Using a quite common decision tree
M=5;N=M-1; I=NROW(datax_A); J=6
Tm = matrix(c(1,0,0,NA,
              1,0,1,NA,
              0,NA,NA,NA,
              1,1,NA,0,
              1,1,NA,1),M,N,byrow = TRUE)
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N)
print(Tm)

## Estimating IRTree parameters for the item component only (dataset A)
Y = irtrees::dendrify(as.matrix(datax_A[,2:7]),Tm)
mod_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Y)
Alpha = matrix(summary(mod_glm)$coefficients$cond[,1],J,N) #easiness/difficulty to answer a given items (across decision nodes)
Sigma_eta = summary(mod_glm)$varcor$cond$person 
R_eta = attr(x = Sigma_eta,which = 'correlation') 
Sigma_eta = matrix(as.numeric(Sigma_eta),N,N) 
Etas = ranef(mod_glm)$cond$person #raters ability

## Table 1
#(a)
Xtab = summary(mod_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:6,],Xtab[7:12,],Xtab[13:18,],Xtab[19:24,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),N)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study: Estimates of item parameters"
attributes(Xtab_tex)$label = "tab1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})
#(b)
Xtab = R_eta; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(Sigma_eta)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 1: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs1_3"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Computing response probabilities (prediction on dataset B) 
Y = irtrees::dendrify(as.matrix(datax_B[,2:7]),Tm)
mod_glm2 = update(mod_glm, . ~ ., data = Y)
Eta = apply(ranef(mod_glm2)$cond$person,2,as.numeric)

## Computing P(Y in {1,2,3,4,5})|Alpha,Eta (raters decision probabilities)
I = NROW(datax_B)
PYy = array(data = NA,dim = c(I,J,M))
Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}
for(i in 1:I){
  for(j in 1:J){
    PYy[i,j,] = mapply(function(m){prod((exp((Eta[i,]+Alpha[j,])*Tm[m,])/(1+exp(Eta[i,]+Alpha[j,])))^Dm[m,])},1:M)
  }
}

## Figure 2
cols = c("#133955","#81BE83","#ECEBBD","#D79FC7","#FF9B49")
tikzDevice::tikz(file='fig2.tex',width=7,height=3,sanitize = TRUE)
par(mar=c(5,8,5,2)+0.1,mfrow=c(1,3)) # Doubles left margin.
i=102; barplot(t(PYy[i,,]),names.arg = paste0("item ",1:6),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = cols,cex.names = 1.35); title(paste0("Rater ",i), line = 1,adj=0,cex=1.35)
i=203; barplot(t(PYy[i,,]),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = cols,cex.names = 1.35); title(paste0("Rater ",i), line = 1,adj=0,cex=1.35)
i=304; barplot(t(PYy[i,,]),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = cols,cex.names = 1.35); title(paste0("Rater ",i), line = 1,adj=0,cex=1.35)
add_legend("bottom",legend = 1:5,fill = cols,border = FALSE,bty = "n",cex=1.5,ncol=5)
dev.off()



# CoDa --------------------------------------------------------------------
# Note: in what follows, we ignore the information due to the order of the different classes/responses

## Preparing dataset (it may contains variables currently not used for the analysis)
datax_coda = data.frame(apply(PYy,c(1,3),mean),datax_B) #averaging across items given a rater
# Note: in this case, the arithm mean equals the Aitchnson-based mean (see the function mean.acomp())

datax_coda$GENDER = factor(datax_coda$GENDER,labels = c("F","M"))
str(datax_coda)
summary(datax_coda)

X = data.frame(y=as.numeric(as.matrix(datax_coda[,1:5])),
          resp=rep(1:5,each=NROW(datax_coda)),     
          gender=rep(datax_coda$GENDER,5),
          depression=rep(datax_coda$DEPRESSION,5));
      

## Figure 3
tikzDevice::tikz(file='fig3.tex',width=7,height=3,sanitize = TRUE)
par(mfrow=c(1,3))
plot(0,0,bty="n",col=0,ylim=c(0,1),xlim=c(0.9,max(X$depression)+0.5),xlab="",ylab="")
for(j in 1:5){points(X$depression[X$resp==j&X$gender=="F"],X$y[X$resp==j&X$gender=="F"],col=cols[j],pch=19,cex.lab=2.35,cex.axis=1.25)}
title("Depression (Female)", line = 1,adj=0,cex=1.25)
plot(0,0,bty="n",col=0,ylim=c(0,1),xlim=c(0.9,max(X$depression)+0.5),xlab="",ylab="")
for(j in 1:5){points(X$depression[X$resp==j&X$gender=="M"],X$y[X$resp==j&X$gender=="M"],col=cols[j],pch=19,cex.lab=1.35,cex.axis=1.25)}
title("Depression (Male)", line = 1,adj=0,cex=1.25)
boxplot(datax_coda[datax_coda$GENDER=="F",1:5],frame=FALSE,col=cols,names = rep("F",5),notch = TRUE,at = seq(1,10,length=5),xlim=c(0,12),border = "#132C2D")
boxplot(datax_coda[datax_coda$GENDER=="M",1:5],frame=FALSE,col=cols,names = rep("M",5),add = TRUE,notch = TRUE,at = seq(2,11,length=5),border = "#132C2D")
title("Gender", line = 1,adj=0,cex=1.25)
add_legend("bottom",legend = 1:5,fill = cols,border = FALSE,bty = "n",cex=1.5,ncol=5)
dev.off()


## Dirichlet regression for compositional responses
datax_coda$Y = DR_data(Y = datax_coda[,1:5])

mod_full = DirichReg(formula = Y~DEPRESSION+GENDER,data = datax_coda,model="alternative")
drop1(mod_full,sort = TRUE)

mod0 = DirichReg(formula = Y~DEPRESSION+GENDER|1,model="alternative",data = datax_coda,base = 1)
mod1 = DirichReg(formula = Y~DEPRESSION+GENDER+DEPRESSION:GENDER|1,model="alternative",data = datax_coda,base = 1)
anova(mod0,mod1)
summary(mod1)

## Table 2
out = summary(mod1)
Xtab = cbind(out$coef.mat[1:4,1:3],out$coef.mat[5:8,1:3],out$coef.mat[9:12,1:3],out$coef.mat[13:16,1:3])
Xtab = rbind(Xtab,c(out$coef.mat[17,1:3],rep(NA,9)))
rownames(Xtab) = c("$\\beta_{0}$","$\\beta_{\\text{depres}}$","$\\beta_{\\text{gender:M}}$","$\\beta_{\\text{depres x gender}}$","\\phi")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$","$z$"),(M-1))
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study: "
attributes(Xtab_tex)$label = "tab2"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


## Figure 4
tikzDevice::tikz(file='fig4.tex',width=7,height=3,sanitize = TRUE)
Xnew = data.frame(DEPRESSION=c(seq(min(datax_coda$DEPRESSION[datax_coda$GENDER=="F"]), max(datax_coda$DEPRESSION[datax_coda$GENDER=="F"]),length.out=100),
                               seq(min(datax_coda$DEPRESSION[datax_coda$GENDER=="M"]), max(datax_coda$DEPRESSION[datax_coda$GENDER=="M"]),length.out=100)),
                  GENDER=rep(c("F","M"),each=100)); Ynew = predict.DirichletRegModel(object = mod1,newdata = Xnew)
par(mfrow=c(1,5))
for(j in 1:5){
  plot(0,0,bty="n",col=0,ylim=c(0,1),xlim=c(0.9,max(X$depression)+0.5),xlab="",ylab="")
  points(X$depression[X$resp==j&X$gender=="F"],X$y[X$resp==j&X$gender=="F"],col="#708090",pch=20,cex.lab=2.35,cex.axis=1.25)
  points(X$depression[X$resp==j&X$gender=="M"],X$y[X$resp==j&X$gender=="M"],col="#93C572",pch=20,cex.lab=2.35,cex.axis=1.25)
  lines(Xnew[Xnew$GENDER=="F",1],Ynew[Xnew$GENDER=="F",j],pch=20,col="#708090",lwd=4)
  lines(Xnew[Xnew$GENDER=="M",1],Ynew[Xnew$GENDER=="M",j],pch=20,col="#93C572",lwd=4)
  title(paste0("(Response) y=",j), line = 1,adj=0,cex=1.25)
}
add_legend("bottom",legend = c("Female","Male"),fill = c("#708090","#93C572"),border = FALSE,bty = "n",cex=1.5,ncol=2)
dev.off()





