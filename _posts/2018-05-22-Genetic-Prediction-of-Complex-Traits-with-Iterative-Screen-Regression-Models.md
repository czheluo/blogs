---
layout: single
author_profile: false
comments: true
toc: true
enable_mathjax: true
output: html_document
title : "Genetic Prediction of Complex Traits with Iterative Screen Regression Models"
output: html_document
tags: [GS, ISR, Bayes alphabetic model, ]
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

## Iterative Screen Regression model (ISR)

As mentioned [before](http://mengluocv.me/blog/2018/05/12/An-Efficient-Iterative-Screen-Regression-model-For-Genome-wide-Association-Study-in-Structured-Populations/) that our method, here, was used to genomic prediction.
<div style="align:center;width:40;"><img src="{{ "/images/Blog/GS/ISRGS.jpg" | prepend: site.baseurl }}"></div>
<p style="text-align: center;"> Schematic overview of model-based iterative screen regression for GS </p>
## Simulation

I used genotypes from an existing cattle GWAS [data set](http://www.g3journal.org/content/5/4/615) with
5024 individuals and 42,551 SNPs and simulated phenotypes. And consider the TWO following different simulation scenarios to a range of possible genetic architectures.
Scenario I, where only set 500 SNPs are causal. Specifically, I randomly selected 50 group-one SNPs, 150 group-two SNPs, 300 group-three SNPs. We simulated SNP effect sizes all from a standard normal distribution. But scaled their effects in each group separately so that the proportion of genetic variance explained by the four groups are 0.15, 0.25, and 0.60, respectively. We set the total proportion of phenotypic variance (PVE; i.e., SNP heritability) to be either 0.2, 0.5, or 0.8, representing low, moderate, and high heritability, respectively.

Scenario II satisfies the BayesR modeling assumption, where a small proportion of SNPs are causal (500 SNPs). These causal SNPs come from three effect-size groups. The simulations were similar to scenario I. Here, the proportion of PVE by the three groups are 0.15, 0.25, and 0.6, respectively. Again, we set the total PVE to be either 0.2, 0.5, or 0.8. This simulation scenario consists of one simulation setting for each PVE.

In each setting, I performed 20 simulation replicates. In each replicate, I randomly split the data into a training data with 80% individuals and a test data with the remaining 20% individuals. I then fitted different methods on the training data and evaluated their prediction performance on the test data(i.e., Monte Carlo cross validation).

## Performance (Prediction accuracy)

Performance was measured by R (Prediction accuracy is the Pearson’s correlation) and the R difference with respect to ISR, where a negative value (i.e., values below the red horizontal line) indicates worse performance than ISR. The sample R differences are obtained from 20 replicates in each scenario. Also performance was measured by MSE difference with respect to ISR, where a positive value (i.e. values above the red horizontal line) indicates worse performance than ISR. The sample MSE differences are obtained from 20 replicates in each scenario.Methods for comparison include ISR(blue), DPR.MCMC (cyan), BayesB (black), BayesA (green), Bayes LASSO (red), DPR.VB (red), BSLMM (yellow), rrBLUP (magenta), and BayesC (grey).

```matlab
% ISR
% ten-fold cross validation
clear,clc
warning('off','all')
 load PHS_all.mat
  %load qtn100_wheat_0.2.mat
 Y=y;X=x;
mn=20;
nfolds=10;foldid=[];[n,p]=size(x);rinf=zeros(nfolds,mn);
rref=zeros(nfolds,mn);test=zeros(n,mn);
tic
lmn=2;
for ny=1:20
   N = size(x,1);
if (isempty(foldid))
    foldid = randsample([repmat(1:nfolds,1,floor(N/nfolds)) 1:mod(N,nfolds)],N);
else
    nfolds = max(foldid);
end
for i=1:nfolds
    which=foldid==i;
    %if verbous, disp(['Fitting fold # ' num2str(i) ' of ' num2str(nfolds)]);end
    x=X(~which,:);y=Y(~which,lmn);
    [b,effect,plm,seblm,ft,lmbx,lm] =ISR(x,y,2,2,0.01,3,25);
    %cvfit =glmnet(x(~which,:), y(~which),family, options);
    %predmat(which,:) = glmnetPredict(cvfit, type,x(which,:),options.lambda);
    x0=[ones(size(X(which,:),1),1) X(which,lm)];
    x1=[ones(size(X(~which,:),1),1) X(~which,lm)];
    %eff=zeros(size(x,2)+1,1);eff(1)=b(1); eff(2:end)=effect;
    predmat1=x0*b;
    predmat2=x1*b;
    rinf(i,ny)=corr(Y(which,lmn),predmat1);
    rref(i,ny)=corr(Y(~which,lmn),predmat2);
    test(which,ny)=x0*b;
end
save PHS_cattle.mat test rref rinf
end
toc
```

```r
library(R.matlab)
library(glmnet)
library(rrBLUP)
library(BGLR)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
#####simulation####
###BayesB
mat<-readMat('sim05_cattle_multple_qtn150.mat')
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20
acc1<-matrix(0,mn,1)
#ans2<-NULL
acc2<-matrix(0,mn,1)
ytest<-matrix(0,1005,mn)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]

  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBB=BGLR(y=ytr,ETA=list( list(X=T,model='BayesB')),
              nIter=nIter,burnIn=burnIn,saveAt='bb')
    acc1[ij-1]<-cor(t%*%fmBB$ETA[[1]]$b+fmBB$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBB$ETA[[1]]$b+fmBB$mu,ytr)
    ytest[,ij-1]<-t%*%fmBB$ETA[[1]]$b+fmBB$mu
    writeMat(paste('BB/cattle/sim05/y_',ij-1,'_cattle_sim05qtn150_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
  }
  writeMat(paste('BB/cattle/sim05/cattle_sim05qtn150_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1

###BayesA
mat<-readMat('sim05_cattle_multple_qtn150.mat')
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20
acc1<-matrix(0,mn,1)
#ans2<-NULL
acc2<-matrix(0,mn,1)
ytest<-matrix(0,1005,mn)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]

  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBA=BGLR(y=ytr,ETA=list( list(X=T,model='BayesA')),
              nIter=nIter,burnIn=burnIn,saveAt='ba')
    acc1[ij-1]<-cor(t%*%fmBA$ETA[[1]]$b+fmBA$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBA$ETA[[1]]$b+fmBA$mu,ytr)
    ytest[,ij-1]<-t%*%fmBA$ETA[[1]]$b+fmBA$mu
    writeMat(paste('BA/cattle/sim05/y_',ij-1,'_cattle_sim05qtn150_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
  }
  writeMat(paste('BA/cattle/sim05/cattle_sim05qtn150_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1

###BayesC
mat<-readMat('sim05_cattle_multple_qtn150.mat')
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20
acc1<-matrix(0,mn,1)
#ans2<-NULL
acc2<-matrix(0,mn,1)
ytest<-matrix(0,1005,mn)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]

  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBC=BGLR(y=ytr,ETA=list( list(X=T,model='BayesC')),
              nIter=nIter,burnIn=burnIn,saveAt='bc')
    acc1[ij-1]<-cor(t%*%fmBC$ETA[[1]]$b+fmBC$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBC$ETA[[1]]$b+fmBC$mu,ytr)
    ytest[,ij-1]<-t%*%fmBC$ETA[[1]]$b+fmBC$mu
    writeMat(paste('BC/cattle/sim05/y_',ij-1,'_cattle_sim05qtn150_BC.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
  }
  writeMat(paste('BC/cattle/sim05/cattle_sim05qtn150_BC.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1

###BayesLASSO
mat<-readMat('sim05_cattle_multple_qtn150.mat')
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20
acc1<-matrix(0,mn,1)
#ans2<-NULL
acc2<-matrix(0,mn,1)
ytest<-matrix(0,1005,mn)
nIter=1500; burnIn=1000
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]

  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBL=BGLR(y=ytr,ETA=list( list(X=T,model='BL')),
              nIter=nIter,burnIn=burnIn,saveAt='bl')
    acc1[ij-1]<-cor(t%*%fmBL$ETA[[1]]$b+fmBL$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBL$ETA[[1]]$b+fmBL$mu,ytr)
    ytest[,ij-1]<-t%*%fmBL$ETA[[1]]$b+fmBL$mu
    writeMat(paste('BL/cattle/sim05/y_',ij-1,'_cattle_sim05qtn150_BL.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
  }
  writeMat(paste('BL/cattle/sim05/cattle_sim05qtn150_BL.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1

###RRBLUP

mat<-readMat('sim05_cattle_multple_qtn150.mat')
Y<-read.csv('cattle_sim0.5_qtn150.csv',header = TRUE)
mn<-20
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]

  for (j in 1:1) {
    M<-as.matrix(mat$x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    RR<- mixed.solve(y=ytr,Z=T)
    acc1[ij-1]<-cor(t%*%RR$u,yte)
    acc2[ij-1]<-cor(T%*%RR$u,ytr)
    ytest[,ij-1]<-t%*%RR$u
    writeMat(paste('RR/cattle/sim05/y_',ij-1,'_cattle_sim05qtn150_RR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
  }
  writeMat(paste('RR/cattle/sim05/cattle_sim05qtn150_RR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1

###DPR

Y<-read.csv('PHS_pheno.csv',header = TRUE)
geno.file<-'PHS_wheat.txt.gz'
pheno.file<-'PHS_wheat.pheno.txt'
out.file<-'PHS_wheat'
mn<-20
acc1<-matrix(0,10,mn)
acc2<-matrix(0,10,mn)
ytest<-matrix(0,185,mn)
nfolds<-10
time1 <- Sys.time()
for (ij in 2:11) {
  luoy<-Y[,c(1,ij)]

  for (j in 1:mn) {

    #luoy<-luoy[sample(nrow(luoy)),]
    folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
nfolds<-10
time1 <- Sys.time()
for (i in 1:nfolds) {
  y<-matrix(NA,185,1)
  testing<- which(folds==i)
  y[-testing]<-luoy[,2][-testing]
  write.table(y,file = 'PHS_wheat.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
  system(sprintf("./DPR -g %s -p %s -dpr 2 -maf 0.00000001 -w 3000 -s 3000 -o %s",
                 geno.file,pheno.file,out.file),ignore.stdout = FALSE)
  system(sprintf("./DPR -g %s -p %s -epm output/PHS_wheat.param.txt -emu output/PHS_wheat.log.txt -predict -o %s",
                 geno.file,pheno.file,out.file),ignore.stdout = FALSE)
  pred<-read.table(file='output/PHS_wheat.prdt.txt',header = FALSE)
  acc1[i,j]<-cor(pred$V1[testing],luoy[,2][testing])
  acc2[i,j]<-cor(pred$V1[-testing],luoy[,2][-testing])
  ytest[testing,j]<-pred$V1[testing]
}
  }
  writeMat(paste('PHS/y_',ij-1,'_PHS_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytest)
}
t1 <- f(time1)
t1
###BSLMM

geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim08qtn500.pheno.txt'
out.file<-'BS_cattle_5024sim08qtn500'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.8_qtn500.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 16:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim08qtn500.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim08qtn500.param.txt -emu output/BS_cattle_5024sim08qtn500.log.txt -ebv output/BS_cattle_5024sim08qtn500.bv.txt -k output/BS_cattle_5024sim08qtn500.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim08qtn500.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim08/sim08qtn500',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim08/sim08qtn500_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 <- f(time1)
t1
```

* Scenario I
<div align="center"><img src="{{ "/images/Blog/GS/I.jpg" | prepend: site.baseurl }}"></div>
* Scenario II
<div align="center"><img src="{{ "/images/Blog/GS/II.jpg" | prepend: site.baseurl }}"></div>
## Real data applications

- Cattle data
  The cattle data consists of 5,024 samples and 42,551 SNPs after removing SNPs that have a HWE p-value < 10.4, a genotype call rate <95%, or an MAF < 0.01. For the remaining SNPs, I imputed missing genotypes with the estimated mean genotype of that SNP. I analyzed three traits: MFP, MY, and SCS. All phenotypes were quantile normalized to a standard normal distribution before analysis.
- Rice data
  The rice data consists of 1,132 samples and 464,831 SNPs after removing SNPs that have a genotype call rate <95%, or an MAF < 0.05. For the remaining SNPs, I also imputed missing genotypes with the estimated mean genotype of that SNP. Only one trait geain length(GL) was used to analysis.

<div align="center"><img src="{{ "/images/Blog/GS/cattleGL.jpg" | prepend: site.baseurl }}"></div>
* Mice data (CFW)
Another outbred CFW (Carworth Farms White) mice population that including a set of 92,734 single-nucleotide polymorphism markers which were genotyped 1,161 individuals were also used to analysis.
<div align="center"><img src="{{ "/images/Blog/GS/mice.jpg" | prepend: site.baseurl }}"></div>
```matlab
% coded for MSE BOXPLOT
clear,clc;
%load GS_02qtn100.mat
%load GS_05qtn100.mat
%load GS_02qtn150.mat
%load GS_02qtn500.mat
load GS_02qtn500_va.mat
%load cl8.mat
load cl8_1.mat
[n,p]=size(x);X=zeros(n,p);
for i=1:8
    for j=1:20
    X(j,i)=mean((x(j,i)-mean(x(:,i))).^2);
    end
end
for i=1:8
    	XX(:,i)=X(:,i)-X(:,1);
end
x=XX;x(:,1)=[];x=x(:);
n=20 ; xx=(1:7)';
r=repmat(xx,1,n)';
g=r(:)';
Y=1:0.04:1.24;
h6=boxplot(x,g,'positions',Y,'Notch','off','Symbol','.','widths',0.03); %1:0.5:4.5
set(h6,'linewidth',1)
%set(gca,'xtick',[mean(y(1:2)) mean(y(3:4)) ...
    %mean(y(5:6)) mean(y(7:8)) mean(y(9:10))])
set(gca,'xticklabel',{'SR','BayesB','BayesA','BayesC','Bayes LASSO','DPR','BSLMM','rrBLUP'},'FontName','Times New Roman',...
    'FontWeight','bold','Fontsize',14,'XTickLabelRotation',45)
   set(gca,'xticklabel',[],'FontName','Times New Roman','FontWeight','bold','Fontsize',14);
   ylabel('R difference',...
       'FontName','Times New Roman','FontWeight','bold','Fontsize',14);
   box on;grid on;
h= findobj(gca,'Tag','Box');
%cl=unifrnd(0,1,8,3);
col=zeros(8,3);
for i=1:8
    if isodd(i)
        col(i,:)=cl(i,:);
    else
        col(i,:)=cl(i,:);
    end
end
for j=1:length(h)
   % if isodd(i)
        h4(j)=patch(get(h(j),'XData'),get(h(j),'YData'),col(j,:),'EdgeColor',col(j,:),'FaceAlpha',0.8);
    %else
       %h4(j)=patch(get(h(j),'XData'),get(h(j),'YData'),col(j,:),'EdgeColor',col(j,:),'FaceAlpha',0.5);
    %end
end
title('h^2=0.2')
set(gca,'xticklabel',[],'FontName','Times New Roman','FontWeight','bold','Fontsize',14);
  % ylabel('MSE difference',...
     %  'FontName','Times New Roman','FontWeight','bold','Fontsize',14);
box on;grid on;
ylabel('MSE difference',...
       'FontName','Times New Roman','FontWeight','bold','Fontsize',14);
hold on
%xlim([0.97 1.35]);
%ylim([0.97 1.35]);
yl=[0,0];
plot(get(gca,'xlim'),yl,'b:','LineWidth',1.5);
c = get(gca, 'Children');
legend(c(2:8),'DPR','BayesB','BayesA','Bayes LASSO','BSLMM','rrBLUP','BayesC');
%clear h
```

#### If you are interested in our model or any question and omission (bug). Please feel free to contact us.

* [Meng Luo](https://github.com/mengluoML)


