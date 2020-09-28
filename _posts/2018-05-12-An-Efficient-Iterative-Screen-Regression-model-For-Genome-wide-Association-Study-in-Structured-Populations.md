---
layout: single
author_profile: true
comments: true
enable_mathjax: true
output: html_document
title : "An Efficient Iterative Screen Regression model For Genome wide Association Study in Structured-Populations"
enable_mathjax: true
output: 
html_document:
   md_document:
       variant: markdown_github
category: [BMD_man,BMD_qq]
tags: [GWAS, ISR, Mult-loci model, Single-locus model,Simulation]
---


<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

## Iterative Screen Regression model (ISR)

Here we introduce a new unique variable selection procedure of regression statistic method, call Iterative screening regression(ISR). Where we formulated a new regression information criterion (RIC) and used this criterion as the objective function of the entire variable screen process. We evaluate various model selection criteria through simulations, which suggest that the proposed ISR method performs well in terms of FDR and power. Finally, we show the usefulness of our approach by applying it to A. thaliana and mouse data.

<div align="center"><img src="{{ "/images/Blog/GWAS/ISRGWAS.jpg" | prepend: site.baseurl }}"></div>

<p style="text-align: center;"> Schematic overview of model-based iterative screen regression for GWAS </p>

## Simulation

Human dataset derived from [PLINK](http://gigadb.org/dataset/view/id/100094/) included two real human genotype datasets, the first dataset included 1000 samples and 100000 makers (SNPs) over all chromosomes. The second included 10000 samples(6000 cases and 4000 control) and 88058 markers (SNPs), and only included in 19, 20, 21, and 22 chromosomes. Also, another outbred [CFW](https://datadryad.org/resource/doi:10.5061/dryad.2rs41) (Carworth Farms White) mice population that including a set of 92,734 single-nucleotide polymorphism markers which were genotyped 1,161 individuals were also used to perform one simulation experiments. well, all simulation both setting the heritability was 0.5. The first 100 phenotypes used all CFW mice dataset that including 100 markers were randomly selected as causal loci, respectively. We also assigned an additive effect randomly drawn from a standard normal distribution and added a random environmental term, where the  $h^2$ of the simulated traits only was 0.5, here. The second and third 100 phenotypes used human dataset that also including 100 markers were randomly selected as causal loci, respectively. Where the $h^2$ of the simulated traits (complex traits) only was 0.5.



<img src="https://render.githubusercontent.com/render/math?math=Y_j=\displaystyle\sum_{i=1}^{100} X_i\beta_i+ \varepsilon,\varepsilon \backsim MVN_n(0,\sigma_g^2((1-h^2)/h^2)),j=1,2,3,....1161.">


Mice:
$$Y_j=\displaystyle\sum_{i=1}^{100} X_i\beta_i+ \varepsilon,\varepsilon \backsim MVN_n(0,\sigma_g^2((1-h^2)/h^2)),j=1,2,3,....1161.$$



Human:
$$Y_j=\displaystyle\sum_{i=1}^{100} X_i\beta_i+ \varepsilon,\varepsilon \backsim MVN_n(0,\sigma_g^2((1-h^2)/h^2)),j=1,2,3,....1000.$$



Human:
$$Y_j=\displaystyle\sum_{i=1}^{100} X_i\beta_i+ \varepsilon,\varepsilon \backsim MVN_n(0,\sigma_g^2((1-h^2)/h^2)),j=1,2,3,....10000.$$


## Power versus FDR and TPR (Type one error)

How to define the power versus FDR and TPR, just saw [here](https://en.wikipedia.org/wiki/Sensitivity_and_specificity). As following was the confusion matrix.

* $$TPR=\frac{TP}{TP+FN}$$,True positive rate (Power & Sensitivity)
* $$FDR=\frac{FP}{TP+FP}$$,False discovery rate
* $$FPR=\frac{FP}{FP+TN}$$,False positive rate (1-Specificity)


<div align="center"><img src="{{ "/images/Blog/GWAS/power.jpg" | prepend: site.baseurl }}"></div>

### Receiver Operating Characteristic curves of simulations
Receiver Operating Characteristic (ROC) curves, demonstrating the trade-off between false positives and true positives, using simulated data based on the RA data. Lines trending toward the upper left corner denote better ROC curves while those in the opposite direction correspond to poorer methods.

### Mice simulation result
A receiver operating characteristic curve was a plot of the statistical power against the controlled FDR and FPR(type I error). Where we also calculated the AUC value, smaller AUC value indicated that method batter perform.


<div align="center"><img src="{{ "/images/Blog/GWAS/FDR_TPIALL.jpg" | prepend: site.baseurl }}"></div>

### Human simulation result

<div align="center"><img src="{{ "/images/Blog/GWAS/humanpower.jpg" | prepend: site.baseurl }}"></div>
Statistical power was defined as the proportion of simulated markers detected at cost defined by either False Discovery Rate (FDR) or False Positive Rate (FPR, Type I error). a The two types of Receiver Operating Characteristic (ROC) curves are displayed separately for TPR (true positive rate, power) versus FDR and FPR (the two simulations of Scenarios Ⅴ (1-2)). b The Area Under the Curves (AUC) are also displayed separately for TPR (true positive rate, power) versus FDR and FPR for 100 simulations. Four GWAS methods (ISR, FarmCPU, FaSTLMM, and PLINK-Fisher) were compared with phenotypes simulated from real genotypes in human. The simulated phenotypes had a heritability of 50%, controlled by 100 SNPs. These markers were randomly sampled from the available 100000 (88025) Single Nucleotide Polymorphism (SNPs). b .To specify the multiple comparison procedures using Least Significant Difference (LSD) after ANOVA. Here, $$*$$ represent a significant level of 0.05; $$**$$ represent a significant level of 0.01; $$***$$ represent a significant level of 0.001.

## Estimated Effect (PVE)

Computing proportion of variance in phenotype explained by a given SNP (PVE) :

$$PVE=\frac{\beta^2Var(X)}{ Var(Y)}=\frac{\beta^2Var(X)}{\beta^2Var(X)+\sigma^2}\approx \frac{Var(X\beta)}{ Var(Y)}$$

### Mice

<div align="center"><img src="{{ "/images/Blog/GWAS/MICEPVE.jpg" | prepend: site.baseurl }}"></div>

a The distribution of all simulated effects (all true effect) and the distribution of effects of loci identified (100 casual loci within 100 simulations, and only true positive) by six methods. The solid line shows the effect size by different methods and the phenotype with 50% of PVE. The bottom boxplot was explained the variance of the loci effect estimated by ISR, FASTmrEMMA, FarmCPU, GEMMA, and CMLM within the 100 simulations. Performance of estimating PVE is measured by the root of mean square error (RMSE), where a lower value indicates better performance. The true PVEs are shown as the horizontal dash lines.
### Human

<div align="center"><img src="{{ "/images/Blog/GWAS/HUMANPVE.jpg" | prepend: site.baseurl }}"></div>
a The explained variance of the casual loci effects estimated by ISR, FarmCPU, FaSTLMM, and PLINK-Fisher within the 100 simulations (The two simulations of Scenarios V (1-2)). b The distribution of all simulated effects (True Effect, black line) and the distribution of effects of loci identified (after 0.05 Bonferroni correction) by ISR (906 loci and 3461 loci), FarmCPU (760 loci and 3037 loci), FaSTLMM (433 loci and 2297 loci) and PLINK-Fisher (446 loci and 4041), respectively (The two simulations of Scenarios V (1-2)).
## Real dataset
### LD(Linkage Disequilibrium)
Locuszoom plots for genome-wide significant SNP both found by ISR and FASTmrEMMA. The locuszoom plots showing the zoom in of the most significant SNP both found by FASTmrEMMA and ISR two methods. The points for each SNP are colored by the level of the linkage disequilibrium (r2) with the index SNP, the SNP with the highest association to the quantitative trait meristem zone length. Following was the matlab code for locuszoom plot.

```matlab
clear,clc
warning off
addpath D:\\MATLAB\\isreg\\color\\cbrewer2\\cbrewer2
addpath D:\\MATLAB\\isreg\\arrow
addpath D:\\MATLAB\\isreg\\color\\colorspace\\colorspace
load meritem_201.mat
%tit='Chromosome 4';
subplot(4,1,1)
bin=2500000;%windon size 
p1=max(po);p2=min(po);
p=round((p1-p2)/bin);pl=zeros(length(po),1);
%for i=1:1
%a=0;
pol=po;
pm=zeros(p,1);
for j=1:p
    [a,b]=find(po<(po(1)+bin));
    pm(j)=length(a);
    col=zeros(p,3);
        %pl(a>pm(j))=j;%
    cl=unifrnd(0,1,1,3);
    col(j,:)=cl;
    stem(po(a),.05*ones(size(po(a))),'Color',cl,'Marker','none');
    %bar(po(a),.05*ones(size(po(a))),'FaceColor',cl,'LineWidth',.2)
        %a=a+0;
    po(a)=[];  
    hold on
end
xlim([9.45*1000000 9.565*1000000]);
pm=num2str(pm);
lgd = legend(pm,'Location','north','Orientation','horizontal');
title(lgd,'Plotted SNPs');
%h2=stem(po,.05*ones(size(po)),'Color','k','Marker','none');
set(gca,'Xtick',[],'Ytick',[]);ylim([0 0.05]);%ylabel('plotted SNPs')
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',14);

subplot(4,1,2)
fpml=-log10(f_pml);cb=40;
scatter(pol/1000000,fpml,cb,ld,'filled','LineWidth',2)
%colormap jet% the catgory of colormap{parula,jet,hsv,hot,cool,spring,summer,
%autumn,winter,gray,bone,copper,pink,lines,colorcube,prism,flag,white}
cm=colorbar('location','east','AxisLocation','in');
cm.Label.String = 'r^2';
%new color type
colormap(cbrewer2('seq','RdBu',6));%PuBu,YIGnBu,RdBu,RdYlBu
box on;
set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
ylabel('-log_{10}(\itP)');%title(tit);
set(gca,'XTick',[]);
xlim([9.45 9.565]);ylim([-0.2 4]);

subplot(4,1,3)
mpml=-log10(m_pml);
scatter(pol/1000000,mpml,cb,ld,'filled','LineWidth',2)
%colormap jet% the catgory of colormap{parula,jet,hsv,
%hot,cool,spring,summer,autumn,winter,gray,bone,copper,
%pink,lines,colorcube,prism,flag,white}
cm=colorbar('location','east');
cm.Label.String = 'r^2';
%new color type
colormap(cbrewer2('seq','RdBu',6));
box on;
set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
ylabel('-log_{10}(\itP)');%title(tit);
set(gca,'XTick',[]);
xlim([9.45 9.565]);ylim([-3 40]);


subplot(4,1,4);
for i=1:length(lab)
    a1=[al(i),al(i)];
    ps=px(i,:)/1000000;
    %line(ps,a1,'LineWidth',7);
    %annotation('arrow',)
    %dim=[ps,a1];
    %annotation(hf, 'arrow', dim);
    cl=unifrnd(0,1,1,3);
    text(mean(ps)-0.005,a1(1)+0.4,lab(i),'Color',cl,'FontName','Times New Roman',...
             'FontWeight','bold','FontSize',10);
    %davinci( 'arrow', 'X', ps, 'Y',a1,'ArrowType', 'double','LineWidth',2);%,'Shaft.Type','rectangle'
    %daspect( [1 1 1] ) 
    davinci( 'arrow', 'X', ps, ...
                  'Y',a1, ...
                  'ArrowType','single',...%double
                  'Shaft.Type','rectangle', ...
                  'Shaft.Width',0.3,... #'Head.Length',0.03,...
                  'Head.Width',0.3, ...# 'Head.Sweep',0, ...
                  'Color',cl, ...
                  'EdgeColor','k', ...
                  'LineWidth',0.5 );
    hold on
    %setting the right or left arrow (if 0 or 1)
    %text(ps,a1,'\rightarrow','FontSize',14,'FontWeight','bold')
end
xlim([9.45 9.565]);
set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
%set(gca,'xlim',[pos(n1),pos(n2)]);
set(gca,'YTick',[]);
xlabel('Position on Chromosome 3 (Mb)');set(gca,'Box','on','XGrid','on');   

```
<div align="center"><img src="{{ "/images/Blog/GWAS/LOCUSZOOM.jpg" | prepend: site.baseurl }}"></div>

### Manhattan plot

Here was how I have coded for BMD GWAS Manhattan plot.

```matlab
load manha.mat
[nsnp,phe]=size(mlssr);nchr=length(unique(chr));
figure("Position",[0 200 1000 500])
%nasnp=label(lm);
%dims = size(nasnp);
k=1;
f(1,1)=mlssr(1,2)+2000000;%if you used the genetic map, just only changed the number 2000000 to 
for i=2:1:nsnp
   if mlssr(i,1)~= mlssr(i-1,1)
    f(i,1)=f(i-1,1)+mlssr(i,2)+2000000;
    x(k)=f(i,1)-1000000;
    if k>1
    xg(k)=(f(i-1,1)-x(k-1))/2+x(k-1);
    ch(k)=mlssr(i-1,1);
    else
    xg(k)=f(i-1,1)/2-mlssr(1,2)/2;
    end
    k=k+1;
   else
   f(i,1)=mlssr(i,2)-mlssr(i-1,2)+f(i-1,1);
   ch(k)=mlssr(i-1,1);
    end
end
xg(k)=(f(i,1)-x(k-1))/2+x(k-1);
ch(k)=mlssr(i,1);
for j=3:1:phe
for i=1:1:nsnp
    f(i,j-1)=-log10(mlssr(i,j));
end
i=1;
end
j=1;
for j=1:1:(phe-2)
    figure(j);
    cl=load('colorchrhg.txt');ncl=randperm(length(cl(:,1)));cl=cl(ncl,:);
    %cl=cl(1:nchr,:);% all kinds of color
    scatter(f(:,1),f(:,j+1),10,cl(chr,:),'o','filled');%cl(mod(mlssr(:,1),nchr)+1,:)/1
    %display the gene name or not
    %text(f(lm,1),f(lm,j+1)+.15*4, nasnp(:,1:dims(2)),'FontSize', 10);
    chrs = unique(chr);
    chrs = chrs(:)';%p0=0;
  for c = chrs
      % Plot the SNPs on the chromosome.
      is = find(chr == c);
      maxpos = max(pos(is));
     % load cl.mat
      if ~isodd(c)
          %clr = unifrnd(0,1,1,3);
          clr=cl(c,:);
      else
          %clr = unifrnd(0,1,1,3);
          clr=cl(c,:);
      end
      % Add the chromosome number.
      %for il=1:nchr
      text(xg(c),-0.15* (max(-log10(plm)) - min(-log10(plm))),num2str(c),...
          'Color',cl(c,:),'FontSize',9,'HorizontalAlignment','center',...
          'VerticalAlignment','bottom','FontName','Times New Roman',...
          'FontWeight','bold','FontSize',14);
       %p0 = p0 + maxpos;
  end
  set(gca,'xlim',[f(1,1),f(nsnp,1)]);
  set(gca,'XTick', []);
  %set(gca,'XTickLabel',ch);
  set(gca,'TickDir','in');
  set(gca,'TickLength',[0.005,0.1]);
  set(gca,'FontName','Times New Roman', 'FontWeight','bold','FontSize',14);
  xlabel('Chromosome');ylabel('-log_{10}(\itP)');
  hold on;
  %add the signifigant line
  ylm=-log10(7.7e-06);yl=[ylm,ylm];
  ylm1=-log10(0.05/nsnp);yl1=[ylm1,ylm1];
  plot(get(gca,'xlim'),yl,'r:','LineWidth',1.5);
  plot(get(gca,'xlim'),yl1,'k:','LineWidth',1.5);
 text(xg(19),-log10(7.7e-06),num2str(7.7e-06))
 text(xg(19),-log10(0.05/nsnp),num2str(5.39E-07))
end
```
## R code for permutation in FarmCPU
We calculated a significance threshold via permutation, which is a standard approach for QTL mapping in mice that controls for the type I error rate (P<0.1, THE RED DASHLINE IN VISUILIZATION PLOT).
```r
p <- matrix(NA,92734,1000)
time <- Sys.time()
for(i in 1:1000){
  n <- length(Y[,1])
  test <- sample(n)
  y <- Y[test,]
  lmFarmCPU <- FarmCPU(Y = y, GD = genm, GM = mapm,threshold.output = 22)
  p[,i] <- lmFarmCPU$GWAS$P.value
  graphics.off()
}
write.table(p,'permutation_BMD_P.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
t <- f(time)
t
p <- read.table('permutation_BMD_p.txt',header = FALSE)
#save(genm,mapm1 file="mice_1160_92734.RData")
P <- NULL
for(i in 1:1000){
  P[i] <- min(p[,i])
}

#histplot
p <- read.csv('permutation_BMD_P.csv',header = FALSE)
p <- c(p)
p$x
hist(p$V1, breaks=30, col="gray")
# P<0.1
abline(v=sort(p$V1)[900],col='red')

```

<div align="center"><img src="{{ "/images/Blog/GWAS/bmdman.png" | prepend: site.baseurl }}"></div>

### QQ plot
Error bars in quantile-quantile plots. The quantile-quantile plot 95% confidence intervals shown in plots were computed by assuming M independent p-values, uniformly distributed on [0,1]. The $$k^{th}$$ largest P value from this distribution is the  $$k^{th}$$ order statistic, which is known to have a Beta $$(k, M-k-1)$$  distribution. The mean of this distribution is $$\frac{k}{M+1}$$. Thus, at point $$log(\frac{k}{M+1})$$ on the horizontal and vertical axis, error bars can be plotted as the $$(log)95%$$ probability interval of this Beta distribution.

Here was my matlab coded for BMD GWAS QQplot.

```matlab
load qq.mat
o = -log10(sort(plm,'ascend'));
alphaLevel=0.05;halflevel=0.5;oneMinusAlphalevel=1-alphaLevel;
x=1:1:length(o);
e = -log10( (x-0.5)./(length(o)))';
cl=unifrnd(0,1,24,3);
nclm=randperm(length(cl));
scatter(e,o,30,cl(nclm(1),:),'o','filled')
hold on
M=length(o);
mSeq=10.^(log10(0.5):0.1:(log10(M - 0.5)+0.1));
n=length(mSeq);n1=length(mSeq);alpha=ones(n,1);half=ones(n,1);oneMinusAlpha=ones(n,1);
if alphaLevel==alphaLevel
    for i=1:n
        alpha(i)=betainv(alphaLevel,mSeq(i),M-mSeq(i));
    end
end
if halflevel==halflevel
    for i=1:n
        half(i)=betainv(halflevel,mSeq(i),M-mSeq(i));
    end
end
if oneMinusAlphalevel==1-alphaLevel
    for i=1:n
        oneMinusAlpha(i)=betainv(1-alphaLevel,mSeq(i),M-mSeq(i));
    end
end
%alpha=[alpha half oneMinusAlpha];
betaDown=half-alpha;betaUp=oneMinusAlpha-half;
theoreticalPval=(mSeq/M)';
lowerBar=-log10(theoreticalPval-betaDown);%lowerBar=real(lowerBar);
upperBar=-log10(theoreticalPval+betaUp);
yBar=-log10(theoreticalPval);
g=[0 0.5 0];
line(yBar,lowerBar,'Color',g,'LineStyle','-.','LineWidth',2.5)
line(yBar,upperBar,'Color',g,'LineStyle','-.','LineWidth',2.5)
line(e,e,'Color','k','LineStyle','-','LineWidth',2.5)
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',16);
xlabel('Expected -log_{10}(\itP)','FontName','Times New Roman','FontWeight','bold','FontSize',16);
ylabel('Observed -log_{10}(\itP)','FontName','Times New Roman','FontWeight','bold','FontSize',16);
lg=legend('ISR','95% CI');
set(lg,'FontName','Times New Roman','FontWeight','bold','FontSize',12);
```

<div align="center"><img src="{{ "/images/Blog/GWAS/bmdqq.png" | prepend: site.baseurl }}"></div>

### Calculated time
Comparison of computing time of ISR and others methods in the third simulation scenarios. The average of computing times using five methods (FarmCPU-R, GEMMA-C++, CMLM-R, MLMM-R, and FASTmrEMMA-R) are compared with ISR-M (MATLAB language) for 100 simulations. The dataset containing 1161 individuals genotyped with 20000 markers. Simulation studies using a computer with an Inter(R) i3-6100 @3.7GHz CPU in windows 10. 

<div align="center"><img src="{{ "/images/Blog/GWAS/bar_time.jpg" | prepend: site.baseurl }}"></div>
 
#### If you are interested in our model or any question and omission (bug). Please feel free to contact us.

* [Meng Luo](https://github.com/mengluoML)
* [Shiliang Gu](http://www.wheatlab-yzu.com/article_show.asp?id=2184)


