---
layout: single
author_profile: true
comments: true
enable_mathjax: true
output: html_document
title : "Optimization by a new improved contraction expansion algorithm and its application"
toc: false
tags: [Nonlinear equation, Parameter estimation, Optimization,Contraction-expansion algorithm, Numerical differentiation, Fitting]
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

## Fitting with a curve and surface

> fitting a statistical/mathematical model to data finds its application in almost all empirical sciences, viz., physics, chemistry, zoology, botany, environmental sciences and economics, etc. That the processing of data in experimental methods, i.e., such as UV spectroscopy, X-ray analysis, IR spectrometry or chromatographic techniques, etc., is a key part of any research work, in which curve and surface ﬁtting has become a standard in the last few years.

> It has four objectives: the first, to describe the observed (or experimentally obtained) dataset by a statistical/mathematical formula that select the appropriate model; the second, to estimate the parameters of the formula so obtained and interpret them, which the interpretation is consistent with the generally accepted principles of the discipline concerned; the third, to predict, interpolate or extrapolate the expected values of the dependent variable with the estimated formula; and the last, to use the formula for designing, controlling or planning.

> Curve or surface fitting, called as parameter estimation of nonlinear equation, is a very complex nonlinear scheme, which hardly can be achieved for global optimization and optimum estimation. With the deep study of variables quantification and the universal application of non-linear regression analysis, it is very necessary to propose new algorithms of curve and surface fitting and nonlinear fitting software.

## Contraction-expansion algorithm

### Contraction stage (C stage)

> Contraction stage constitutes an initial searching domain begin with initial values and shrinks the searching range constantly to get the optimal parameters.


<div align="center"><img src="{{ "/images/Blog/fitting/C.jpg" | prepend: site.baseurl }}"></div>

> Calculate the searching points in turn to find the step points with minimum Q-value as the next searching center points$$b_j^{(1)}$$.

<div align="center"><img src="{{ "/images/Blog/fitting/C1.jpg" | prepend: site.baseurl }}"></div>

- Because Q usually is not the monotone function of b, it usually cannot get the optimal parameters after contraction stage, and it always turns up $$Q_{min}^{(C)}>Q_(min)$$.

<div align="center"><img src="{{ "/images/Blog/fitting/4.gif" | prepend: site.baseurl }}"></div>

- Step points method’s improvement

<div align="center"><img src="{{ "/images/Blog/fitting/improve.jpg" | prepend: site.baseurl }}"></div>

- Numerical differentiation
  It is hard for us to get the partial derivative function, so it isn’t suitable to put it become the first condition of curve and surface fitting. If we use numerical differentiation, we can get the estimation value of the partial derivative function in b .

<div align="center"><img src="{{ "/images/Blog/fitting/ND.jpg" | prepend: site.baseurl }}"></div>

When the surface and curve fitting involve multi-parameter, the approximate partial derivative of the j-th parameter of the formula in the optimal parameter point $$b^{(0)}$$：

<div align="center"><img src="{{ "/images/Blog/fitting/b.jpg" | prepend: site.baseurl }}"></div>

## Improved Gauss-Newton:

$$ \vartriangle = A^{-1}K \Longrightarrow \vartriangle \Longrightarrow b = b^{(0)} + \vartriangle $$

## 8 datasets with higher level of difficulty

<div align="center"><img src="{{ "/images/Blog/fitting/table.jpg" | prepend: site.baseurl }}"></div>

### The certified results of NIST and estimated by improved C-E algorithm for 7 datasets

<div align="center"><img src="{{ "/images/Blog/fitting/table1.jpg" | prepend: site.baseurl }}"></div>

### The curve fitting figure of result MGH10, Thurber and Eckerle4

<div align="center"><img src="{{ "/images/Blog/fitting/FIG1.jpg" | prepend: site.baseurl }}"></div>

### Global and local optimal for the problem (7)
  <div align="center"><img src="{{ "/images/Blog/fitting/surfacetable.jpg" | prepend: site.baseurl }}"></div>

<div align="center"><img src="{{ "/images/Blog/fitting/surface.jpg" | prepend: site.baseurl }}"></div>

### The eighth test problem

<div align="center"><img src="{{ "/images/Blog/fitting/test8.jpg" | prepend: site.baseurl }}"></div>

<div align="center"><img src="{{ "/images/Blog/fitting/test8_1.jpg" | prepend: site.baseurl }}"></div>

- The surface fitting plot of MathWorks Peaks function simulation data
  <div align="center"><img src="{{ "/images/Blog/fitting/sda.jpg" | prepend: site.baseurl }}"></div>

<div align="center"><img src="{{ "/images/Blog/fitting/sda1.jpg" | prepend: site.baseurl }}"></div>

> _Some of those test data are very hard, and may never get right answers without using 1stOpt(Auto2Fit). Even for 1stOpt(Auto2Fit), it does not ensure every run will be successful. … In some cases, you may try to change the control parameter of ‘Population Size’ (1stOpt team said)_ They have suggested that to solve these problems one should use Global Levenberg-Marquard or Global BFGS method.

## Others application

- A phenotype (trait) that result from Mutiple normal distribution superposition (multiple gene result in phenotype variation).

```matlab
% Coded for Mutiple normal distribution superposition
clear,clc
x=xlsread('litian','b1:b2170');
g=input('number of groups, g= ')
rg=max(x)-min(x);
ii=rg/(g-1);
xi=linspace(min(x)-.5*ii,max(x)+.5*ii,g+1);
subplot(2,2,1);
xi1=xi-ii/4;
hist(x,xi)
axis tight
subplot(2,2,2)
xi2=xi;
hist(x,xi2)
axis tight
subplot(2,2,3)
xi3=xi+ii/4;
hist(x,xi3)
axis tight
subplot(2,2,4)
xi4=xi+ii/2;
hist(x,xi4)
axis tight
sel=input('put your selection, sel= ')
figure(1),clf
if sel==1
    xi=xi1;
elseif sel==2
    xi=xi2;
elseif sel==3
    xi=xi3;  
elseif sel==4
    xi=xi4;
end
f=hist(x,xi);
if f(1)==0;
    f(1)=[];xi(1)=[];
end
g=size(f,2);
if f(g)==0
    f(g)=[];xi(g)=[];
end
g=size(f,2);
f,xi
hist(x,xi)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
%axis tight
pause(5)
hold on
%f1=[0,f,0];
%x1=[xi(1)-ii,xi,xi(g)+ii];
%plot(x1,f1,'r-','linewidth',3)
%axis tight
%stop
fx=@(b,x)b(1)*exp(-b(2)*(x-b(3)).^2)+b(4)*exp(-b(5)*(x-b(6)).^2)+b(7)*exp(-b(8)*(x-b(9)).^2)+b(10)*exp(-b(11)*(x-b(12)).^2)+b(13)*exp(-b(14)*(x-b(15)).^2);
fx1=@(b,x)b(1)*exp(-b(2)*(x-b(3)).^2);
fx2=@(b,x)b(4)*exp(-b(5)*(x-b(6)).^2);
fx3=@(b,x)b(7)*exp(-b(8)*(x-b(9)).^2);
fx4=@(b,x)b(10)*exp(-b(11)*(x-b(12)).^2);
fx5=@(b,x)b(13)*exp(-b(14)*(x-b(15)).^2);
bf24=[175.4829698,66.74574337,0.276553824,202.8109964,3123.739144,0.032487556,33.22388584,2264.479333,0.391977911,319.6439187, 956.9020534,0.064620109,57.21979178,212.6365602,0.157921054];
b=bf24;
%figure(1),clf
%hold on
x1=-0.05:.000001:0.55;
y=fx(b,x1);
y1=fx1(b,x1);
y2=fx2(b,x1);
y3=fx3(b,x1);
y4=fx4(b,x1);
y5=fx5(b,x1);
plot(x1,y,'k-','linewidth',3)
plot(x1,y1,'r-','linewidth',2)
plot(x1,y2,'b-','linewidth',2)
plot(x1,y3,'m-','linewidth',2)% five
plot(x1,y4,'g-','linewidth',2)%m ,two
plot(x1,y5,'c-','linewidth',3)%three
title('FHB')
axis([-0.00010,0.55,0,370])
N=sum(f);
disp([ '    miu,          sigma2,            sigma ,           N'])
rst(1,1)=b(3);rst(1,2)=.5/b(2);rst(1,3)=sqrt(rst(1,2));rst(1,4)=N*b(1)/(b(1)+b(4)+b(7)+b(10)+b(13));
rst(2,1)=b(6);rst(2,2)=.5/b(5);rst(2,3)=sqrt(rst(2,2));rst(2,4)=N*b(4)/(b(1)+b(4)+b(7)+b(10)+b(13));
rst(3,1)=b(9);rst(3,2)=.5/b(8);rst(3,3)=sqrt(rst(3,2));rst(3,4)=N*b(7)/(b(1)+b(4)+b(7)+b(10)+b(13));
rst(4,1)=b(12);rst(4,2)=.5/b(11);rst(4,3)=sqrt(rst(4,2));rst(4,4)=N*b(10)/(b(1)+b(4)+b(7)+b(10)+b(13));
rst(5,1)=b(15);rst(5,2)=.5/b(14);rst(5,3)=sqrt(rst(5,2));rst(5,4)=N*b(13)/(b(1)+b(4)+b(7)+b(10)+b(13));
disp(rst)
```

<div align="center"><img src="{{ "/images/Blog/fitting/litian.jpg" | prepend: site.baseurl }}"></div>

## All kinds of Models perfectly solved by CE-Algorithm 

<div align="center"><img src="{{ "/images/Blog/fitting/5.gif" | prepend: site.baseurl }}"></div>


#### If you are interested in our algorithm, or Found any omission (bug) . Please feel free to contact us.

- [Meng Luo](https://github.com/mengluoML)
- [Shiliang Gu](http://www.wheatlab-yzu.com/article_show.asp?id=2184)


