---
layout: single
author_profile: false
comments: true
enable_mathjax: true
output: html_document
title : "Solve traveling salesman problem using EMF CE algorithm"
tags: [Contraction-Expansion algorithm (CE), Exchange-Move-Flip (EMF), Combinatorial Optimization,traveling salesman problem (TSP)]
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

## Traveling salesman problem (TSP)

The traveling salesman problem (TSP) is well known the classical and fundamental NP-hard combinatorial optimization problems. The classical TSP that can be described as following: find a path through a weighted graph that starts and ends at the same city, includes every other city exactly once, and minimizes the total distance tour of n cities. [MORE ABOUT](https://en.wikipedia.org/wiki/Travelling_salesman_problem)

## TSPs Instance

### Art TSPs

| [![da Vinci's Mona Lisa](http://www.math.uwaterloo.ca/tsp/data/ml/mona-lisa100K.gif)](http://www.math.uwaterloo.ca/tsp/data/ml/mona-lisa100K.gif) | [![van Gogh's Self Portrait 1889](http://www.math.uwaterloo.ca/tsp/data/art/vangogh120K.gif)](http://www.math.uwaterloo.ca/tsp/data/art/vangogh120K.gif) | [![Botticelli's The Birth of Venus](http://www.math.uwaterloo.ca/tsp/data/art/venus140K.gif)](http://www.math.uwaterloo.ca/tsp/data/art/venus140K.gif) |
| :-----------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------------------------------------------: |
|                                  [mona-lisa100K.tsp](http://www.math.uwaterloo.ca/tsp/data/ml/mona-lisa100K.tsp)                                  |                                       [vangogh120.tsp](http://www.math.uwaterloo.ca/tsp/data/art/vangogh120K.tsp)                                        |                                        [venus140K.tsp](http://www.math.uwaterloo.ca/tsp/data/art/venus140K.tsp)                                        |

| [![Velazquez's Juan de Pareja](http://www.math.uwaterloo.ca/tsp/data/art/pareja160K.png)](http://www.math.uwaterloo.ca/tsp/data/art/pareja160K.png) | [![Courbet's The Desperate Man](http://www.math.uwaterloo.ca/tsp/data/art/courbet180K.png)](http://www.math.uwaterloo.ca/tsp/data/art/courbet180K.png) | [![Vermeer's Girl with a Pearl Earring](http://www.math.uwaterloo.ca/tsp/data/art/earring200K.gif)](http://www.math.uwaterloo.ca/tsp/data/art/earring200K.gif) |
| :-------------------------------------------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|                                      [pareja160K.tsp](http://www.math.uwaterloo.ca/tsp/data/art/venus140K.tsp)                                      |                                      [courbet180K.tsp](http://www.math.uwaterloo.ca/tsp/data/art/courbet180K.tsp)                                      |                                          [earring200K.tsp](http://www.math.uwaterloo.ca/tsp/data/art/earring200K.tsp)                                          |

### National TSPs

| [![China - 71,009 Cities](/images/Blog/TSP/national/chpoints.gif)](/images/Blog/TSP/national/chpoints.gif) | [![Egypt - 7,146 Cities](/images/Blog/TSP/national/egpoints.gif)](/images/Blog/TSP/national/egpoints.gif) | [![Greece - 9,882 Cities](/images/Blog/TSP/national/grpoints.gif)](/images/Blog/TSP/national/grpoints.gif) |
| :--------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------: |
|                [China - 71,009 Cities](http://www.math.uwaterloo.ca/tsp/world/ch71009.tsp)                 |                 [Egypt - 7,146 Cities](http://www.math.uwaterloo.ca/tsp/world/eg7146.tsp)                 |                 [Greece - 9,882 Cities](http://www.math.uwaterloo.ca/tsp/world/gr9882.tsp)                 |

### The best known results for the TSP Art instances are given in the table below. I would be happy to post any improvements you find!

<div align="center"><img src="{{ "/images/Blog/TSP/bestartTSP.jpg" | prepend: site.baseurl }}"></div>

## EMF-CE algorithm

Here, I will introdutes a novel search algorithm that based on Contraction-Expansion algorithm and integrated three operators Exchange, Move and Flip (EMF-CE) is proposed for the traveling salesman problem (TSP). EMF-CE uses a negative exponent function to generate critical value as the feedback regulation of algorithm implementation. Also, combined Exchange Step, Move step with Flip step and constitute of more than twenty combinatorial optimization of program elements. It has been shown that the integration of local search operators can significantly improve the performance of EMF-CE for TSPs. We test small and medium scale (51-1000 cities) TSPs were taken from the TSPLIB online library. The experimental results show the efficiency of the proposed EMF-CE for addressing TSPs in comparison with other state-of-the-art algorithms.

## METHODOLOGIES FOR SOLVING TSPs

#### _A. Evolutionary Algorithms_

Evolutionary algorithms (EAs) forebode a very promising direction. However, only general problem-independent EAs are frequently inefficient in solving TSPs, especially large TSPs. Tsai et al proposes an evolutionary algorithm, called the heterogeneous selection evolutionary algorithm (HeSEA), which integrates edge assembly crossover (EAX) and Lin–Kernighan (LK) local search, through family competition and heterogeneous pairing selection. Where the HeSEA is especially efficient to cope with the large TSPs.

#### _B. Simulated Annealing_

Simulated annealing (SA) is a probabilistic technique for approximating the global optimum of a given function. Specifically, it is a metaheuristic to approximate global optimization in a large search space. It is often used when the search space is discrete (e.g., all tours that visit a given set of cities). But the major shortage of SA is that it can be extremely slow and require significantly more processing time than other meta-heuristics. Several studies have tried to improve SAs performance by changing the generation and the acceptance mechanisms.

#### _C. Ant Colony Optimization Algorithms_

Ant colony optimization algorithm (ACO) is a probabilistic technique for solving computational problems which can be reduced to discovery good paths through graphs. However, for multifaceted large-scale optimization problems, the drawback of the ACO that was leisurely premature convergence.

#### _D. Heuristics And Neural Network_

Local search with k-change neighborhoods, k-opt, are the most widely used heuristic method for the traveling salesman problem. Helsgaun proposed a new LKH called LKH-2, and based on the LKH-1. which eliminates many of the limitations and shortcomings LHK-1. Combined use of general K-opt moves with six partitioning schemes, for obtaining high-quality solutions and reducing the complexity of solving large-scale problem instances, respectively. As we all known, LKH-2 is the best search algorithm to achieve the best optimal tour so far. Rego et al surveyed leading heuristics for the TSPs, which gathered with Lin–Kernighan (LK) and stem-and-cycle (S&C) methods, as well known the most effective and efficient local search ejection chain (EC) methods. Proposing six and two variants of the Lin–Kernighan (LK) and stem-and-cycle (S&C) ejection chain method, respectively.

### PROPOSED EXCHANGE-MOVE-FLIP BASED ON CE ALGORITHM

we proposed a novel search algorithm that based on we proposed Contraction-Expansion algorithm is introduced to provide the best compromise between the convergence speed and solution quality. The main idea of the proposed EMF-CE is taking advantage of integrated three operators Move, Exchange and Flip techniques and cooperated the CE to search the optimal solution of TSPs. We will use ten nodes as an example to illustrate the process of the EMF-CE algorithm, as following Fig show.

<div align="center"><img src="{{ "/images/Blog/TSP/CE/CE.png" | prepend: site.baseurl }}"></div>

## EXPERIMENTAL STUDY

In this section, we describe the test problems and implementation details of the algorithms that we use in our computational study. We present and discuss the results generated by others different algorithms that have reported from literatures on the same problems in TSPLIB. Also, in our experiments, which we presented new evaluation index described as following, and it was used to evaluate the performance of EMF-CE algorithm. And used the following index that were described as the efficiency index of EMF-CE algorithm.

## CONCLUSION

we presented a new method, which Contraction-expansion algorithm based on Move, Change and Flip three operators called EMF-CE for TSPs. The main purpose of designing EMF-CE is to achieve near to optimal solution quality and the most accurate algorithm while the method would give comparatively slower of convergence speed. Its computational intensity is O(n3). But it is now can achieve global optimal for the small scale TSPs. For the generalized traveling salesman problems that round number and round integer may affect the end optimal result. Usually can take 4 to 6 significant digit to calculate the distance and TSP path length, in order to prevent the loss of precision and affect the realization of the optimal solution.

The EMF-CE performance was compared with those of several computationally comparable including ASA-GS, GCGA, HGA, CONN and LBSA etc.. in terms of the qualitative comparison between those algorithms. However, we only show that five algorithm compare with the EMF-CE(have compared and reported in etc.). Cause the others state-of-the-art (GATM, GSM, SA, ITS, GSA-ACS-PSO, RABNET-TSP and MSA-IBS etc. have described and overviewed in section II) algorithm that can’t give the best tour lengths (reach to global optimal) for the small scale TSPs and let alone the medium scale TSPs.

<div align="center"><img src="{{ "/images/Blog/TSP/CE/time.png" | prepend: site.baseurl }}"></div>
 EMF-CE average CPU time versus the number of cities for 43 benchmark TSPs from TSPLIB
<div align="center"><img src="{{ "/images/Blog/TSP/CE/result.png" | prepend: site.baseurl }}"></div>
 Average shortest tour lengths percent above optimality for EMF-CE and five comparable algorithms

## EMF-CE CODE

```matlab
clear,clc
#inputdata
[Dimension,NodeCoord,NodeWeight,Name]=FileInput('lin318.tsp');%gr120,rat783
x=NodeCoord(:,2:3);
x1=x;
tspm=42220;
[n,k]=size(x);
cp=randperm(n);
x=x(cp,:);
rx=range(x(:,1));ry=range(x(:,k));
md=round((1.5+.27*n)/(1+.008*n));me=md+7;
mc=(2+.0063*n)/(1+.0007*n);
str=num2str([1:n]');
%figure(1),clf
% if k==2
%     plot(x(:,1),x(:,2),'o','markersize',8,'markerfacecolor','c','markeredgecolor','r','linewidth',2)
%     text(x(:,1)-.0075*rx,x(:,2)+.015*ry,str,'fontsize',10,'fontweight','b')
% elseif k==3
%     plot3(x(:,1),x(:,2),x(:,3),'o','markersize',10,'markeredgecolor','c','linewidth',2.5)
%     hold on
%     stem3(x(:,1),x(:,2),x(:,3),'filled')
%     text(x(:,1),x(:,2),x(:,3)+.02*ry,str,'fontsize',10,'fontweight','b')
%     grid on
% end
x(n+1,:)=x(1,:);Td=0;
%d=pdist2(x,x);d=round(d);
d=pdist2(x,x);d=round(d*100)/100;
for i=1:n
    Td=Td+d(i,i+1);
end
crt=.15*rand*Td/n; N=n*(n+1)/2; nc=0;

X=zeros(n,2,5);tp=zeros(1,5);
for i=1:5
    X(:,:,i)=x1;tp(i)=Td;
end
tmx=zeros(1,k);tmh=zeros(1,n+1);tml=zeros(n+1,1);
tpx=zeros(3,k);tph=zeros(3,n+1);tpl=zeros(n+1,3);
TSP=Td; disp(['TSP0=',num2str(Td)])
tsp=zeros(1,30*n);tsp(1:2)=Td;Xm=x;Dm=d;v1=2;rep=0;
pause(.01), t1=tic;
tic
while  rep<=150+1.5*n && tsp(v1)>tspm

    if rep>11 && mod(v1,md)==0
        n1=random('unif',n/29,n/12,1,md+6);
        while sum(n1)>1.04*n, n1=n1-1; end
        while sum(n1)>.96*n, n1(n1==min(n1))=[]; end
        n1=round(n1);
        if sum(n1)<.8*n
            n2=length(n1); n2=n2+2;n1(n2-1)=fix((n-sum(n1))/2);n1(n2)=n-sum(n1);
        else
            n2=length(n1); n2=n2+1;n1(n2)=n-sum(n1);
        end
        n1(n1<=0)=[];n2=length(n1);
        ncp=randperm(n2);m=1;b=zeros(n2,max(n1));
        for i=1:n2
            b(i,1:n1(i))=m:m+n1(i)-1;
            m=m+n1(i);
        end
        a=b(ncp,:);
        fp=randi(2,1,n2)-1;
        for i=1:n2
            if fp(i)==1
                fliplr(a(i,:));
            end
        end
        a=a'; a=a(:); a(a==0)=[];
        x=X(:,:,randi(5,1,1));
        x=x(a,:);
       % mi=mod(v1,5);
       % if mod(mi,3)==1
         %   y=[]; ind=kmeans(x,n2);
         %   for i=1:n2
            %    y=[y;x(ind==i,:)];
           % end
          %  x=y;
       % end
        x(n+1,:)=x(1,:); Td=0;
        %d=pdist2(x,x);d=round(d);
        d=pdist2(x,x);d=round(d*100)/100;
        for i=1:n
            Td=Td+d(i,i+1);
        end
        TSP=Td; crt=(.15+.3*rand)*crt; nc=0;
        disp([md,TSP,crt,n2])
    elseif mod(v1,me)==0 && rep>12
        crt=crt/(4+6*rand); x=Xm; d=Dm; TSP=tsp(v1);
    elseif mod(v1,5+fix(log(v1)))==4
        if mod(v1,2)==1, x=X(:,:,randi(5,1,1)); B=[x;x]; else B=[x0;x0]; end
        in=fix(random('unif',.3*n,.7*n,1));
        x=B(in:in+n-1,:); x(n+1,:)=x(1,:);
        %Td=0; d=pdist2(x,x);d=round(d);
       Td=0; d=pdist2(x,x);d=round(d*100)/100;
        for i=1:n
            Td=Td+d(i,i+1);
        end
        TSP=Td;
    end

    if mod(v1,9)==0
        rd=randperm(20);
        for i1=1:20
            switch rd(i1)
                case 1
                    ts1
                case 2
                    ts2
                case 3
                    ts3
                case 4
                    ts4
                case 5
                    ts5
                case 6
                    ts6
                case 7
                    ts7
                case 8
                    ts8
                case 9
                    ts9
                case 10
                    ts10
                case 11
                    ts11
                case 12
                    ts12
                case 13
                    ts13
                case 14
                    ts14
                case 15
                    ts15
                case 16
                    ts16
                case 17
                    ts17
                case 18
                    ts18
                case 19
                    ts17
                otherwise
                    ts18
            end
        end
        ts19
        pause(.1)
    else
        ts1,ts2,ts3,ts4,ts5,ts6,ts7,ts8,ts9,ts10,ts11,ts12,ts13,ts14,ts15,ts16,ts17,ts18,ts19
        pause(.1)
    end

    if mod(v1,2)==0
        crt=crt*exp(-.53*(sqrt(nc)-sqrt(.03*N+.15*n*mod(v1,n)+1.5*n*mod(v1,fix(n/mc))))/n);
        crt=.1*TSP/n*(crt<.1*TSP/n)+(crt>=.1*TSP/n)*crt*(crt<=6*TSP/n)+6*TSP/n*(crt>6*TSP/n);
        %disp([1,nc crt])
        ne=nc;
    else
        crt=crt*exp(-.5*(sqrt(nc)-sqrt(.05*N+.2*n*mod(v1,n)+3*n*mod(v1,fix(n/mc))))/n);
        crt=.1*TSP/n*(crt<.1*TSP/n)+(crt>=.1*TSP/n)*crt*(crt<=6*TSP/n)+6*TSP/n*(crt>6*TSP/n);
        %disp([2,nc,crt])
        ne=nc;nc=0;
    end
    x0=Xm;x0(n+1,:)=[];
    if tsp(v1)<tsp(v1-1)
        for i=5:-1:2
            X(:,:,i)=X(:,:,i-1);tp(i)=tp(i-1);
        end
        X(:,:,1)=x0;tp(1)=tsp(v1);rep=0;
    end
    rep=rep+1;v1=v1+1;tsp(v1)=tsp(v1-1);
    disp([v1-2,tsp(v1),crt,ne])
    pause(.1)
%     figure(2),clf
%     plot(x0(:,1),x0(:,2),'o','markersize',6,'markerfacecolor','c','markeredgecolor','r','linewidth',2)
%     text(x0(:,1)-.0075*rx,x0(:,2)+.017*ry,str,'fontsize',8,'fontweight','b')
%     hold on
%     plot(Xm(:,1),Xm(:,2),'b-','linewidth',2.5)
%     title(['tsp=',num2str(tsp(v1))],'fontsize',14)
%     pause(.1)
end
toc
t2=toc(t1)
disp(tsp(v1))
figure(2),clf
if k==2
    x0=Xm;x0(n+1,:)=[];
    plot(x0(:,1),x0(:,2),'o','markersize',8,'markerfacecolor','c','markeredgecolor','r','linewidth',2)
    text(x0(:,1)-.0075*rx,x0(:,2)+.015*ry,str,'fontsize',10,'fontweight','b')
    hold on
    plot(Xm(:,1),Xm(:,2),'b-','linewidth',2.5)
elseif k==3
    x0=Xm;x0(n+1,:)=[];
    plot3(x0(:,1),x0(:,2),x0(:,3),'o','markersize',10,'markeredgecolor','c','linewidth',3)
    stem3(x0(:,1),x0(:,2),x0(:,3),'filled')
    text(x0(:,1),x0(:,2),x0(:,3)+.017*ry,str,'fontsize',10,'fontweight','b')
    hold on
    plot3(Xm(:,1),Xm(:,2),Xm(:,3),'r-','linewidth',2)
    grid on
end
title(['tsp=',num2str(tsp(v1))],'fontsize',16)
axis tight
filename='ch150.xlsx';
%result={'Xih','bi','SEb','Up','F','p'; Tr(i,:); num2str(b(i+1)); num2str(SEb(i)); num2str(up(i)); num2str(f(i))};
x0;
xlswrite(filename,x0);
```
#### If you have any question for our algorithm, and Found any omission (bug) . Please feel free to contact us.

- [Meng Luo](https://github.com/mengluoML)
- [Shiliang Gu](http://www.wheatlab-yzu.com/article_show.asp?id=2184)




