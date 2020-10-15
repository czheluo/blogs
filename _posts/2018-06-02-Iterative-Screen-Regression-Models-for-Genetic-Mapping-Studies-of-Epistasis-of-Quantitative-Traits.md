---
layout: single
author_profile: false
comments: true
toc: true
enable_mathjax: true
output: html_document
title : "Iterative Screen Regression Models for Genetic Mapping Studies of Epistasis of Quantitative Traits"
enable_mathjax: true
output: html_document
tags: [GWAS, ISR, Epistasis, Genetic Mapping, Gene–gene interaction, Exhaustive Search]
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/latest.js?config=TeX-MML-AM_CHTML">
</script>

## Iterative Screen Regression model (ISR)

As mentioned [before](http://mengluocv.me/blog/2018/05/12/An-Efficient-Iterative-Screen-Regression-model-For-Genome-wide-Association-Study-in-Structured-Populations/) that our method, here, was used to Epistasis GWAS.

<div align="center"><img src="{{ "/images/Blog/EGWAS/ISRWGWAS.jpg" | prepend: site.baseurl }}"></div>

<p style="text-align: center;"> Schematic overview of model-based iterative screen regression for Epistasis GWAS </p>

## Simulations: Type I error control

To test the power of MAPIT, I again consider simulation designs similar to those proposed by previous [epistatic analysis studies](https://www.sciencedirect.com/science/article/pii/S0002929710003782). First,I assume that the broad-sense heritability is known ($$H^2 = 0.6$$). Next, I randomly chosed 1000~1008 from the chromosome of all control cases from the human data set (i.e. n ≈ 1,000~2,000 and p ≈ 1,000~1,008) to simulate continuous phenotypes that mirror genetic architectures affected by a combination of additive and pairwise epistatic effects. Specifically, I randomly choose different causal SNPs that directly affect the phenotype and classify the causal variants into TWO groups: (1) a set of additive SNPs; (2) a another set of interaction SNPs.
The additive effect sizes of all causal SNPs again come from a standard normal distribution or β ∼ MVN(0, I). Next, I create a separate matrix $$W$$ which holds the pairwise interactions of all causal the SNPs. These SNPs have effect sizes also drawn as α ∼ MVN(0, I). We scale both the additive and pairwise genetic effects so that collectively they explain a fixed proportion of genetic variance. Namely, the additive effects make up ρ%, while the pairwise interactions make up the remaining (1 − ρ)%. Once we obtain the final effect sizes for all causal SNPs, we draw errors to achieve the target $$H^2$$. The phenotypes are then created by summing all effects using simulation models: y = Xβ + Wα + ε, where $$X$$ is the additive genotype matrix, $$W$$ is the epistasis genotype matrix.

### I consider a few scenarios that depend on two parameters:

> (1 − ρ), which measures the portion of $$H^2$$ that is contributed by the interaction effects of the second groups of causal SNPs. Specifically, the phenotypic variance explained (PVE) by the additive genetic effects is said to be $$V(Xβ) = ρH^2$$, while the PVE of the pairwise epistatic genetic effects is given as $$V(Wα) = (1-ρ)H^2$$.

> P1/P2, which are the number of causal SNPs in each of the TWO groups respectively.

Specifically, we set ρ = {0.5, 0.8} and choose P1/P2 = 10/10 (scenario I), 50/10 (scenario II), 90/10 (scenario III), 10/50 (scenario IV), and 10/90 (scenario V).

### Estimating and identifying epistatic effects

<div align="center"><img src="{{ "/images/Blog/EGWAS/HUMAN.jpg" | prepend: site.baseurl }}"></div>
Empirical power to detect simulated causal interacting makers and estimating their marginal PVE. Groups 1 and 2 causal markers are colored in light red and light blue, respectively. These figures are based on a broad-sense heritability level of H2 = 0.6 and parameter ρ = 0.8, estimated with 100 replicates. Here, ρ = 0.8 was used to determine the portion of broad-sense heritability contributed by interaction effects. (A) shows the power of MAPIT to identify SNPs in each causal group under significance level α = 0.05. The lines represent 95% variability due to resampling error. (B) shows boxplots of the marginal PVE estimates for the group 1 and 2 causal SNPs from MAPIT for the four simulation scenarios. The true PVEs per causal SNP (0.012 for the group 1 SNPs; 0.012, 0.006, 0.0024, and 0.0012 for the Group 2 SNPs) are shown as dashed grey horizontal lines.

### Identifying pairwise interactions

<div align="center"><img src="{{ "/images/Blog/EGWAS/tepyoneerror.jpg" | prepend: site.baseurl }}"></div>
We compare the mapping abilities of MAPIT (solid line) to the exhaustive search procedure in PLINK (dotted line) in scenarios I, II, III, IV, and V under broad-sense heritability level $$H^2 = 0.6$$ and $$ρ = {0.5,0.8}$$. Here, $$ρ$$ was used to determine the portion of broad-sense heritability contributed by interaction and additive effects. Group 1 and group 2 causal SNPs with two different color. The x-axis shows the Type one error, while the y-axis gives the rate at which true causal variants were identified. Results are based on 100 replicates in each case.

## Real data applications

I assess ISE ability to detect epistasis in a quantitative trait loci (QTL) association mapping study for RICE.

### Detecting epistasis In Rice

Here, Manhattan 3D plot showed that the addtive and domain effects, and the interacted each other.

```matlab
%Coded 3D manhattan plot (sort of)
clear,clc
warning off
%map=readtable('all_yieldperplant1998(10chose)5.0.00001.txt','Delimiter','space');
load binmap_rice_add.mat
load yieldperplant1998(5chose)3.0.00001.d.all.mat
%start=start*1000000;stop=stop*1000000;
chl=table2array(all_result(:,1));
x1 = table2array(map(:,3));x2 = table2array(map(:,4));
x4 = table2array(map(:,5));x5 = table2array(map(:,6));
label1=table2cell(map(:,1));label2=table2cell(map(:,2));
label3=[label1 label2]; clear label1 label 2
x3 = table2array(map(:,8));an=find(x3==0);
x3(an)=0.1;plm=-log10(x3);%clear map
chr=unique(x1);ch=num2str(chr);
%for i=1:chr

% two type of PLOT (STEM OR SCATTER)
chose=input('Chose a type of epistasis plot(stem or scatter(1)) = ');
if chose==1
    for i=1:length(chr)
        [g,h]=find(x1==i);
        s=plm*5;c=x2;
        scatter3(x1(g),x2(g),plm(g),s(g),c(g),'filled');
        %scatter3(x4(g),x5(g),plm(g),s(g),c(g),'filled');
        hold on
    end
else
    %cl=zeros(length(chr),3);
    cl=load('colorchrhg.txt');ncl=randperm(length(cl(:,1)));cl=cl(ncl,:);
    for i=1:length(chr)
        [c,d]=find(x1==i);
        %cl1=cl(i,:);
        %cl1=unifrnd(0,1,1,3);
        %cl2=unifrnd(0,1,1,3);
        stem3(x1(c),x2(c),plm(c),'filled','Color',cl(i,:));clear x4 x5%,'MarkerEdgeColor',cl2);
        %stem3(x4(c),x5(c),plm(c),'filled','Color',cl1)%,'MarkerEdgeColor',cl2);
        hold on
        %cl(i,:)=cl1;
    end
    save col1.mat cl
end
set(gca,'XTick', chr);set(gca,'XTickLabel',ch);set(gca,'xlim',[1-0.5,length(chr)+0.5]);
set(gca,'YTick', chr);set(gca,'YTickLabel',ch);set(gca,'ylim',[1-0.5,length(chr)+0.5]);%ylim([1 length(chr)]);
set(gca,'FontName','Times New Roman','FontWeight','bold','FontSize',16);
xlabel('Chromosome');ylabel('Chromosome');zlabel('-log_{10}(\itP)');box on
[a,b]=find(x3<0.05/length(plm));
label=join(label3(a,:),'X');
 hold on
 if chose==1
     scatter3(x2(a),x1(a),plm(a),s(a),c(a),'filled')
     hold on
     stem3(x2(a),x1(a),plm(a),'Marker','none')
    for i=1:length(a)
        cl=unifrnd(0,1,1,3);
        text(x2(a(i)),x1(a(i)),plm(a(i))+1.5,label(i),'Color',cl,'FontName','Times New Roman',...
            'FontWeight','bold','FontSize',10);
    end
 else
     for i=1:length(a)
         stem3(x2(a(i)),x1(a(i)),plm(a(i)),'filled','Color',cl(x1(a(i)),:))
         hold on
     end
     for i=1:length(a)
         text(x2(a(i)),x1(a(i)),plm(a(i))+1.5,label(i),'Color',cl(x1(a(i)),:),'FontName','Times New Roman',...
             'FontWeight','bold','FontSize',10);
    end
 end

%export_fig epis_stem_5million.bmp
```

<div align="center"><img src="{{ "/images/Blog/EGWAS/manEpi.jpg" | prepend: site.baseurl }}"></div>

#### Coded for circle perl software (config file)

```conf
# config file for circle software

<<include etc/colors_fonts_patterns.conf>>
#<<include bands.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>
<image>
dir*=/Perl/new/circos-0.69-3/bin
file*=grainwe.png#rice7-7.png
#angle_orientation* = counterclockwise
#chromoso11mes_reverse = /chr1/
#image_map_use      = yes
#image_map_strict   = removeparam
background*= white
#svg = no
<<include etc/image.conf>>
</image>
#karyotype   = data/karyotype/chr1-chr12rice.txt
karyotype   = data/karyotype/bin1-bin12.txt
chromosomes_units = 1000000
#chromosomes       = chr1;chr2;chr3;chr4;chr5
#chromosomes_color   = chr1=red,chr2=orange,chr3=green,chr4=blue
#chromosomes_color = chr1=vvblack,chr2=gnbu-6-seq-4,chr3=orrd-5-seq-1,chr4=oranges-9-seq-7,chr5=pubu-4-seq-4
chromosomes_display_default = yes
<links>
<link>
file          = data/grainw1999.txt #interm/yieldperplant1998.a1.txt#link1.txt
radius        = 0.47r
#color         = red
#ribbon = yes
#flat   = yes
bezier_radius = 0.2r
#crest         = 0.2
thickness     = 15
</link>
<link>
file          = data/grainw1998.txt #interm/yieldperplant1998.a1.txt#link1.txt
radius        = 0.47r
#color         = red
#ribbon = yes
#flat   = yes
bezier_radius = 0.2r
#crest         = 0.2
thickness     = 15
</link>
<backgrounds>
<background>
color = vvlporange
</background>
</backgrounds>
</links>
<plots>
<plot>
type  = histogram
fill_under = yes
thickness = 4
file = data/grainw1998.heatmap.txt#bin1-bin12hist.txt
min = -3#-3
max = 2#4
r0 = 0.5r
r1 = 0.57r
<rules>
<rule>
condition=1
fill_color=eval(var(chr1))
color=eval(var(chr1)) #to make mutilple color
flow=continue
</rule>
</rules>
<backgrounds>
<background>
color = vlgrey
</background>
</backgrounds>
#url = script?type=label&value=[value]
</plot>
#</plot>
<plot>
type  = scatter
glyph = circle
glyph_size = 12p
#type = heatmap#histogram
file = data/grainw1998.scatter.txt#grainw1998.heatmap.txt#bin1-bin12hist.txt
min = 1#-3
max = 47.5#4
r0 = 0.6r
r1 = 0.67r
<rules>
<rule>
condition=1
fill_color=eval(var(chr1))
color=eval(var(chr1)) #to make mutilple color
flow=continue
</rule>
</rules>
<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>
#url = script?type=label&value=[value]
#url   = script?type=label&value=[value]&color=[color]
</plot>
<plot>
type  = text
color = blue
file  = data/rice_grainw_gene.txt
r0    = 0.678r
r1    = 0.85r
label_size = 30p
label_font = bold
show_links     = yes
link_dims      = 0p,10p,15p,10p,10p
link_thickness = 8p
link_color     = red
padding        = 0p
rpadding       = 0p
label_snuggle         = yes
max_snuggle_distance  = 0.8r
snuggle_tolerance     = 0.25r
snuggle_sampling      = 5
#url = script?type=label&id=[id]&text=[value]
<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>
</plot>
<plot>
type  = histogram
thickness = 4
fill_under = yes
fill_color = black
color = black
file = data/grainw1999.heatmap.txt#bin1-bin12hist.txt
min = -2.6#-3
max = 2.6#4
r0 = 0.86r
r1 = 0.92r
<rules>
<rule>
condition=1
fill_color=eval(var(chr1))
color=eval(var(chr1)) #to make mutilple color
flow=continue
</rule>
</rules>
<backgrounds>
<background>
color = vlgrey
</background>
</backgrounds>
#url = script?type=label&value=[value]
</plot>
<plot>
type  = scatter
fill_color = red
glyph = circle
glyph_size = 12p
file  = data/grainw1999.scatter.txt
r0    = 0.93r
r1    = 0.99r
min   = 0
max   = 43
<rules>
<rule>
condition=1
color=eval(var(chr1)) #to make mutilple color
flow=continue
</rule>
</rules>
<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>
#url = script?type=scatter-square&value=[value]&start=[start]&end=[end]
</plot>
</plots>
<<include etc/housekeeping.conf>>
```

<div align="center"><img src="{{ "/images/Blog/EGWAS/RICE.jpg" | prepend: site.baseurl }}"></div>



#### If you are interested in our model or any question and omission (bug). Please feel free to contact me.

* [Meng Luo](https://github.com/mengluoML)


