function FiveStagesCornPlotToshi(sv,m,mycmap,gene,setzeros,min1,max1,min2,max2,min3,max3,min4,max4,min5,max5)
% function FiveStagesCornPlot(sv,m,mycmap,gene,setzeros,min1,max1,min2,max2,min3,max3,min4,max4,min5,max5)
%  This function is designed to generate the corn plot for E5.5-E7.5 embryos.
%  ----------------------------------------------------
%  Input: 
%     sv: a struct which contains sample information and the corresponding values for E5.5, E6.0, E6.5, E7.0 and 
%         E7.5 embryos (sv.sample and sv.value, sv(1).sample means samples from E5.5 stage, and sv(1).value means 
%         values from E5.5 stage). The value for each stage could a one line vector or a matrix. Usually the 
%         sv.value contains the gene expression levels. Each line in sv.value across all five stages should be consistent.
%     m: which line of sv.value should be plot, if sv.value only contain one line, m=1.
%     gene: the name of this plot, usually the gene symbol.
%     setzeros: whether set the minimum of color bar of all stages as 0, '1' means set the minimum of color bar as 0 and '0' means not.
%     min1,max1: minimum and maxinumof the colorbar for E5.5 (optional).
%     min2,max2: minimum and maxinumof the colorbar for E6.0 (optional).
%     min3,max3: minimum and maxinumof the colorbar for E6.5 (optional).
%     min4,max4: minimum and maxinumof the colorbar for E7.0 (optional).
%     min5,max5: minimum and maxinumof the colorbar for E7.5 (optional).

% ---------------------------------------------------------------

if nargin==5
    if setzeros
        if max(sv(1).value(m,:))==0
            min1=0; max1=0.01;
        else
            min1=0; max1=max(sv(1).value(m,:));
        end
        if max(sv(2).value(m,:))==0
            min2=0; max2=0.01;
        else
            min2=0; max2=max(sv(2).value(m,:));
        end
        if max(sv(3).value(m,:))==0
            min3=0; max3=0.01;
        else
            min3=0; max3=max(sv(3).value(m,:));
        end
        if max(sv(4).value(m,:))==0
            min4=0; max4=0.01;
        else
            min4=0; max4=max(sv(4).value(m,:));
        end
        if max(sv(5).value(m,:))==0
            min5=0; max5=0.01;
        else
            min5=0; max5=max(sv(5).value(m,:));
        end
    else
        min1=min(sv(1).value(m,:)); max1=max(sv(1).value(m,:));
        min2=min(sv(2).value(m,:)); max2=max(sv(2).value(m,:));
        min3=min(sv(3).value(m,:)); max3=max(sv(3).value(m,:));
        min4=min(sv(4).value(m,:)); max4=max(sv(4).value(m,:));
        min5=min(sv(5).value(m,:)); max5=max(sv(5).value(m,:));
    end
elseif nargin==15
    if setzeros
        min1=0; min2=0; min3=0; min4=0; min5=0;
    end
end

min1=0.65
max1=0.98
min2=0.65
max2=0.98
min3=0.65
max3=0.98
min4=0.65
max4=0.98
min5=0.65
max5=0.98


min1=-0.05
max1=0.05
min2=-0.05
max2=0.05
min3=0-0.05
max3=0.05
min4=-0.05
max4=0.05
min5=-0.05
max5=0.05

sample1=sv(1).sample;
tvalue1=sv(1).value(m,:);
sample2=sv(2).sample;
tvalue2=sv(2).value(m,:);
sample3=sv(3).sample;
tvalue3=sv(3).value(m,:);
sample4=sv(4).sample;
tvalue4=sv(4).value(m,:);
sample5=sv(5).sample;
tvalue5=sv(5).value(m,:);
figure;
subplot(1,5,1);
if isempty(tvalue1)
    min1=0; max1=0.1;
end
E55CornPlot(tvalue1,sample1,'E5.5',mycmap,min1,max1)
set(gca,'position',[0.05,0.25,0.06,0.2]);
subplot(1,5,2);
if isempty(tvalue2)
    min2=0; max2=0.1;
end
E60CornPlot(tvalue2,sample2,'E6.0',mycmap,min2,max2);
set(gca,'position',[0.2 0.25 0.06 0.3]);
subplot(1,5,3);
if isempty(tvalue3)
    min3=0; max3=0.1;
end
E65CornPlot(tvalue3,sample3,'E6.5',mycmap,min3,max3);
set(gca,'position',[0.35 0.25 0.08 0.4]); 
subplot(1,5,4);
if isempty(tvalue4)
    min4=0; max4=0.1;
end
E70CornPlot(tvalue4,sample4,'E7.0',mycmap,min4,max4);
set(gca,'position',[0.52 0.25 0.13 0.53]); 
subplot(1,5,5);
if isempty(tvalue5)
    min5=0; max5=0.1;
end
E75CornPlot(tvalue5,sample5,'E7.5',mycmap,min5,max5);
set(gca,'position',[0.75 0.25 0.14 0.54]); 



set(gcf,'position',[200,200,750,300],'PaperPositionMode','manual','PaperUnits','inches','PaperPosition', [0.1 0.1 7.5 2.8]);
H=suptitle(gene);

h.Interpreter = 'none';