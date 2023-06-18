function E65CornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
% function E65CornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
%  This function is just suitable for the LMD-seq data
%  ----------------------------------------------------
%  Input: 
%     tvalue: expression value of other type of score fo each sample in this stage (the length must be the same as the sample)
%     sample: the name of all sample (cell array), the order much be consistent with the tvalue (eg. 1A, 2P...)
%     name: the name of this plot, usually the gene symbol
%     mycmapnew: set current colormap
%     CLimMin: minimum of the colorbar (optional)
%     CLimMax: maximum of the colorbar (optional)

% ---------------------------------------------------------------


if nargin==4
    CLimset=0; CLimMin=0; CLimMax=0;
else
    CLimset=1;
end

label={'EA','A','P','EP'};
makrt={'o','o','o','o'};
colorindex={'k','flat','flat','k'};
colorindex={'[0.2 0.2 0.2]','flat','flat','[0.2 0.2 0.2]'};
linewidth=[0.05,0.01,0.01,0.05];

matrix=zeros(8,length(label));
for i=1:8
    if i==1
        label={'EAP','AP','AP','EAP'}; 
    else
        label={'EA','A','P','EP'}; 
    end
    for j=1:length(label)
        tlabel=strcat(num2str(i),label(j));
        tpos=find(strcmp(sample,tlabel));
        if ~isempty(tpos)
            matrix(i,j)=tvalue(tpos);
        else
            matrix(i,j)=NaN;
        end
    end
end
label={'EA','A','P','EP'}; 
weight=[0,0,0.2,0.35,0.4,0.46,0.51,0.52,0.53];

nodesize=55;
for i=1:size(matrix,1)
    if i==1
        x=linspace(1.8,3.2,3);
        y=repmat(i,length(x));
        m=matrix(i,[1 2 4]);
        makr=makrt([1 2 4]);
        colors=colorindex([1 2 4]);
        lwidth=linewidth([1 2 4]);
    else
        x=linspace(1.2-weight(i),3.8+weight(i),4);
        y=repmat(i,length(x));
        m=matrix(i,:);
        makr=makrt;
        colors=colorindex;
        lwidth=linewidth;
    end
    for j=1:length(m)
        if ~isnan(m(j))
            scatter(x(j),y(j),nodesize,m(j),makr{j},'fill','MarkerEdgeColor',colors{j},'LineWidth',lwidth(j)); hold on
        else
            scatter(x(j),y(j),nodesize,0,makr{j},'MarkerEdgeColor',colors{j},'LineWidth',lwidth(j)); hold on
        end
    end
end

box on
set(gca,'Xlim',[0,size(matrix,2)+1],'XTick',1:size(matrix,2),'XAxisLocation','top','XTickLabel',label,'XTickLabelRotation',90);
set(gca,'Ylim',[0,size(matrix,1)+1],'YTick',1:size(matrix,1));

if CLimset
    set(gca,'CLim',[CLimMin,CLimMax]);
else
    if isempty(tvalue)
        set(gca,'CLim',[0,0.1]);
    else
        if max(matrix(~isinf(matrix)))==0
            set(gca,'CLim',[0,0.1]);
        else
            set(gca,'CLim',[0,max(tvalue)]);
        end
    end
end
set(gcf,'Colormap',mycmapnew); 
set(gca,'position',[0.1,0.6,0.15,0.185],'FontSize',8,'LineWidth',1); 
colorbar('FontSize',8);
title(['\fontsize{8}' name],'position',[2.5,-1.3]);
