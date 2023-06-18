function E75CornPlot(svalue,sample,genename,mycmapnew,CLimMin,CLimMax)
% function E75CornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
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

label={'EA','MA','A','L1','R1','L2','R2','P','MP','EP'}; 
makrt={'o','d','o','o','o','o','o','o','d','o'};
colorindex={'[0.2 0.2 0.2]','[0.2 0.2 0.2]','flat','flat','flat','flat','flat','flat','[0.2 0.2 0.2]','[0.2 0.2 0.2]'};
linewidth=[0.05,0.05,0.01,0.01,0.01,0.01,0.01,0.01,0.05,0.05];
matrix=zeros(10,length(label));
for i=1:10
    if i==1
        label={'EAP','MA','A','L1','R1','L2','R2','P','MP','EAP'}; 
    elseif i==2||i==9
        label={'EA','MA','A','L','R','L','R','P','MP','EP'}; 
    else
        label={'EA','MA','A','L1','R1','L2','R2','P','MP','EP'}; 
    end
    for j=1:length(label)
        tlabel=strcat(num2str(i),label(j));
        tpos=find(strcmp(sample,tlabel));
        if ~isempty(tpos)
            matrix(i,j)=svalue(tpos);
        else
            matrix(i,j)=NaN;
        end
    end
end
label={'EA','MA','A','L1','R1','L2','R2','P','MP','EP'};
weight=[0,0.08,0.2,0.35,0.425,0.46,0.51,0.52,0.53,0.535]*2.1;
nodesize=55;
for i=1:size(matrix,1)
    if i==1
        x=linspace(4,7,4);
        y=repmat(i,length(x));
        m=matrix(i,[1 3 8 10]);
        makr=makrt([1 3 8 10]);
        colors=colorindex([1 3 8 10]);
        lwidth=linewidth([1 3 8 10]);
    elseif i==2
        x=linspace(2.3,8.7,8);
        y=repmat(i,length(x));
        m=matrix(i,[1 2 3 4 7 8 9 10]);
        makr=makrt([1 2 3 4 7 8 9 10]);
        colors=colorindex([1 2 3 4 7 8 9 10]);
        lwidth=linewidth([1 2 3 4 7 8 9 10]);
    % elseif i==9
    %     x=linspace(1.1,9.9,8);
    %     y=repmat(i,length(x));
    %     m=matrix(i,[1 2 3 4 7 8 9 10]);
    %     makr=makrt([1 2 3 4 7 8 9 10]);
    %     colors=colorindex([1 2 3 4 7 8 9 10]);
    else
        x=linspace(2.1-weight(i),8.9+weight(i),10);
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
set(gca,'Ylim',[0.2,size(matrix,1)+0.8],'YTick',1:size(matrix,1));
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
colorbar('FontSize',8);
set(gca,'position',[0.1,0.45,0.17,0.28],'FontSize',8,'LineWidth',1);
title(['\fontsize{8}' genename],'position',[5.5,-1]);

