function E60EpiCornPlot(tvalue,sample,genename,mycmapnew,CLimMin,CLimMax)
% function E60EpiCornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
%  This function is just suitable for the LMD-seq data
%  ----------------------------------------------------
%  Input: 
%     tvalue: expression value of other type of score fo each sample in this stage (the length must be the same as the sample)
%     sample: the name of all sample (cell array), the order much be consistent with the tvalue (eg. 1A, 2P...)
%     genename: the name of this plot, usually the gene symbol
%     mycmapnew: set current colormap
%     CLimMin: minimum of the colorbar (optional)
%     CLimMax: maximum of the colorbar (optional)

% ---------------------------------------------------------------

if nargin==4
    CLimset=0; CLimMin=0; CLimMax=0;
else
    CLimset=1;
end

label={'A','P'};
for i=1:5
    for j=1:length(label)
        if i==1
            tpos=find(strcmp(sample,strcat(num2str(i),char(label)')));
            if ~isempty(tpos)
                matrix1(i,j)=tvalue(tpos);
            else
                matrix1(i,j)=min(tvalue);
            end
        else
            tlabel=strcat(num2str(i),label(j));
            tpos=find(strcmp(sample,tlabel));
            if ~isempty(tpos)
                matrix1(i,j)=tvalue(tpos);
            else
                matrix1(i,j)=min(tvalue);
            end
        end
    end
end
caxis([min(min(matrix1)),max(max(matrix1))]);
weight=[0,0.24,0.36,0.44,0.46];
nodesize=90;
for i=1:size(matrix1,1)
    clear ydata
    if i==1
        scatter(1.5,i,nodesize,matrix1(i,1),'o','fill');
        hold on;
    else
        scatter([1-weight(i),2+weight(i)],[i i],nodesize,matrix1(i,[1,2]),'o','fill'); 
        hold on;

    end
end
box on;
if CLimset
    set(gca,'CLim',[CLimMin,CLimMax]);
else
    if max(max(matrix1))==0
        set(gca,'CLim',[0,0.1]);
    else
        set(gca,'CLim',[0,max(max(matrix1))]);
    end
end
set(gca,'Xlim',[0,3],'XTick',[0.5,2.5],'XAxisLocation','top','XTickLabel',label);
set(gca,'Ylim',[0,size(matrix1,1)+1],'YTick',1:size(matrix1,1));
set(gca,'FontSize',9,'FontName','Arial','LineWidth',1,'FontWeight','bold');
set(gcf,'Colormap',mycmapnew); 
colorbar('FontSize',6);
title(['\fontsize{9}' genename],'position',[1.5,-1.2]);

set(gca,'position',[0.1,0.1,0.32,0.535]);
set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition', [1 1 2 2.4]); 


