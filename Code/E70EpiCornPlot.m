function E70EpiCornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
% function E70EpiCornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
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
nodesize=90;
label={'A','L','R','P'};
matrix=zeros(11,length(label));
for i=11:-1:1
    for j=1:length(label)
        clear tlabel tpos
        tlabel=strcat(num2str(i),label(j));
        tpos=find(strcmp(sample,tlabel));
        if ~isempty(tpos)
            matrix(i,j)=tvalue(tpos);
        else
            matrix(i,j)=min(tvalue);
        end
    end
end
% caxis([min(min(matrix)),max(max(matrix))]);
weight=[0,0,0.2,0.35,0.4,0.47,0.51,0.52,0.53,0.54,0.55];
label={'A','L','R','P'};
% subplot(5,5,1)
for i=1:size(matrix,1)
    clear ydata
    if i==1
        scatter([1.8,3.2],[i i],nodesize,matrix(i,[1 4]),'o','fill');
        hold on
    else
%         ydata=linspace(1.2-weight(i),3.8+weight(i),4);
        scatter(linspace(1.2-weight(i),3.8+weight(i),4),[i i i i],nodesize,matrix(i,:),'o','fill')
        hold on
    end
end
box on
set(gca,'Xlim',[0,size(matrix,2)+1],'XTick',linspace(0.65,4.35,4),'XAxisLocation','top','XTickLabel',label);
set(gca,'Ylim',[0,size(matrix,1)+1],'YTick',1:size(matrix,1));
%  set(gca,'YTick',1:size(cvalue,1),'YTickLabel',gene);
%  xlabel('Region')
%  ylabel('Layer');
% set(gca,'FontSize',5);
% load('MyColormaps3','mycmapnew')
%load('MyColormaps.mat','mycmapnew');

% colorbar
if CLimset
    set(gca,'CLim',[CLimMin,CLimMax]);
else
    if max(max(matrix))==0
        set(gca,'CLim',[0,0.1]);
    else
        set(gca,'CLim',[0,max(max(matrix))]);
    end
end
%set(gcf,'Colormap',mycmapnew)
set(gcf,'Colormap',colormap)
set(gca,'position',[0.1,0.1,0.5,0.75]); 
set(gca,'FontSize',9,'FontName','Arial','LineWidth',1,'FontWeight','bold');
colorbar('FontSize',6);
title(['\fontsize{9}' name],'position',[2.6,-1.2]);
% set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition', [1 1 2.3 2.8]);  % originally used 
set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition', [1 1 3.5 2.7]); 
% % below code ara for the output format of .tif
% % setting the position of the figure in the output files and not display the figures
% set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition', [1 1 1.2 1]); %[1 1 0.65 1]
% % if not need to display the colorbar, set the below commend as: colorber off
% colorbar 
% set(gcf,'visible','off');
