function E75EpiCornPlot(tvaluen,samplename,genename,mycmapnew,CLimMin,CLimMax)
% function E75EpiCornPlot(tvalue,sample,name,mycmapnew,CLimMin,CLimMax)
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

label={'A','L1','R1','L2','R2','P'}; 
value=zeros(9,length(label));
for i=1:9
    for j=1:length(label)
        if i==2||i==9
            tlabel=strcat(num2str(i),label{j}(1));
            tpos=find(strcmp(samplename,tlabel));
            if ~isempty(tpos)
                value(i,j)=tvaluen(tpos);
            else
                value(i,j)=min(tvaluen);
            end
        else
            tlabel=strcat(num2str(i),label(j));
            tpos=find(strcmp(samplename,tlabel));
            if ~isempty(tpos)
                value(i,j)=tvaluen(tpos);
            else
                value(i,j)=min(tvaluen);
            end
        end
    end
end
% plot
matrix=value;
caxis([min(min(matrix)),max(max(matrix))]);
weight=[0,0,0.2,0.35,0.425,0.46,0.51,0.52,0.53,0.53];
nodesize=100;
for i=1:size(matrix,1)
    clear ydata
    if i==1
        scatter([1.8,3.2],[i i],nodesize,matrix(i,[1 6]),'o','fill');
        hold on
    elseif i==2||i==9
        scatter(linspace(1.2-weight(i),3.8+weight(i),4),[i i i i],nodesize,matrix(i,[1 2 5 6]),'o','fill');
        hold on
    else
        temppos=linspace(1.2-weight(i),3.8+weight(i),6);
        scatter(temppos([1,6]),[i i],nodesize,matrix(i,[1,6]),'o','fill'); 
        hold on;
        scatter(temppos([2:5]),[i i i i],nodesize*0.9,matrix(i,2:5),'o','fill');
        hold on
    end
end
box on
set(gca,'Xlim',[0,5],'XTick',[linspace(0.35,4.65,6)],'XAxisLocation','top','XTickLabel',label);
set(gca,'Ylim',[0,size(matrix,1)+1],'YTick',1:size(matrix,1));
if CLimset
    set(gca,'CLim',[CLimMin,CLimMax]);
else
    if max(max(matrix))==0
        set(gca,'CLim',[0,0.1]);
    else
        set(gca,'CLim',[0,max(max(matrix))]);
    end
end
set(gca,'FontSize',9,'FontName','Arial','LineWidth',1,'FontWeight','bold');
set(gcf,'Colormap',mycmapnew); 
colorbar('FontSize',6);
set(gca,'position',[0.1,0.1,0.55,0.8]);
title(['\fontsize{9}' genename],'position',[2.5,-1]);
% below code ara for the output format of .png
% setting the position of the figure in the output files and not display the figures
set(gcf, 'PaperPositionMode', 'manual','PaperUnits', 'inches','PaperPosition', [1 1 2.3 2.5]); %[1 1 0.65 1]
% print(gcf,'-dpng',[output lower(genename) '.png']);

