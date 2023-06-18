
% genetage corn plot of interested genes for all samples from E5.5-E7.5 stages
clear
inputpath='./';
colormap(jet)
mycmapnew = colormap %(jet)
inputfilename='AverageExpSubs.txt';
file=importdata([inputpath inputfilename]);
gene=file.textdata(2:end,1); 
sample=file.textdata(1,2:end); 
value=log10(file.data+1);
% load([inputpath 'E5.5-E7.5.exprs.log10.mat']);
%load('MyColormaps.mat','mycmapnew');
mycmapnew = colormap(parula)
for i=1:length(sample)
    tmp=sample{i};
    pos=strfind(tmp,'.');
    nsample(i)={tmp(1:pos(1)-1)};
    stage(i)={tmp(pos(1)+1:end)};
end
ustage={'E5.5','E6.0','E6.5','E7.0','E7.5'};
for i=1:length(ustage)
    pos=find(strcmp(stage,ustage{i}));
    sv(i).sample=nsample(pos);
    sv(i).value=value(:,pos);
end
output=[pwd '/corn.plot/'];
if ~exist(output,'dir')
    mkdir(output);
end

for i=1:length(gene)
    FiveStagesCornPlot(sv,i,mycmapnew,gene{i},1);
    print(gcf,'-dpdf',[output gene{i} '.pdf']); close all;
end

% 
% % genetage corn plot of interested genes for only epiblast/ectoderm samples from E5.5-E7.5 stages
% clear
% inputpath='./';
% inputfilename='AverageExpSubs.txt';
% file=importdata([inputpath inputfilename]);
% gene=file.textdata(2:end,1); 
% sample=file.textdata(1,2:end); 
% %value=log10(file.data+1);
% value=log10(file.data+1);
% 
% % load([input 'E5.5-E7.5.exprs.log10.mat']);
% %load('MyColormaps.mat','mycmapnew');
% mycmapnew = colormap(parula)
% 
% %mycmapnew = colormap %(jet)
% 
% for i=1:length(sample)
%     tmp=sample{i};
%     pos=strfind(tmp,'.');
%     nsample(i)={tmp(1:pos(1)-1)};
%     stage(i)={tmp(pos(1)+1:end)};
% end
% ustage={'E5.5','E6.0','E6.5','E7.0','E7.5'};
% for i=1:length(ustage)
%     pos=find(strcmp(stage,ustage{i}));
%     sv(i).sample=nsample(pos);
%     sv(i).value=value(:,pos);
% end
% output=[pwd '/corn.plot/'];
% if ~exist(output,'dir')
%     mkdir(output);
% end
% for i=1:length(gene)
%     FiveStagesEpiCornPlot(sv,i,mycmapnew,gene{i},1);
%     print(gcf,'-dpdf',[output gene{i} '.Epi.pdf']); close all;
% end
% 
