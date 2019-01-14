%-------------------------------Script Information-------------------------
% This script is used for ploting top3 Principle Components by Tensor PCA
% and trying to find clustering information
% Author: Wenlin Wu
% Last change:2018-08-20
% Change: add savefig

%-------------------------------Output Information-------------------------
%Please check all output picture for detailed information

%-------------------------------Setup Parameter----------------------------
mypath = '/Users/alex/code/Matlab/popNet_HOPCA/results/Dipy/';
outpath = '/Users/alex/code/Matlab/popNet_HOPCA/results/Dipy/';
response_path = '/Users/alex/code/Matlab/popNet_HOPCA/data/';

%---------PC score infor---------
%input information to describe which PC you'd like to plot
num_con = 332; %Please pick [332, 33]
trial_res = '*1'; %input which trial of results you'd like to check
date_res = '190107'; % input which date of results you'd like to check

type_res = 1; %Please pick [1-'All', 2-'noYoung', 3-'onlyYoung'] 
orth = 1; %Please pick [1-'V_orth', '2-U,V,W_orth']note:we usually orth = 1
k = 15; % # of factors to extract

%-----visulization parameter-----
%indicate which group(s) of genotype you want to play with
gen_type = 1;%1-gen0, 2-gen3, 3-gen4, 4-gen3&gen4, 5-ALl animal
age_type = 3;%1-all, 2-old, 3-young

%------------------------Main Part of Plot PC scores-----------------------
if orth == 1
    orth_type = 'V';
else
    orth_type = 'UVW';
end

%load PC score information
% load([mypath date_res trial_res '_TNPCA_' orth_type '_' num2str(num_con) '_'...
%     num2str(type_res) '_K' num2str(k) '_OutputsNor.mat'], 'U', 'V', 'percent_store')
load([mypath date_res trial_res 'genALL_TNPCA_' orth_type '_' num2str(num_con) '_'...
    num2str(type_res) '_K' num2str(k) '_Outputs.mat'], 'U', 'V', 'percent_store')
%load Animal information e.g. gender, genotype
%col1-index; col2-genotype(revised:3-2); col3-sex(0-male,1-female);
%col4-age(revised old-2/young-1);
load([response_path 'response_all.mat'])

if type_res == 2
    response = response(response(:,4)==2,:);
elseif type_res == 3
    response = response(response(:,4)==1,:);
end

x = U(:,1:3);
y = response(:,1:4);

name_gen0 = 'all';
if age_type == 2
    temp = y(:,4)==2;
    y = y(temp,:);
    x = x(temp,:);
    name_gen0 = 'old';
    name_ = 'old';
elseif age_type == 3
    temp = y(:,4)==1;
    y = y(temp,:);
    x = x(temp,:);
    name_gen0 = 'young';
    name_ = 'young';
end

%divide data by different gennotype
idx_gen0 = find(y(:,2)==0);
idx_gen3 = find(y(:,2)==2);
idx_gen4 = find(y(:,2)==4);

if gen_type == 1
    idx = idx_gen0;
    name = ['gen0' name_gen0];
elseif gen_type == 2
    idx = idx_gen3;
    name = 'gen3';
elseif gen_type == 3
    idx = idx_gen4;
    name = 'gen4';
elseif gen_type == 4
    idx = cat(1,idx_gen3,idx_gen4);
    name = 'gen3&4';
else
    idx = [1:size(x,1)].';
    name = ['all' name_];
end

data_x = x(idx,:);
data_y = y(idx,:);

%-------plot--------
%cmap = [0 1 0; 1 0 0; 0 0 1];
cmapsex = [0/255 204/255 255/255;255/255 153/255 255/255];
cmapgen = [0 1 0; 0 0 1; 1 0 0];
cmapage = [0 1 0; 1 0 1];

%----------sex-----------
%plot 3D scatter points colored by sex
f1 = figure(1);
scatter3(data_x(:,1),data_x(:,2),data_x(:,3),[],data_y(:,3),'*')
colormap(cmapsex);
title(['Top 3 PC scores for ' name ' animals by sex'],'FontSize',...
    16, 'FontWeight', 'bold');
xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
f1name = [name '_Top3_PCs' '_colorbysex'];
colorbar
saveas(f1, [outpath f1name],'png');
savefig([outpath f1name]);

%plot 2D scatter points colored by sex
f2 = figure(2);
scatter(data_x(:,1),data_x(:,2),[],data_y(:,3),'*')
colormap(cmapsex);
title(['PC1&PC2 for ' name ' animals by sex'],'FontSize',...
    16, 'FontWeight', 'bold');
xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
f2name = [name '_PC1&PC2'  '_colorbysex'];
colorbar
saveas(f2, [outpath f2name],'png');
savefig([outpath f2name]);

f3 = figure(3);
scatter(data_x(:,2),data_x(:,3),[],data_y(:,3),'*')
colormap(cmapsex);
title(['PC2&PC3 for ' name ' animals by sex'],'FontSize',...
    16, 'FontWeight', 'bold');
xlabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
f3name = [name '_PC2&PC3'  '_colorbysex'];
colorbar
saveas(f3, [outpath f3name],'png');
savefig([outpath f3name]);

f4 = figure(4);
scatter(data_x(:,1),data_x(:,3),[],data_y(:,3),'*')
colormap(cmapsex);
title(['PC1&PC3 for ' name ' animals by sex'],'FontSize',...
    16, 'FontWeight', 'bold');
xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
f4name = [name '_PC1&PC3'  '_colorbysex'];
colorbar
saveas(f4, [outpath f4name],'png');
savefig([outpath f4name]);

%plot age effect in gen0
if gen_type == 1
    
    %plot 3D scatter points colored by age in gen0
    f5 = figure(5);
    scatter3(data_x(:,1),data_x(:,2),data_x(:,3),[],data_y(:,4),'*')
    colormap(cmapage);
    title(['Top 3 PC scores for ' name ' animals by age'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    zlabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f5name = [name '_Top3_PCs'  '_colorbyage'];
    colorbar
    saveas(f5, [outpath f5name],'png');
    savefig([outpath f5name]);

    %plot 2D scatter points colored by age in gen0
    f6 = figure(6);
    scatter(data_x(:,1),data_x(:,2),[],data_y(:,4),'*')
    colormap(cmapage);
    title(['PC1&PC2 for ' name ' animals by age'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    f6name = [name '_PC1&PC2'  '_colorbyage'];
    colorbar
    saveas(f6, [outpath f6name],'png');
    savefig([outpath f6name]);

    f7 = figure(7);
    scatter(data_x(:,2),data_x(:,3),[],data_y(:,4),'*')
    colormap(cmapage);
    title(['PC2&PC3 for ' name ' animals by age'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f7name = [name '_PC2&PC3'  '_colorbyage'];
    colorbar
    saveas(f7, [outpath f7name],'png');
    savefig([outpath f7name]);

    f8 = figure(8);
    scatter(data_x(:,1),data_x(:,3),[],data_y(:,4),'*')
    colormap(cmapage);
    title(['PC1&PC3 for ' name ' animals by age'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f8name = [name '_PC1&PC3'  '_colorbyage'];
    colorbar
    saveas(f8, [outpath f8name],'png');
    savefig([outpath f8name]);
end

%plot genotype effect
if gen_type == 4 || gen_type==5

    %plot 3D scatter points colored by gennotype in gen3&4
    f9 = figure(9);
    scatter3(data_x(:,1),data_x(:,2),data_x(:,3),[],data_y(:,2),'*')
    colormap(cmapgen);
    title(['Top 3 PC scores for ' name ' animals by gennotype'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    zlabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f9name = [name '_Top3_PCs'  '_colorbygenno'];
    colorbar
    saveas(f9, [outpath f9name],'png');
    savefig([outpath f9name]);

    %plot 2D scatter points colored by gennotype
    f10 = figure(10);
    scatter(data_x(:,1),data_x(:,2),[],data_y(:,2),'*')
    colormap(cmapgen);
    title(['PC1&PC2 for ' name ' animals by gennotype'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    f10name = [name '_PC1&PC2'  '_colorbygenno'];
    colorbar
    saveas(f10, [outpath f10name],'png');
    savefig([outpath f10name]);

   f11 = figure(11);
    scatter(data_x(:,2),data_x(:,3),[],data_y(:,2),'*')
    colormap(cmapgen);
    title(['PC2&PC3 for ' name ' animals by gennotype'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC2', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f11name = [name '_PC2&PC3'  '_colorbygenno'];
    colorbar
    saveas(f11, [outpath f11name],'png');
    savefig([outpath f11name]);
    
    f12 = figure(12);
    scatter(data_x(:,1),data_x(:,3),[],data_y(:,2),'*')
    colormap(cmapgen);
    title(['PC1&PC3 for ' name ' animals by gennotype'],'FontSize',...
        16, 'FontWeight', 'bold');
    xlabel('PC1', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('PC3', 'FontSize', 16, 'FontWeight', 'bold');
    f12name = [name '_PC1&PC3'  '_colorbygenno'];
    colorbar
    saveas(f12, [outpath f12name],'png');
    savefig([outpath f12name]);
end

    
