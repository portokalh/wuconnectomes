%-------------------------------Script Information-------------------------
% This script is used for extrcting Principle Components by Tensor PCA and
% also help you decide how many PCs to use
% Author: Wenlin Wu
% Last change:2018-11-4
% Change: add variables for different data(dipy or dsi studio)(2018-10-28)
% add gen_type(2018-11-4)

%-------------------------------Output Information-------------------------
%U: subject mode, store PC score for the data;
%V: network modem store network basis;
%percent_store: store variation information that top PCs can reflect

%-------------------------------Setup Parameter----------------------------
mypath = '/Users/alex/code/Matlab/popNet_HOPCA/data/';
outpath = '/Users/alex/code/Matlab/popNet_HOPCA/results/';

data_source = 1; %1-dipy, 2-DSIstudio
trial = '*1'; %today 1st try **important otherwise will overwrite other trials today
num_con = 332; %Please pick [332, 33]
type = 1; %Please pick [1-'All', 2-'noYoung', 3-'onlyYoung']
orth = 1; %Please pick [1-'V_orth', '2-U,V,W_orth']note:we usually orth = 1
gen_type = 1;%Please pick [1-'All', 2-'gen3&4', 3-'only gen0']

k = 15; % # of factors to extract

%----------------Main Part of TensorPCA&Variation identify-----------------
%load response to help filter data group
load([mypath 'response_all.mat']);

gen34_idx = response(response(:,2)~=0,1);
gen0_idx = response(response(:,2)==0,1);
young_idx = response(response(:,4)==1,1);

%load connectivity data
if data_source == 1
    outpath = [outpath 'Dipy/'];
    load([mypath 'connectivity_all' num2str(num_con) 'Dipy.mat'])
else
    outpath = [outpath 'DSI/'];
    load([mypath 'connectivity_all' num2str(num_con) 'DSI.mat'])
end

date = datestr(now,30);

if num_con == 332
    connectivity = connectivity332;
else
    connectivity = connectivity33;
end


%filter out the Young samples
if type == 2
    connectivity(:,:,young_idx) = [];
end

if type == 3
    connectivity = connectivity(:,:,young_idx);
end

%filter gennotype
if gen_type == 2
    connectivity = connectivity(:,:,gen34_idx);
    gen_outName = 'gen3&4';
elseif gen_type==3
    connectivity = connectivity(:,:,gen0_idx);
    gen_outName = 'gen0';
else
    gen_outName = 'genALL';
end

% selectForold = [7;26;28;29;30;33;34;46;47;51;36;37;38;39;40;41;42;43;44;45];
% connectivity = connectivity(:,:,selectForold);

% selectForgen = [15;16;17;18;19;20;21;22;23;24;3;4;8;9;10;11;14;53;54;55];
% connectivity = connectivity(:,:,selectForgen);

X = tensor(connectivity);

%PC_num = 5; %indicate number  of PCs you'd like to use

%Tensor PCA
%     rng(2333);
if orth == 1
    [V,D,U,W,Xhat,obj] = hopca_popNet(X,k,foptions);
    orth_type = 'V';
else
    [V,D,U,W,Xhat,obj] = hopca_popNet_new(X,k,foptions);
    orth_type = 'UVW';
end

%Calculate variation percentage
% U subject mode
%V network mode
PCs.U = U;
PCs.V = V;
percent_store = []; %help select # of PCs
for i = 1:k
    [percent] = var_explained(X,i,PCs);
    percent_store = cat(2, percent_store, [i;percent]);
    disp(['top ' num2str(i) '/' num2str(k) ' PCs explain ' num2str(percent) '% variance of X'])
end

save([outpath date(3:8) trial gen_outName '_TNPCA_' orth_type '_' num2str(num_con) '_' num2str(type) '_K' num2str(k) '_Outputs.mat'],'U','V','D','percent_store','obj','connectivity');
