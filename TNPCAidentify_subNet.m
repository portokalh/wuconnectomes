load '/Users/alex/code/Matlab/popNet_HOPCA/data/response_all.mat'; %load response information
%disp(response_infor)
%load 'results/Dipy/181104*1gen3&4_TNPCA_V_332_1_K15_Outputs';
load '/Users/alex/code/Matlab/popNet_HOPCA/results/Dipy/190107*1genALL_TNPCA_V_332_1_K15_OutputsForgen';

selectForgen = [15;16;17;18;19;20;21;22;23;24;3;4;8;9;10;11;14;53;54;55];
response = response(selectForgen,:);
idx = response(response(:,2)~=0,2);
idx(idx==4)=1; %gen4
idx(idx==2)=0; %gen3

% selectForold = [7;26;28;29;30;33;34;46;47;51;36;37;38;39;40;41;42;43;44;45];
% response = response(selectForold,:);
% idx = response(response(:,2)==0,4);
% idx(idx==1)=0; %young
% idx(idx==2)=1; %old

data = cat(2,U(:,1:10),idx);
[w,t,fp]=fisher_training(data(:,1:end-1),data(:,end));
w = w/norm(w);
u0 = mean(data(data(:,end)==0,1:end-1));
u1 = mean(data(data(:,end)==1,1:end-1));
s = norm(u1-u0);

Net_change = zeros(332,332);
for i = 1:size(data(:,1:end-1),2)
    base_net = D(i)*w(i)*(V(:,i)*V(:,i)');
    Net_change = Net_change+base_net;
    %disp(i)
end
Net_change = s*Net_change;

%remove diagnal
Net_change = Net_change-diag(diag(Net_change));

k = 200;
maxcon = maxk(Net_change(:),k);
mincon = mink(Net_change(:),k);
mostcon = cat(1,maxcon,mincon);
mostcon = sort(mostcon,'descend','ComparisonMethod','abs');

%find the ROI index
ridx = [];
cidx = [];
for j = 1:2:k
    disp(mostcon(j))
    [r,c] = find(Net_change==mostcon(j));
    ridx = cat(1,ridx,r);
    cidx = cat(1,cidx,c);
end
mostcon_idx = [ridx cidx];
mostcon_idx = mostcon_idx(1:2:k,:);

Autonomy_data = readtable('/Users/alex/code/Matlab/popNet_HOPCA/data/atonomyInfo.csv');

%disp(' ')
%disp(' ')
%disp(['the top' num2str(k/2) ' connected subnetwork contribute to distinguish gen3 and gen4 '])
disp(['the top' num2str(k/2) ' connected subnetwork contribute to distinguish old and young '])
%disp(' ')

for h = 1:size(mostcon_idx,1)
    i1 = mostcon_idx(h,1);
    i2 = mostcon_idx(h,2);
    con1 = char(Autonomy_data(i1,2).Variables);
    d1 = char(Autonomy_data(i1,4).Variables);
    con2 = char(Autonomy_data(i2,2).Variables);
    d2 = char(Autonomy_data(i2,4).Variables);
    disp([num2str(h) ' ' con1 '_' d1 '---' con2 '_' d2 ' ' num2str(i1) ' ' num2str(i2)])
    %disp(' ')
end
    
    
