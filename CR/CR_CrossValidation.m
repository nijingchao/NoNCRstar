%% CR leave-one-out cross validation
function CR_CrossValidation(alpha, c, MaxIter, epsilon, DataType)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for cross-network consistency.
% c: a regularization parameter for seed preference.
% MaxIter: the maximal number of iteration for updating ranking vector.
% epsilon: a convergence parameter.
% DataType: an indicator of which dataset to use, 1: P_G_NoSN_PPICenter.mat,
% 2: P_G_NoSN_GCNCenter.mat.

%% Parameter initialization
if nargin < 5
    DataType = 1;
end
if nargin < 4
    epsilon = 1e-6;
end
if nargin < 3
    MaxIter = 1000;
end
if nargin < 2
    c = 0.85;
end
if nargin < 1
    alpha = 0.5;
end

%% Load NoSN data (CR only uses the center networks in the NoSN)
if DataType == 1
    load('../Datasets/P_G_NoSN_PPICenter.mat');
elseif DataType == 2
    load('../Datasets/P_G_NoSN_GCNCenter.mat');
else
    error('The value of DataType is invalid');
end

%% Extract NoN by considering center networks only
h = length(TSGeneNets);
G = cell(h,1);
G_ID = cell(h,1);
G_Seeds = cell(h,1);

for i = 1:h
    
    tmpTSGeneNets = TSGeneNets{i};
    tmpTSGeneNetsID = TSGeneNetsID{i};
    tmpSeeds = Seeds{i};
    G{i} = tmpTSGeneNets{1};
    G_ID{i} = tmpTSGeneNetsID{1};
    G_Seeds{i} = tmpSeeds{1};
    
end

%% CR precomputation
if DataType == 1
    PrecompFileName = 'CR_Precomp_Values_PPICenter.mat';
else
    PrecompFileName = 'CR_Precomp_Values_GCNCenter.mat';
end

if exist(PrecompFileName, 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('CR precomputation starts ...');
    CR_Precomputation(G, G_ID, PhenotypeSimNet, PrecompFileName);
end

disp('Load the precomputation file ...');
load(PrecompFileName);

%% Leave-one-out cross validation
% Expand test genes s.t. test gene one by one
ExpandSeeds = vertcat(G_Seeds{:});

% Leave-one-out cross validation loop
RankScoreRecord = cell(1,length(ExpandSeeds));
RankRecord = cell(1,length(ExpandSeeds));
TotalCounter = 0;

disp('Leave-one-out cross validation starts ...');

for j = 1:h
    
    for t = 1:length(G_Seeds{j})
        
        TotalCounter = TotalCounter + 1;
        
        % Initialize query vector
        e = [];
        
        for i = 1:h
            
            if i == j
                head = length(e) + 1; % Record head position in e/r for final evaluation
            end
            
            tmp_e = zeros(length(G_ID{i}),1);

            if ismember(G_Seeds{j}(t),G_Seeds{i})
                subn = length(G_Seeds{i}) - 1;
                if subn ~= 0
                    [Fia, seedidx] = ismember(G_Seeds{i},G_ID{i});
                    tmp_e(seedidx) = 1/subn;
                    [Fia1,seedidx1] = ismember(G_Seeds{j}(t),G_ID{i});
                    tmp_e(seedidx1) = 0;
                end
            else
                subn = length(G_Seeds{i});
                [Fia, seedidx] = ismember(G_Seeds{i},G_ID{i});
                tmp_e(seedidx) = 1/subn;
            end
            
            e = [e; tmp_e];
            
            if i == j
                tail = length(e); % Record tail position in e/r for final evaluation
            end
            
        end
        
        % CR
        e = sparse(e);
        [r, Objs, Deltas] = CR(Gnorm, Ynorm, I_n, e, alpha, c, MaxIter, epsilon);
        
        % Record results
        RankScore = r(head:tail);
        RankScore = (round(RankScore*1e16))/1e16;
        FullRankScore = zeros(length(AllGeneID),1);
        [Fia, idx] = ismember(G_ID{j}, AllGeneID);
        FullRankScore(idx) = RankScore;
        
        SeedGeneList = setdiff(G_Seeds{j}, G_Seeds{j}(t));
        [proj, ia, idx] = intersect(SeedGeneList, AllGeneID);
        FullRankScore(idx) = 0;
        
        RankScoreRecord{TotalCounter} = FullRankScore;
        
        % Sort results
        [SortFullRankScore, IX] = sort(FullRankScore, 'descend');
        RankRecord{TotalCounter} = IX;
        
        disp(['Finished Number of Folds/Total Number of Folds: ' num2str(TotalCounter) '/' num2str(length(ExpandSeeds))]);
        
    end
    
end

%% Save results
disp('Leave-one-out cross validation finishes, save results ...');

save('CRResults.mat','RankScoreRecord','RankRecord','ExpandSeeds','AllGeneID');

%% Evaluation
disp('AUC value evaluation starts ...');

AUCEvaluation(RankRecord, ExpandSeeds, AllGeneID);

end