%% CRstar leave-one-out cross validation
function CRstar_CrossValidation(alpha, beta, c, MaxIter, epsilon, DataType)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for center-auxiliary network consistency.
% beta: a regularization parameter for cross-network consistency.
% c: a regularization parameter for query preference.
% MaxIter: the maximal number of iteration for updating ranking vector.
% epsilon: a convergence parameter.
% DataType: an indicator of which dataset to use, 1: P_G_NoSN_PPICenter.mat,
% 2: P_G_NoSN_GCNCenter.mat, 3: P_G_NoSN_PPICenter_2GCNs.mat.

%% Parameter initialization
if nargin < 6
    DataType = 1;
end
if nargin < 5
    epsilon = 1e-6;
end
if nargin < 4
    MaxIter = 1000;
end
if nargin < 3
    c = 0.85;
end
if nargin < 2
    beta = 0.5;
end
if nargin < 1
    alpha = 0.3;
end

%% Load NoSN data
if DataType == 1
    load('../Datasets/P_G_NoSN_PPICenter.mat');
elseif DataType == 2
    load('../Datasets/P_G_NoSN_GCNCenter.mat');
elseif DataType == 3
    load('../Datasets/P_G_NoSN_PPICenter_2GCNs.mat');
else
    error('The value of DataType is invalid');
end

%% CRstar precomputation
if DataType == 1
    PrecompFileName = 'CRstar_Precomp_Values_PPICenter.mat';
elseif DataType == 2
    PrecompFileName = 'CRstar_Precomp_Values_GCNCenter.mat';
else
    PrecompFileName = 'CRstar_Precomp_Values_PPICenter_2GCNs.mat';
end

if exist(PrecompFileName, 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('CRstar precomputation starts ...');
    CRstar_Precomputation(TSGeneNets, TSGeneNetsID, PhenotypeSimNet, PrecompFileName);
end

disp('Load the precomputation file ...');
load(PrecompFileName);

%% Extract test genes of the center networks
h = length(Seeds);
TestSeeds = cell(1,h);

for i = 1:h
    
    tmpSeeds = Seeds{i};
    TestSeeds{i} = tmpSeeds{1};
    
end

%% Leave-one-out cross validation
% Expand test genes s.t. test gene one by one
ExpandSeeds = vertcat(TestSeeds{:});

% Leave-one-out cross validation loop
RankScoreRecord = cell(1,length(ExpandSeeds));
RankRecord = cell(1,length(ExpandSeeds));
TotalCounter = 0;

disp('Leave-one-out cross validation starts ...');

for i = 1:h
    
    TSGeneNetsID_i = TSGeneNetsID{i};
    TestSeeds_i = TestSeeds{i};
    
    for j = 1:length(TestSeeds_i)
        
        TotalCounter = TotalCounter + 1;
        
        % Initialize query vector
        e = [];
        for p = 1:h
            
            % Initialize seeds for individual TPPINs
            tmpSeeds = Seeds{p};
            tmpTSGeneNetsID = TSGeneNetsID{p};
            
            for q = 1:length(tmpSeeds)
                
                if p == i && q == 1
                    head = length(e) + 1; % Record head position in e/r for final evaluation
                end
                
                tmp_e = zeros(length(tmpTSGeneNetsID{q}),1);
                
                if ismember(TestSeeds{i}(j),tmpSeeds{q})
                    subn = length(tmpSeeds{q}) - 1;
                    if subn ~= 0
                        [Fia, seedidx] = ismember(tmpSeeds{q},tmpTSGeneNetsID{q});
                        tmp_e(seedidx) = 1/subn;
                        [Fia1,seedidx1] = ismember(TestSeeds{i}(j),tmpTSGeneNetsID{q});
                        tmp_e(seedidx1) = 0;
                    end
                else
                    subn = length(tmpSeeds{q});
                    [Fia, seedidx] = ismember(tmpSeeds{q},tmpTSGeneNetsID{q});
                    tmp_e(seedidx) = 1/subn;
                end
                
                e = [e; tmp_e];
                
                if p == i && q == 1
                    tail = length(e); % Record tail position in e/r for final evaluation
                end
                
            end
            
        end
        
        % CRstar
        e = sparse(e);
        [r, Objs, Deltas] = CRstar(Gnorm, Snorm, Ynorm, e, alpha, beta, c, MaxIter, epsilon);
        
        % Record results
        RankScore = r(head:tail);
        RankScore = (round(RankScore*1e16))/1e16;
        FullRankScore = zeros(length(AllGeneID),1);
        [Fia, idx] = ismember(TSGeneNetsID_i{1}, AllGeneID);
        FullRankScore(idx) = RankScore;
        
        SeedGeneList = setdiff(TestSeeds{i}, TestSeeds{i}(j));
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

save('CRstarResults.mat','RankScoreRecord','RankRecord','ExpandSeeds','AllGeneID');

%% Evaluation
disp('AUC value evaluation starts ...');

AUCEvaluation(RankRecord, ExpandSeeds, AllGeneID);

end