%% WCRstar leave-one-out cross validation
function WCRstar_CrossValidation(beta, gamma, c, MaxIter, epsilon)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% beta: a regularization parameter for cross-network consistency.
% gamma: a parameter controlling the complexity of learnt weights alpha
% c: a regularization parameter for query preference.
% MaxIter: the maximal number of iteration for alternating minimization approach.
% epsilon: a convergence parameter.

%% Parameter initialization
if nargin < 5
    epsilon = 1e-3;
end
if nargin < 4
    MaxIter = 1000;
end
if nargin < 3
    c = 0.85;
end
if nargin < 2
    gamma = 0.005;
end
if nargin < 1
    beta = 0.5;
end

%% Load NoSN data (WCRstar uses P_G_NoSN_PPICenter_2GCNs.mat as default)
load('../Datasets/P_G_NoSN_PPICenter_2GCNs.mat');

%% WCRstar precomputation
if exist('WCRstar_Precomp_Values.mat', 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('WCRstar precomputation starts ...');
    WCRstar_Precomputation(TSGeneNets, TSGeneNetsID, PhenotypeSimNet);
end

disp('Load the precomputation file ...');
load('WCRstar_Precomp_Values.mat');

%% Extract test genes of the center networks
h = length(Seeds);
TestSeeds = cell(1,h);

for i = 1:h
    
    tmpSeeds = Seeds{i};
    TestSeeds{i} = tmpSeeds{1};
    
end

%% Leave-one-out cross validation
% Expand underlying seeds s.t. test seed one by one
ExpandSeeds = vertcat(TestSeeds{:});

% Leave-one-out cross validation loop
RankScoreRecord = cell(1,length(ExpandSeeds));
RankRecord = cell(1,length(ExpandSeeds));
TotalCounter = 0;
Allalphas = cell(1,length(ExpandSeeds));
Allphis_aux = cell(1,length(ExpandSeeds));

disp('Leave-one-out cross validation starts ...');

for i = 1:h
    
    TSGeneNetsID_i = TSGeneNetsID{i};
    TestSeeds_i = TestSeeds{i};
    
    for j = 1:length(TestSeeds_i)
        
        TotalCounter = TotalCounter + 1;
        
        % Construct query vectors
        es = cell(h,1);
        
        for p = 1:h
            
            % Initialize seeds for individual tissue-specific molecular networks
            tmpSeeds = Seeds{p};
            tmpTSGeneNetsID = TSGeneNetsID{p};
            tmp_es = cell(length(tmpSeeds), 1);
            
            for q = 1:length(tmpSeeds)
                
                tmp_e = sparse(length(tmpTSGeneNetsID{q}),1);
                
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
                
                tmp_es{q} = tmp_e;
                
            end
            
            es{p} = tmp_es;
        end
        
        % WCRstar
        [r, alphas, phis_aux, Objs, Deltas] = WCRstar(Gnorms, Snorms, Ynorms, Ss, PhenotypeNeighbors, es, beta, gamma, c, MaxIter, epsilon);
        
        % Record results
        RankScore = r{i}{1};
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
        
        % Record learnt weights and center-auxiliary network inconsistencies
        Allalphas{TotalCounter} = alphas;
        Allphis_aux{TotalCounter} = phis_aux;
        
        disp(['Finished Number of Folds/Total Number of Folds: ' num2str(TotalCounter) '/' num2str(length(ExpandSeeds))]);
        
    end
    
end

%% Save results
disp('Leave-one-out cross validation finishes, save results ...');

save('WCRStarResults.mat','RankScoreRecord','RankRecord','ExpandSeeds','AllGeneID','Allalphas','Allphis_aux');

%% Evaluation
disp('AUC value evaluation starts ...');

AUCEvaluation(RankRecord, ExpandSeeds, AllGeneID);

end