%% WCRstar precomputation
function WCRstar_Precomputation(G, G_ID, A)

%%% Input parameters
%
% G: the tissue-specific molecular networks after mapping.
% G_ID: the corresponding gene IDs of molecular networks in G.
% A: the disease similarity network.

h = length(G);

%% Construct normalized G
Gnorms = cell(h,1);
for i = 1:h
    
    k_i = length(G{i});
    Gnorms{i} = cell(k_i,1);
    
    for p = 1:k_i
        
        D = sum(G{i}{p},2);
        D = D.^(-0.5);
        D = diag(D);
        Gnorms{i}{p} = D*G{i}{p}*D;
        
    end
    
end

%% Construct normalized S
Ss = cell(h,1);
Snorms = cell(h,1);

for i = 1:h
    
    tmpIDs = G_ID{i};
    ns = cellfun(@length,tmpIDs);
    k_i = length(tmpIDs) - 1;
    S_i = cell(k_i,1);
    Snorms_i = cell(k_i,1);
    
    for p = 1:k_i
        
        [proj, I1, I2] = intersect(tmpIDs{1}, tmpIDs{p+1});
        S_i1_ip = sparse(I1,I2,ones(length(proj),1),ns(1),ns(p+1));
        S_i{p} = S_i1_ip;
        Snorms_i{p} = (1/k_i^(0.5))*S_i1_ip;
        
    end
    
    Ss{i} = S_i;
    Snorms{i} = Snorms_i;
    
end

%% Construct normalized Y
Os = cell(h,1);
Ynorms = cell(h,1);
PhenotypeNeighbors = cell(h,1);
dA = sum(A,2);

for i = 1:h
    
    tmpIDs_i = G_ID{i};
    ns_i1 = length(tmpIDs_i{1});
    PhenotypeNeighbors_i = find(A(i,:)>0)';
    Ynorms_i = cell(length(PhenotypeNeighbors_i),1);
    Os_i = cell(length(PhenotypeNeighbors_i),1);
    
    for j = 1:length(PhenotypeNeighbors_i)
        
        t = PhenotypeNeighbors_i(j);
        tmpIDs_t = G_ID{t};
        ns_t1 = length(tmpIDs_t{1});
        [proj, I1, I2] = intersect(tmpIDs_i{1}, tmpIDs_t{1});
        O_it = sparse(I1, I2, ones(length(proj),1), ns_i1, ns_t1);
        Os_i{j} = O_it;
        Ynorms_i{j} = (A(i,t)/((dA(i)^(0.5))*(dA(t)^(0.5))))*O_it;
        
    end
    
    Os{i} = Os_i;
    Ynorms{i} = Ynorms_i;
    PhenotypeNeighbors{i} = PhenotypeNeighbors_i;
    
end

%% Save precomputed matrices
save('WCRstar_Precomp_Values.mat','Gnorms','Snorms','Ynorms','Ss','Os','PhenotypeNeighbors');

end