%% CRstar precomputation
function CRstar_Precomputation(G, G_ID, A, PrecompFileName)

%%% Input parameters
%
% G: the tissue-specific molecular networks after mapping.
% G_ID: the corresponding gene IDs of molecular networks in G.
% A: the disease similarity network.
% PrecompFileName: the file name to store precomputation results.

%% Initialization
h = length(G);
ns = cell(h,1);
n = 0;

for i = 1:h
    
    ns{i} = cellfun(@length, G_ID{i});
    n = n + sum(ns{i});
    
end

%% Construct normalized G
Gnorms = cell(size(G));

for i = 1:h
    
    G_i = G{i};
    G_i = blkdiag(G_i{:});
    D = sum(G_i,2);
    D = D.^(-0.5);
    D = diag(D);
    G_i = D*G_i*D;
    Gnorms{i} = G_i;
    
end

Gnorm = blkdiag(Gnorms{:});

%% Construct normalized S
S = cell(h,1);
Du = cell(h,1);

for i = 1:h
    
    G_ID_i = G_ID{i};
    ns_i = ns{i};
    k_i = length(G_ID_i);
    S_i = cell(k_i,k_i);
    
    for p = 1:k_i
        
        for q = 1:k_i
            
            S_i{p,q} = sparse(ns_i(p),ns_i(q));
            
        end
        
    end
    
    Du_i = cell(k_i,1);
    
    if k_i-1 > 0
        Du_i{1} = (k_i-1)*speye(ns_i(1));
    else
        Du_i{1} = speye(ns_i(1));
        S_i{1,1} = speye(ns_i(1));
    end
    
    for p = 2:k_i
        
        [proj, I1, I2] = intersect(G_ID_i{1}, G_ID_i{p});
        S_i1_ij = sparse(I1,I2,ones(length(proj),1),ns_i(1),ns_i(p));
        S_i{1,p} = S_i1_ij;
        S_i{p,1} = S_i1_ij';
        Du_i{p} = speye(ns_i(p));
        
    end
    
    S_i = cell2mat(S_i);
    S{i} = S_i;
    Du_i = blkdiag(Du_i{:});
    Du{i} = Du_i;
    
end

S = blkdiag(S{:});
Du = blkdiag(Du{:});
Du = diag(Du);
Du = Du.^(-0.5);
Du = diag(Du);
Snorm = Du*S*Du;

%% Construct normalized Y
Y = sparse(n,n);
Dv = cell(h,1);
dA = sum(A,2);

for i = 1:h
    
    G_ID_i = G_ID{i};
    ns_i = ns{i};
    ns_sums = cellfun(@sum, ns);
    k_i = length(G_ID_i);
    
    Dv_i = cell(k_i,1);
    Dv_i{1} = dA(i)*speye(ns_i(1));
    
    for p = 2:k_i
        
        Dv_i{p} = speye(ns_i(p));
        
    end
    
    Dv_i = blkdiag(Dv_i{:});
    Dv{i} = Dv_i;
    
    Y_i = sparse(ns_sums(i),n);
    
    for j = i:h
        
        if i == j
            Y_ij = speye(ns_sums(i),ns_sums(j));
        else
            Y_ij = sparse(ns_sums(i),ns_sums(j));
        end
        
        G_ID_j = G_ID{j};
        ns_j = ns{j};
        [proj, I1, I2] = intersect(G_ID_i{1}, G_ID_j{1});
        O_ij = sparse(I1, I2, A(i,j)*ones(length(proj),1), ns_i(1), ns_j(1));
        Y_ij(1:size(O_ij,1),1:size(O_ij,2)) = O_ij;
        Y_i(:, sum(ns_sums(1:j-1))+1:sum(ns_sums(1:j))) = Y_ij;
        
    end
    
    Y(sum(ns_sums(1:i-1))+1:sum(ns_sums(1:i)), :) = Y_i;
    
end

Y = triu(Y) + triu(Y)' - diag(diag(Y));
Dv = blkdiag(Dv{:});
Dv = diag(Dv);
Dv = Dv.^(-0.5);
Dv = diag(Dv);
Ynorm = Dv*Y*Dv;

%% Save precomputed matrices
save(PrecompFileName,'Gnorm','Snorm','Ynorm');

end