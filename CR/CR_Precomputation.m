%% CR precomputation
function CR_Precomputation(G, G_ID, A, PrecompFileName)

%%% Input parameters
%
% G: the tissue-specific molecular networks after mapping.
% G_ID: the corresponding gene IDs of molecular networks in G.
% A: the disease similarity network.
% PrecompFileName: the file name to store precomputation results.

%% Initialization
h = length(G);
ns = cellfun(@length, G_ID);
n = sum(ns);
I_n = speye(n);

%% Normalize G
Gnorms = cell(size(G));

for i = 1:h
    
    D = sum(G{i},2);
    D = D.^(-0.5);
    D = diag(D);
    Gnorms{i} = D*G{i}*D;
    
end

Gnorm = blkdiag(Gnorms{:});

%% Common node mapping
Y = sparse(n,n);
Dv = cell(h,1);
dA = sum(A,2);

for i = 1:h
    
    Dv{i} = dA(i)*speye(ns(i));
    Y_i = sparse(ns(i),n);
    
    for j = i:h
        
        [proj, I1, I2] = intersect(G_ID{i}, G_ID{j});
        Oij = sparse(I1, I2, ones(length(proj),1), ns(i), ns(j));
        Y_i(:, sum(ns(1:j-1))+1:sum(ns(1:j))) = A(i,j)*Oij;
        
    end
    
    Y(sum(ns(1:i-1))+1:sum(ns(1:i)), :) = Y_i;
    
end

Y = triu(Y) + triu(Y)' - diag(diag(Y));

%% Construct normalized Y
Dv = blkdiag(Dv{:});
Dv = diag(diag(Dv).^(-0.5));
Ynorm = Dv*Y*Dv;

%% Save precomputed matrices
save(PrecompFileName,'Gnorm','Ynorm','I_n');


end