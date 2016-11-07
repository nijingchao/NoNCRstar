%% WCRstar
function [r, alphas, phis_aux, Objs, Deltas] = WCRstar(Gnorms, Snorms, Ynorms, Ss, PhenotypeNeighbors, es, beta, gamma, c, MaxIter, epsilon)

%% Initialization
% Initialize ranking vectors
r = es;

% Initialize parameters
h = length(Gnorms);
ks = cellfun(@length,Gnorms);
ks = ks - 1;
alphas = cell(h,1); % alphas is a gx1 cell of vectors
phis_aux = cell(h,1); % phis_aux is a gx1 cell of vectors

for i = 1:h
    
    alphas{i} = (1/ks(i))*ones(ks(i),1);
    
end

omegas = 1 + (1./ks) + 2*beta;
omegas(ks == 0) = 1 + 2*beta;

% Convergence analysis parameters
J1 = cell(h,1);

for i = 1:h
    
    J1{i} = cell2mat(r{i});
    
end

J1_vec = cell2mat(J1);
J1_vec = (round(J1_vec*1e16))/1e16;

% J1 = 99999;
% J1 = (round(J1*1e16))/1e16;

delta = 99999;
Objs = [];
Deltas = [];
iter = 1;

%% Alternating minimization approach update loop
while delta > epsilon && iter <= MaxIter
    
    % Update r
    
    for i = 1:h
        
        omega_i = omegas(i);
        alphas_i = alphas{i};
        ks_i = ks(i);
        r_i = r{i};
        es_i = es{i};
        Snorms_i = Snorms{i};
        Ynorms_i = Ynorms{i};
        Gnorms_i = Gnorms{i};
        Ss_i = Ss{i};
        PhenotypeNeighbors_i = PhenotypeNeighbors{i};
        
        % Update star r
        
        Ssum = sparse(length(r_i{1}),1);
        
        for p = 1:ks_i
            
            Ssum = Ssum + (alphas_i(p)/omega_i)*(Snorms_i{p}*r_i{p+1});
            
        end
        
        Ysum = sparse(length(r_i{1}),1);
        
        for j = 1:length(PhenotypeNeighbors_i)
            
            r_j = r{PhenotypeNeighbors_i(j)};
            Ysum = Ysum + ((2*beta)/omega_i)*(Ynorms_i{j}*r_j{1});
            
        end
        
        r{i}{1} = (c/omega_i)*(Gnorms_i{1}*r_i{1}) + ((1-c)/omega_i)*es_i{1} + Ssum + Ysum;
    
        % Update auxiliary r
        
        r_i = r{i};
        r_i_1 = r_i{1};
        phis_i = zeros(ks_i,1);
        
        for p = 1:ks_i
            
            denom = 1 + alphas_i(p);
            r_i{p+1} = (c/denom)*(Gnorms_i{p+1}*r_i{p+1}) + ((1-c)/denom)*es_i{p+1} + (alphas_i(p)/denom)*((Snorms_i{p}')*r_i_1);
            
            % Compute auxiliary objectives
            phis_i(p) = J_aux(Ss_i{p}, r_i_1, r_i{p+1},  ks_i);
            
        end
        
        r{i} = r_i;
        
        % Update alphas
        
        [sort_phis_i, idx] = sort(phis_i, 'ascend');
        t = ks_i;
        lambda_i = 0;
        
        while t > 0
            
            lambda_i = (2*gamma + sum(sort_phis_i(1:t)))/t;
            
            if lambda_i - sort_phis_i(t) > 0
                break;
            else
                t = t - 1;
            end
            
        end
        
        for p = 1:ks_i
            
            if p <= t
                alphas_i(p) = (lambda_i - sort_phis_i(p))/(2*gamma);
            else
                alphas_i(p) = 0;
            end
            
        end
        
        alphas_i(idx) = alphas_i;
        alphas{i} = alphas_i;
        
        phis_aux{i} = phis_i;
        
    end
    
    
    % convergence analysis
    J2_vec = J1_vec;
    
    for i = 1:h
        
        J1{i} = cell2mat(r{i});
        
    end
    
    J1_vec = cell2mat(J1);
    J1_vec = (round(J1_vec*1e16))/1e16;
    
    delta = norm(J2_vec - J1_vec, 1);
    Objs = [Objs, norm(J1_vec, 1)];
    Deltas = [Deltas, delta];
    
%     J2 = J1;
%     J1 = J_WCRstar(A, Gnorms, Os, PhenotypeNeighbors, r, es, phis_aux, alphas, beta, gamma, c, h, ks);
%     J1 = (round(J1*1e16))/1e16;
%     
%     delta = J2 - J1;
%     Objs = [Objs, J1];
%     Deltas = [Deltas, delta];
    
    iter = iter + 1;
    
end

end