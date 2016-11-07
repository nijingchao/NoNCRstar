%% Weighted CrossRankStar objective function value
function ObjValue = J_WCRstar(A, Gnorms, Os, PhenotypeNeighbors, r, es, phis_aux, alphas, beta, gamma, c, h, ks)

ObjValue = 0;
dA = sum(A,2);

for i = 1:h
    
    ks_i = ks(i);
    alphas_i = alphas{i};
    phis_aux_i = phis_aux{i};
    r_i = r{i};
    es_i = es{i};
    Gnorms_i = Gnorms{i};
    Os_i = Os{i};
    PhenotypeNeighbors_i = PhenotypeNeighbors{i};
    
    ObjValue = ObjValue + J_within(Gnorms_i{1}, r_i{1}, es_i{1}, c);
    
    for p = 1:ks_i
        ObjValue = ObjValue + J_within(Gnorms_i{p+1}, r_i{p+1}, es_i{p+1}, c);
    end
    
    ObjValue = ObjValue + (alphas_i')*phis_aux_i;
    
    for j = 1:length(PhenotypeNeighbors_i)
        
        t = PhenotypeNeighbors_i(j);
        r_j = r{t};
        ObjValue = ObjValue + beta*A(i,t)*J_cross(Os_i{j}, r_i{1}, r_j{1}, dA(i), dA(t));
        
    end
    
    ObjValue = ObjValue + gamma*(norm(alphas_i, 'fro')^2);
end

end

%% J_within
function Obj_within = J_within(Gnorm_ip, r_ip, e_ip, c)

n_ip = size(Gnorm_ip,1);
I_ip = speye(n_ip);

term1 = (I_ip - Gnorm_ip)*r_ip;
term1 = (r_ip')*term1;

term2 = norm(r_ip - e_ip, 'fro');
term2 = term2^2;

Obj_within = c*term1 + (1-c)*term2;

end

%% J_cross
function Obj_cross = J_cross(O_ij, r_istar, r_jstar, d_i, d_j)

D_istar = O_ij*(O_ij');
D_jstar = (O_ij')*O_ij;

term1 = (1/(d_i^(0.5)))*((O_ij')*r_istar) - (1/(d_j^(0.5)))*(D_jstar*r_jstar);
term1 = norm(term1, 'fro');
term1 = term1^2;

term2 = (1/(d_i^(0.5)))*(r_istar - D_istar*r_istar);
term2 = norm(term2, 'fro');
term2 = term2^2;

term3 = (1/(d_j^(0.5)))*(r_jstar - D_jstar*r_jstar);
term3 = norm(term3, 'fro');
term3 = term3^2;

Obj_cross = term1 + term2 + term3;

end