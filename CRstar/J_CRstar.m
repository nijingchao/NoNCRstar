%% CrossRankStar objective function value
function ObjValue = J_CRstar(Gnorm, Snorm, Ynorm, r, e, alpha, beta, c)

n = size(Gnorm,1);
In = speye(n);

AggrMat = c*(In - Gnorm) + alpha*(In - Snorm) + 2*beta*(In - Ynorm);
term1 = AggrMat*r;
term1 = r'*term1;

term2 = norm(r-e,'fro');
term2 = term2^2;
term2 = (1 - c)*term2;

ObjValue = term1 + term2;

end