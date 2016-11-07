%% CrossRank objective function value
function Obj = J_CR(Gnorm, Ynorm, I_n, r, e, alpha, c)

X = I_n - Ynorm;
Obj = r'*c*(I_n-Gnorm)*r + (1-c)*norm(r-e, 'fro')^2 + 2*alpha*(r'*X*r);

end