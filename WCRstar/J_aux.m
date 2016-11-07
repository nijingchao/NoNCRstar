%% Weighted CrossRankStar center-auxiliary network inconsistency value
function Obj_aux = J_aux(S_istar_ip, r_istar, r_ip, ks_i)

D_istar = S_istar_ip*(S_istar_ip');
D_ip = (S_istar_ip')*S_istar_ip;

Obj_aux_1 = (1/(ks_i^(0.5)))*((S_istar_ip')*r_istar) - D_ip*r_ip;
Obj_aux_1 = norm(Obj_aux_1, 'fro');
Obj_aux_1 = Obj_aux_1^2;

Obj_aux_2 = (1/(ks_i^(0.5)))*(r_istar - D_istar*r_istar);
Obj_aux_2 = norm(Obj_aux_2, 'fro');
Obj_aux_2 = Obj_aux_2^2;

Obj_aux_3 = r_ip - D_ip*r_ip;
Obj_aux_3 = norm(Obj_aux_3, 'fro');
Obj_aux_3 = Obj_aux_3^2;

Obj_aux = Obj_aux_1 + Obj_aux_2 + Obj_aux_3;

end