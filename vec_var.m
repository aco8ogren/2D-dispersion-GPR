V_temp = vfcn(wv_e,'gridded');
t1_comp = reshape(V_temp',[],1);
t2_comp = kfcn(wv_e,wv_s,'scattered');
t3_comp = kfcn(wv_s,wv_s,'scattered');
t4_comp = permute(t2_comp,[2 1]);

t2_comp = permute(t2_comp,[3 2 1]);

size(t1_comp)
size(t2_comp)
size(t3_comp\t4_comp)
t5_comp = t3_comp\t4_comp;
t5_comp = permute(t5_comp,[1 3 2]);
size(t5_comp)
% size(t2_comp.*(t3_comp\t4_comp))
% size(t2_comp.*t5_comp)
size(pagemtimes(t2_comp,t5_comp))

t_subtract_comp = permute(pagemtimes(t2_comp,t5_comp),[3 2 1]);
variances_comp = t1_comp - t_subtract_comp;
V_comp = reshape(variances_comp,N_evaluate(2),N_evaluate(1));

