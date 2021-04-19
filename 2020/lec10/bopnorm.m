function value = bopnorm(A,tail_es,nu)% the operator norm of bounded operators with tail
value = max(wnorm_mat(A,nu),tail_es);
