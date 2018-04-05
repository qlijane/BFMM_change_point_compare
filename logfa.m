function log_fa=logfa(a0,b0,K,pi_v,a,m)
% log of P(a_0|pi)
pi_vp=pi_v(pi_v>0.000000001);
log_fa=(a0-1)*log(a)-b0*a+m*log(gamma(a))-m*K*log(gamma(a/K))+(a/K-1)*sum(sum(log(pi_vp)));