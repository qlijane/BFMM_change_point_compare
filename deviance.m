function dev=deviance(lam1,lam2,tauj,C,Nj,z,pi_v)% Compute the deviance of multiple individuals
% tauj is a vector, C is vector, z is a matrix with each row for a individual
% Nj is a vector recording the # of events for each driver, pi_v is a 38*2
% vector
 tau_s=unique(tauj);
 K=length(tau_s);% # of clusters
m=length(tauj);% # of drivers
idx=ones(1,m);
for j=2:K
idx(tauj==tau_s(j))=j;
end
dev=0;
for j=1:K
dev=dev-2*(log_L_multi(lam1(idx==j),lam2(idx==j),tau_s(j),C(idx==j),Nj(idx==j),z(idx==j,:),sum(idx==j),pi_v(idx==j,j)));
end