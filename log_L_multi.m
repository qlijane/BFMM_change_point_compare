function logr= log_L_multi(lam1,lam2,u,cj,Nj_r,t,K,pi_v)% Compute the log likelihood of multiple individuals
%in the same cluster, pi_v is K*1
% u is the component mean, cj is vector, t is a matrix with each row for a individual with the same change-point 
% Nj is a vector recording the # of events for each driver
% K is the # of drivers in this cluster
%tauj=normrnd(u,sigma,K,1);
Nj_1=zeros(K,1); Nj_2 = zeros(K,1);
for j=1:K
Nj_1(j)=sum(t(j,1:Nj_r(j))<=u);
Nj_2(j) = Nj_r(j)-Nj_1(j);
end
logr=-sum(lam1)*u+u*sum(lam2)-sum(lam2.*cj)+sum(Nj_1.*log(lam1))+sum(Nj_2.*log(lam2))+sum(log(pi_v));
%r=exp(logr);
