% implement parametric clustering method for change-point detection on
% simulated data set
% when the number of clusters are fixed, different lam for different
% clusters
clear
clc
rng('shuffle')
tic
parpool(24)
%***********************************************************
%simulate data
u_v=[150 300];% two possible values of change-points
K_d=2;m=40;  % true value of lam1, lam2
B=200;% generate B datasets 
l1=0.25; l2=0.1; %kth row is for kth cluster,1st row is for before, 2ed row is for after
%-------
estimate_all=zeros(K_d+3,B);% each row is a parameter, u1,u2,u3, l1, l2,  each column 
cover_p=zeros(K_d+2,1);% whether the CI covers the true paramete
true_v=[u_v l1 l2];
true_v1=[u_v 1000*[l1 l2]];
%-----------
parfor s_b=1:B
  [u_v_h,lam_1h,lam_2h,memship,cover_p_f]=FMM_MCMC_f(u_v,m,l1,l2,K_d,@latent_simu_f_lamj2);   
  estimate_all(:,s_b)=[u_v_h;1000*lam_1h;1000*lam_2h;100*memship/m];
  cover_p=cover_p+cover_p_f;
end
%***********************************************************
%find the smallest DIC
%***********************************************************
% %inference
RMSE=zeros(K_d + 2,1);Bias=zeros(K_d+2,1);
for i=1:K_d + 2
RMSE(i)= sqrt((estimate_all(i,:)-true_v1(i))*(estimate_all(i,:)-true_v1(i))'/B);
Bias(i)=sum(estimate_all(i,:)-true_v1(i))/(B*true_v1(i))*100;
end
result=[true_v1' mean(estimate_all(1:K_d+2,:),2) RMSE abs(Bias) 100*cover_p/B]
per_cor_cluster=mean(estimate_all(end,:))
save untitled2.mat result per_cor_cluster;
toc