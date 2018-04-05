% implement parametric clustering method for change-point detection
% find the percentage of correct K
clear
clc
tic
parpool(24)
%***********************************************************
%simulate data
u_v=[150 300];% two possible values of change-points
K_d=2;m=80;  % true value of lam1, lam2
B=200;% generate B datasets 
no_c=0;indi=zeros(1,B);
l1=0.25; l2=0.1; %kth row is for kth cluster,1st row is for before, 2ed row is for after
%-------
parfor s_no=1:B
no_c =no_c + FMM_para_f(u_v,m,l1,l2,K_d,@latent_simu_f_lamj);
end
result=100*no_c/B
toc
save untitled_1.mat result;