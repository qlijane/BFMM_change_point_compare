function no_c = FMM_para_f(u_v,m,l1,l2,K_d,latent_simu_f_lamj)
[z,Nj,C]=latent_simu_f_lamj(u_v,m,l1,l2,K_d);
%________________________Gibbs start__________________________
% initialize
DIC=zeros(1,7);
for K=1:4;
[DIC(K),indi(K)]=FMM_f(z,Nj,C,K,m);
end
%***********************************************************
%find the smallest DIC
[min_DIC,k_f]=min(DIC(1:4));
no_c = 0;
if k_f==K_d 
    no_c=1;
end