function [u_v_h,lam_1h,lam_2h,memship,cover_p]=FMM_MCMC_f(u_v,m,l1,l2,K_d,latent_simu_f_lamj)
Bt=20000;%# of samples for Gibbs
Burn_in=10000;% burn in time
lag=5;%every 5th result is used
true_v=[u_v  [l1 l2]];
cover_p=zeros(K_d+2,1);% whether the CI covers the true paramete
a=[10 10 ]; b=[430 430 ];
b1=1;b2=1;%prior of lam_b,lam_a
 a1=0.25;a2=0.1;
alpha0=1.8;
[z,Nj,C,tau_index]=latent_simu_f_lamj(u_v,m,l1,l2,K_d);%___________data simulation end_______________
 %________________________Gibbs start__________________________
lam_1_gib=zeros(m,Bt);lam_2_gib=zeros(m,Bt);
lam_1_gib(:,1)=0.05; lam_2_gib(:,1)=0.1; 
 % # of uniquce values for model
pi_v_gib=zeros(m,K_d,Bt);% start from uniform, mixture proportion
u_v_gib=zeros(K_d,Bt);% unique values
u_v_gib(:,1)=[100 350];
%---------------------------
% initialize randomly
zj=unidrnd(K_d,m,1);tau_gib=u_v_gib(1,1)*ones(m,Bt);
% zj is the index of tau_j in uk;
i_gib=zeros(m,K_d,Bt);
for k_t=2:K_d
   tau_gib(zj==k_t,1)=u_v_gib(k_t,1);% some drivers change at u1, others change at u2;
end
    for j=1:m
i_gib(j,zj(j),1)=1;
    end
%-------------------------
%_________initialize end__________________
%___________
    for t_g=2:Bt
     t_g
        % a_gib=alpha0(t_g-1)/K+i_gib(:,:,t_g-1);
       a_gib=alpha0/K_d+i_gib(:,:,t_g-1);
        for j=1:m
        pi_v_gib(j,:,t_g)=drchrnd(a_gib(j,:),1);%update pi
        end
% update u

           for k_u=1:K_d
           % k_u=1;
             lam_1=lam_1_gib(zj==k_u,t_g-1);lam_2=lam_2_gib(zj==k_u,t_g-1);
             k_n=sum(zj==k_u);% # of drivers in the kth component
             C_u=C(zj==k_u);Nj_u=Nj(zj==k_u);z_u=z(zj==k_u,:);  pi_u= pi_v_gib(zj==k_u,k_u,t_g) ;                          
         % use rejection accptance sampling from P(uk|data)
           u_v_gib(k_u,t_g)=ra_sampl_u(a(k_u), b(k_u),lam_1,lam_2,C_u,Nj_u,z_u,k_n,pi_u);
           end 
           u_v_gib(:,t_g)=sort(u_v_gib(:,t_g));
% update tau
          for j=1:m             
               log_L=zeros(1,K_d);
               for k=1:K_d
               log_L(k)=log_L_multi(lam_1,lam_2,u_v_gib(k,t_g),C(j),Nj(j),z(j,:),1,pi_v_gib(j,k_u,t_g));
               end 
               adj_c=mean(log_L);
               L_adj=exp(log_L-adj_c);
               pi_Lj=pi_v_gib(j,:,t_g).*L_adj;
               rand_p=sum(pi_Lj)*rand;
               idx = sum(rand_p>cumsum(pi_Lj)) + 1;
               tau_gib(j,t_g)=u_v_gib(idx,t_g);      
               zj(j)=idx;
           end
% update ik
    for j=1:m
i_gib(j,zj(j),t_g)=1;
    end
% update lam_1
           tau =  tau_gib(:,t_g);
            for j=1:m
                Nj_1=sum(z(j,1:Nj(j))<=tau(j));            
            lam_1_gib(j,t_g)=gamrnd(a1+Nj_1,1/(b1+tau(j)),1);
% update lam_2
            lam_2_gib(j,t_g)=gamrnd(a2+Nj(j)-Nj_1,1/(b2+C(j)-tau(j)),1);
            end
    end
%________________________Gibbs end__________________________
%compute DIC
lam_1_use=lam_1_gib(:,Burn_in:lag:end);
  lam_2_use=lam_2_gib(:,Burn_in:lag:end);
  lam_1_use_v=lam_1_use(:);
  lam_2_use_v=lam_2_use(:);
 lam_1h=mean(lam_1_use_v);
lam_2h=mean(lam_2_use_v);
u_v_use=u_v_gib(:,Burn_in:lag:end);
pi_v_use=pi_v_gib(:,:,Burn_in:lag:end);
u_v_h=mean(u_v_use,2);
pi_v_h=mean(pi_v_use,3);

CI=[quantile(u_v_use(1,:),[0.025,0.975]) quantile(u_v_use(2,:),[0.025,0.975]) quantile(lam_1_use_v,[0.025,0.975])  quantile(lam_2_use_v,[0.025,0.975]) ];
for i_ci=1:K_d+2
    if true_v(i_ci)<=CI(2*i_ci) && true_v(i_ci)>=CI(2*i_ci-1)
        cover_p(i_ci)=1;
    end
end 
% the percentage of correct membership
% membership of each driver, # of drivers in each group
memship = 0;
for j=1:m
    [pi_m,I]=max(pi_v_h(j,:));
    if I==tau_index(j)
        memship=memship+1;
    end
end