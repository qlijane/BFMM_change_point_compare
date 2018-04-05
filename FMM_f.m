% finite mixture for simulation
function [DIC,indi]=FMM_f(z,Nj,C,K,m)
Bt=10000;%7000;%# of samples for Gibbs
Burn_in=5000;%5000;% burn in time
lag=2;%every 5th result is used
%-----------
u_v_i=[100 180 260 300 ];
a=10*ones(1,4); b=430*ones(1,4);
b1=1;b2=1;%prior of lam_b,lam_a
 a1=0.25;a2=0.1;
alpha0=1.8;memship=zeros(m,1);
 %________________________Gibbs start__________________________
lam_1_gib=zeros(m,Bt);lam_2_gib=zeros(m,Bt);
lam_1_gib(:,1)=0.05; lam_2_gib(:,1)=0.01; 
% # of uniquce values for model
pi_v_gib=zeros(m,K,Bt);% start from uniform, mixture proportion
u_v_gib=zeros(K,Bt);% unique values
u_v_gib(:,1)=u_v_i(1:K);
%---------------------------
% initialize randomly
zj=unidrnd(K,m,1);tau_gib=u_v_gib(1,1)*ones(m,Bt);
% zj is the index of tau_j in uk;
i_gib=zeros(m,K,Bt);
for k_t=2:K
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
       a_gib=alpha0/K+i_gib(:,:,t_g-1);
        for j=1:m
        pi_v_gib(j,:,t_g)=drchrnd(a_gib(j,:),1);%update pi
        end
% update u

           for k_u=1:K
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
               log_L=zeros(1,K);
               for k=1:K
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
 lam_1h=mean(lam_1_use,2);
lam_2h=mean(lam_2_use,2);
len_use=length(Burn_in:lag:Bt);
u_v_use=u_v_gib(:,Burn_in:lag:end);
tau_use=tau_gib(:,Burn_in:lag:end);
pi_v_use=pi_v_gib(:,:,Burn_in:lag:end);
u_v_h=mean(u_v_use,2);
pi_v_h=mean(pi_v_use,3);
Dev_use=zeros(1,len_use);
indi=0;
for j=1:m
    [pi_m,memship(j)]=max(pi_v_h(j,:));
end
tau_h=u_v_h(memship);
%***********************************************************
    for i=1:len_use
   Dev_use(i)=deviance(lam_1_use(:,i),lam_2_use(:,i),tau_use(:,i),C,Nj,z,pi_v_use(:,:,i)) ;
    end 
    DIC=2*mean(Dev_use)-deviance(lam_1h,lam_2h,tau_h,C,Nj,z,pi_v_h);
% u_v_h_sort=sort(u_v_h);
% u_v_diff=abs(diff(u_v_h_sort));
% lam_diff=abs(lam_1h-lam_2h);
%     if sum(u_v_diff<55)>0 || sum(lam_diff<0.035)>0
%     DIC=100000;
%     end
    %--------------
