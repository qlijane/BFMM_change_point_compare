function u_s=ra_sampl_u(a,b,lam_1,lam_2,C_u,Nj_u,z_u,k_n,pi_u)
% rejection and acceptance sampling from u
ginv=b-a;% for rejection sampling of mu
              x1=a:b;%u is between a and b
              i=1;
                logf=zeros(1,length(x1));
            for u=x1
                logf(i)=log_L_multi(lam_1,lam_2,u,C_u,Nj_u,z_u,k_n,pi_u);
                i=i+1;
            end
            adju_c= log(1/ginv)-(max(logf));
           % logf_a=logf+adju_c;% adjust by constant for computational purpose
            %c_u=exp(max(logf_a))*ginv;% how to choose c? 
            c_u=1;
                  accept = false;
                while accept == false
                  u = rand();
                  v = a + (b-a)*rand();             
                    f1=log_L_multi(lam_1,lam_2,v,C_u,Nj_u,z_u,k_n,pi_u)+adju_c;
                    if c_u*u <= exp(f1)*ginv
                        u_s= v;
                        accept = true;
                    end
                end