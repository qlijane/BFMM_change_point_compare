function ci= hpdr(x)% compute 95% hpd CI
[f_lam_1,lam_1_xi] = ksdensity(x);% Kernal density estimate
% figure(6)
% hist(lam_1_use);
% hold on
% plot(lam_1_xi, f_lam_1)
%md=lam_1_xi(f_lam_1==max(f_lam_1));%mode
px=f_lam_1/sum(f_lam_1); %standardized frequency
pxs=sort(px,'descend'); 
ct=min(pxs(cumsum(pxs)<0.95));
if isempty(ct)==1 
    ct=0;
end
ci=[min(lam_1_xi(px>=ct)) max(lam_1_xi(px>=ct))];