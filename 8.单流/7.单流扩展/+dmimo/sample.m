function s=sample(net,sim,bank,sid,c)
% Cache scalar paths; generate every polarization scenario without new paths.
N=net.Nap;p=bank.positions(:,sid);s.position=p;s.L=zeros(4,N,10);s.A=zeros(4,N,6,10);s.Z=zeros(8,2*N,10);
scatter=sim.scatterer_positions_m.';
for b=1:10
 ap=net.ap_elements(:,:,b);
 for u=1:4
  user=p+net.ue_offsets(:,u);r=vecnorm(ap-user);s.L(u,:,b)=exp(-1j*net.k*r)./r;
  for ell=1:6
   len=vecnorm(ap-scatter(:,ell))+norm(scatter(:,ell)-user);
   s.A(u,:,ell,b)=bank.zeta(sid,ell)*exp(-1j*net.k*len)./len;
  end
 end
 rng(bank.noise_seed(sid,b)+c.noise_seed_offset,'twister');s.Z(:,:,b)=(randn(8,2*N)+1j*randn(8,2*N))/sqrt(2);
end
end
