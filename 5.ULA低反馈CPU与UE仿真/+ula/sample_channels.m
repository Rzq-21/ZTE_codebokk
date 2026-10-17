function sample=sample_channels(net,sim,bank,sample_id)
% Exact element-distance spherical-wave channels and frozen standard noise.
p=bank.positions(:,sample_id);sample.position=p;sample.LOS=zeros(4,64,10);sample.MIX=sample.LOS;sample.Z=sample.LOS;K=10^(sim.Rician_K_dB/10);scatter=sim.scatterer_positions_m.';
for b=1:10
 ap=net.ap_elements(:,:,b);Hl=zeros(4,64);A=Hl;
 for u=1:4
  user=p+net.ue_offsets(:,u);r=vecnorm(ap-user);Hl(u,:)=exp(-1j*net.k*r)./r;
  for ell=1:6,L=vecnorm(ap-scatter(:,ell))+norm(scatter(:,ell)-user);A(u,:)=A(u,:)+bank.zeta(sample_id,ell)*exp(-1j*net.k*L)./L;end
 end
 assert(norm(A,'fro')>0);A=A*norm(Hl,'fro')/norm(A,'fro');sample.LOS(:,:,b)=Hl;sample.MIX(:,:,b)=sqrt(K/(K+1))*Hl+sqrt(1/(K+1))*A;
 rng(bank.noise_seed(sample_id,b),'twister');sample.Z(:,:,b)=(randn(4,64)+1j*randn(4,64))/sqrt(2);
end
end
