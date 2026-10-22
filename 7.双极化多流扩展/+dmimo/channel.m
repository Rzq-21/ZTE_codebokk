function blocks=channel(net,sim,s,c,scenario)
tilt=0;phi=zeros(1,10);
if ismember(scenario,{'tilt20','combined'}),tilt=20;end
if ismember(scenario,{'phase15','combined'}),phi=c.port_phase_rad;end
blocks.LOS=zeros(8,2*net.Nap,10);blocks.MIX=blocks.LOS;K=10^(sim.Rician_K_dB/10);
for b=1:10
 J=dmimo.jones(net,c,s.position,b,0,tilt);L=kron(J,s.L(:,:,b));A=zeros(size(L));
 for ell=1:6
  J=dmimo.jones(net,c,s.position,b,ell,tilt);A=A+kron(J,s.A(:,:,ell,b));
 end
 A=A*norm(L,'fro')/norm(A,'fro');M=sqrt(K/(K+1))*L+sqrt(1/(K+1))*A;
 ix=net.Nap+(1:net.Nap);L(:,ix)=L(:,ix)*exp(1j*phi(b));M(:,ix)=M(:,ix)*exp(1j*phi(b));
 blocks.LOS(:,:,b)=L;blocks.MIX(:,:,b)=M;
end
end
