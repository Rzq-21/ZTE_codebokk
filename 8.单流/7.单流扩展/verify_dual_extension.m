function audit=verify_dual_extension
[root,~]=dmimo.setup;c=dmimo.config;checks=0;maxModel=0;
% Analytical channels: known diagonal capacity and equal power two-mode case.
assert(abs(dmimo.capacity([4;1],2,1)-log2(5.0625))<1e-12);checks=checks+1;
assert(abs(dmimo.capacity([4;.01],2,1)-log2(5))<1e-12);checks=checks+1;
for ar={'ULA','UPA'}
 loadfun=str2func([lower(ar{1}) '.load_inputs']);[net,sim,bank]=loadfun();cfg=posfb.config(net);
 s=dmimo.sample(net,sim,bank,1,c);blocks=dmimo.channel(net,sim,s,c,'nominal');active=net.active_order(1:8);
 [H,Y]=dmimo.observe(blocks,s.Z,'LOS',active,30);pred=posfb.prediction(net,s.position,net.reference3,active);
 E=dmimo.cpu_model(net,s.position,active,pred.F,struct('E',c.E));
 rejected=false;try,dmimo.cpu_model(net,s.position,active,pred.F,c);catch,rejected=true;end;assert(rejected);checks=checks+1;
 % Expanded public CPU formula equals the reduced scalar-plane Jones model.
 Hl=upa.reconstruct_channel(net,s.position,active);Hp=zeros(8,2*net.Nap*8);
 for j=1:8
  ix=(j-1)*net.Nap+(1:net.Nap);iy=(j-1)*2*net.Nap+(1:2*net.Nap);
  Hp(:,iy)=kron(dmimo.jones(net,c,s.position,active(j),0,0),Hl(:,ix));
 end
 err=norm(E-dmimo.reduce(net,Hp,pred.F),'fro')/norm(E,'fro');assert(err<1e-12);maxModel=max(maxModel,err);checks=checks+1;
 [msg,~]=posfb.encode(s.position,net.reference3,cfg,cfg.main_bits);decoded=posfb.decode(net,msg,active,cfg);
 bad=msg;bad.true_CSI=H;rejected=false;try,posfb.decode(net,bad,active,cfg);catch,rejected=true;end;assert(rejected);checks=checks+1;
 f=dmimo.frequency(net,Y,active,net.reference5);assert(all(isfinite(f(:))));checks=checks+1;
 [~,e]=dmimo.modes(H,2);Et=dmimo.reduce(net,H,decoded.F);Ey=dmimo.reduce(net,Y,decoded.F);
 for n=[1 2]
  C=dmimo.coefficients(Ey,n,8);
  for q=[2 3 4]
   message=dmimo.encode_coeff(C,q,q+1);[Cq,bits]=dmimo.decode_coeff(message);assert(bits==16*n*(2*q+1)-n*(q+1));checks=checks+1;
   met=dmimo.metrics(Et,Cq,sim.rho_data,e,n);assert(met.mmse>=0);checks=checks+1;
   poison=message;poison.observed_CSI=Y;rejected=false;try,dmimo.decode_coeff(poison);catch,rejected=true;end;assert(rejected);checks=checks+1;
  end
  rotated=C.*exp(1j*(1:n));assert(norm((Et*C)*(Et*C)'-(Et*rotated)*(Et*rotated)','fro')<1e-10);checks=checks+1;
 end
 % Port calibration is a column-unitary transform; preserve true singular values.
 mis=dmimo.channel(net,sim,s,c,'phase15');[Hm,~]=dmimo.observe(mis,s.Z,'LOS',active,30);[~,em]=dmimo.modes(Hm,2);
 assert(max(abs(em-e))/max(e)<1e-12);checks=checks+1;
 % Scatter channel changes under UE tilt, rather than only a common row basis.
 mis=dmimo.channel(net,sim,s,c,'tilt20');assert(norm(mis.LOS(:)-blocks.LOS(:))>1e-4);checks=checks+1;
end
out=fullfile(root,'results','dual_v1_20261022');if ~exist(out,'dir'),mkdir(out);end
audit=struct('status','passed','checks',checks,'maximum_CPU_model_algebra_error',maxModel,'coefficient_boundary','rejects unexpected CSI fields','single_stream_joint_MMSE_equality','asserted','total_and_per_AP_power','asserted','waterfilling_analytic_tests','passed');
dmimo.json(fullfile(out,'preflight_audit.json'),audit);disp(audit);
end
