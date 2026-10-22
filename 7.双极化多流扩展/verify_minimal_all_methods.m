function verify_minimal_all_methods
% In the rank-one, one-polarization limit, reproduce old Type-I/II messages.
[root,~]=dmimo.setup;checks=0;maxerr=0;
for ar={'ULA','UPA'}
 loadfun=str2func([lower(ar{1}) '.load_inputs']);oldfun=str2func([lower(ar{1}) '.design_baselines']);[net,~,~]=loadfun();B=8;
 rng(2026102201,'twister');v=(randn(B*net.Nap,1)+1j*randn(B*net.Nap,1));v=v/norm(v);
 H=zeros(4,B*net.Nap);H(1,:)=v';Y=zeros(8,2*B*net.Nap);
 for b=1:B,ix=(b-1)*net.Nap+(1:net.Nap);iy=(b-1)*2*net.Nap+(1:net.Nap);Y(1:4,iy)=H(:,ix);end
 old=oldfun(net,H,H,net.active_order(1:B),4);messages=min9.baseline_messages(net,Y,B);
 for kind={'TypeI','TypeII'}
  if strcmp(kind{1},'TypeI'),j=3;else,j=2;end
  msg=messages.(kind{1});[F,db,ib]=min9.decode_baseline(net,msg,B);
  err=norm(F-old{j}.F,'fro');assert(err<1e-10);maxerr=max(maxerr,err);checks=checks+1;
  assert(db==old{j}.direction_bits&&ib==old{j}.internal_bits);checks=checks+1;
 end
end
out=fullfile(root,'results','minimal_all_methods_20261022');dmimo.json(fullfile(out,'legacy_reduction_audit.json'),struct('status','passed','checks',checks,'maximum_old_baseline_beam_error',maxerr,'meaning','rank-one and one-polarization limit reproduces original TypeI/II spatial messages and bit counts'));
disp('Original Type-I/II rank-one reduction checks passed');
end
