function report=verify_upa_ue
% Verify UE phase quantization, exact beam agreement and shared bit accounting.
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);[net,sim,bank]=upa.load_inputs;timer=tic;names={};
sample=upa.sample_channels(net,sim,bank,1);ids=net.active_order(1:8);
for channel={'LOS','LOS_NLOS_K10'}
 [Htrue,Hobs]=upa.observe(sample,channel{1},ids,30);[messages,~]=upa.encode_directions(net,Hobs,ids,net.reference5,sim.frequency_bits_per_dimension);
 [~,~,vobs]=upa.principal(Hobs);
 for R=[3 5]
  if R==3,refs=net.reference3;else,refs=net.reference5;end
  for method={'HAQSE','DFT'}
   msg=upa.subset_message(messages.(method{1}),refs);d=upa.ue_phase_from_message(net,msg,ids,Hobs,sim.UE_phase_bits);cpuPrediction=upa.predict_directions(net,msg,ids);
   check(isequal(d.beam_indices,cpuPrediction.beam_indices),sprintf('%s R%d beam agreement',method{1},R));
   check(norm(d.position-cpuPrediction.position)<1e-12,sprintf('%s R%d position agreement',method{1},R));
   check(strcmp(d.phase_source,'UE_estCSI'),sprintf('%s R%d phase source',method{1},R));
   check(d.phase_indices(1)==0&&numel(d.phase_indices)==8&&all(d.phase_indices>=0&d.phase_indices<16),sprintf('%s R%d phase integer range',method{1},R));
   a=zeros(8,1);for b=1:8,ix=(b-1)*net.Nap+(1:net.Nap);a(b)=d.F(:,b)'*vobs(ix);end
   theta=angle(exp(1j*(angle(a)-angle(a(1)))));k=mod(round(theta*16/(2*pi)),16);
   check(isequal(d.phase_indices(:),k(:)),sprintf('%s R%d CSI projection phase',method{1},R));
   check(abs(norm(d.w)^2-1)<1e-10,sprintf('%s R%d total power',method{1},R));
   power=reshape(abs(d.w).^2,net.Nap,8);check(max(abs(sum(power,1)-1/8))<1e-10,sprintf('%s R%d equal AP power',method{1},R));
   row=bank.budget(strcmp(bank.budget.array,net.array)&bank.budget.B==8&strcmp(bank.budget.method,sprintf('%s_R%d',method{1},R)),:);check(height(row)==1&&d.feedback_bits==row.UE_total_bits,sprintf('%s R%d manifest cost',method{1},R));
   check(d.AP_phase_bits==28&&d.feedback_bits==d.direction_bits+28,sprintf('%s R%d phase bit count',method{1},R));
   [truth.u,truth.s2]=upa.principal(Htrue);m=upa.evaluate(net,sim,Htrue,truth,d,sample.position);check(isfinite(m.SE)&&m.power_gain>=0&&m.power_gain<=1+1e-8,sprintf('%s R%d physical metric',method{1},R));
  end
 end
end
report=struct('passed',true,'checks',numel(names),'test_names',{names},'elapsed_seconds',toc(timer),'MATLAB_version',version,'note','UE phase branch interface and accounting checks; not full-sample performance');
fid=fopen(fullfile(root,'verification_ue.json'),'w');fwrite(fid,jsonencode(report),'char');fclose(fid);fprintf('PASS: %d UE checks, %.1f seconds\n',numel(names),toc(timer));
 function check(ok,name)
  assert(ok,'upa:UEVerification','%s',name);names{end+1}=name;
 end
end
