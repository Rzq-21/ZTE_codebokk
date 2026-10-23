function report=verify_ula
% Frozen-input, geometry, beam, phase and feedback-accounting checks for ULA64.
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);[net,sim,bank]=ula.load_inputs;timer=tic;names={};errors=[];
gold=load(fullfile(root,'shared','single_sample_alignment.mat'));gold=gold.out.ULA64;sample=ula.sample_channels(net,sim,bank,1);ids=net.active_order(1:8);
for name={'LOS','LOS_NLOS_K10'}
 channel=name{1};[H,Hobs,Z,nv]=ula.observe(sample,channel,ids,30);
 if strcmp(channel,'LOS'),H0=gold.Htrue_LOS;O0=gold.Hobs_LOS_30dB;nv0=gold.noise_variance_LOS;else,H0=gold.Htrue_MIX;O0=gold.Hobs_MIX_30dB;nv0=gold.noise_variance_MIX;end
 check([channel ' truth'],max(abs(H(:)-H0(:))),1e-12);check([channel ' noisy CSI'],max(abs(Hobs(:)-O0(:))),1e-12);check([channel ' shared noise'],max(abs(Z(:)-gold.Z_standard(:))),1e-12);check([channel ' noise variance'],abs(nv-nv0),1e-16);
end
for R=[3 5]
 if R==3,refs=net.reference3;else,refs=net.reference5;end;freq=zeros(1,R);
 for j=1:R,b=refs(j);q=net.R(:,:,b)'*(net.pb(:,b)-sample.position)/norm(net.pb(:,b)-sample.position);freq(j)=-q(2)/2;end
 [p,~]=ula.localize(net,refs,freq);check(sprintf('Exact geometry R%d',R),norm(p-sample.position),1e-9);
end
[H,Hobs]=ula.observe(sample,'LOS',ids,30);[msg,diag]=ula.encode_directions(net,Hobs,ids,net.reference3,sim.frequency_bits_per_dimension);check('HAQSE frozen reference',max(abs(diag.frequency_unquantized(:)-gold.HAQSE_LOS_reference_frequencies(:))),1e-10);
L=2^14;k=mod(round(L*diag.frequency_unquantized),L);check('Circular frequency quantizer',max(abs(msg.HAQSE.integer_codes(:)-k(:))),0);
[~,~,v]=upa.principal(Hobs(:,1:64));[m,g]=ula.spectrum(net,v);f=ula.beams(net,m);check('FFT vs direct selected correlation',abs(f'*v-g(m+1)),1e-11);
for Q=2:6,[c,k]=upa.phase_weights(exp(1j*[-3 -.2 2.2].'),Q);check(sprintf('Phase quantizer Q%d',Q),max(abs(angle(c)-angle(exp(1j*2*pi*k/2^Q)))),1e-12);assert(k(1)==0&&all(k>=0&k<2^Q));end
bad=msg.DFT;bad.Htrue=H;rejected=false;try,ula.cpu_from_message(net,bad,ids);catch ex,rejected=strcmp(ex.identifier,'ula:UnexpectedMessageField');end;assert(rejected);check('CPU rejects CSI in message',0,0);
for B=[5 8 10]
 ids=net.active_order(1:B);
 for channel={'LOS','LOS_NLOS_K10'}
  [H,Hobs]=ula.observe(sample,channel{1},ids,30);[truth.u,truth.s2]=upa.principal(H);designs=ula.design_baselines(net,Hobs,H,ids,4);messages=ula.encode_directions(net,Hobs,ids,net.reference5,14);
  for refs={net.reference3,net.reference5}
   for scheme={'HAQSE','DFT'}
    message=upa.subset_message(messages.(scheme{1}),refs{1});cpu=ula.cpu_from_message(net,message,ids);ue=ula.ue_phase_from_message(net,message,ids,Hobs,4);check(sprintf('%s R%d beam parity',scheme{1},numel(refs{1})),max(abs(cpu.beam_indices(:)-ue.beam_indices(:))),0);check(sprintf('%s R%d position parity',scheme{1},numel(refs{1})),norm(cpu.position-ue.position),0);
    if strcmp(scheme{1},'DFT'),for j=1:numel(refs{1}),slot=find(ids==refs{1}(j));check('Received DFT reference index preserved',abs(cpu.beam_indices(slot)-message.integer_codes(j)),0);end,end
    check(sprintf('%s R%d UE phase rule',scheme{1},numel(refs{1})),abs(ue.phase_indices(1))+max([0;ue.phase_indices(:)-15]),0);check(sprintf('%s R%d power norm',scheme{1},numel(refs{1})),abs(norm(ue.w)^2-1),1e-10);
    power=reshape(abs(ue.w).^2,net.Nap,B);check('Equal AP power',max(abs(sum(power,1)-1/B)),1e-10);
    row=bank.budget(strcmp(bank.budget.array,net.array)&bank.budget.B==B&strcmp(bank.budget.method,sprintf('%s_R%d',scheme{1},numel(refs{1}))),:);check(sprintf('%s R%d UE feedback manifest',scheme{1},numel(refs{1})),abs(ue.feedback_bits-row.UE_total_bits),0);
    check(sprintf('%s R%d CPU phase cost',scheme{1},numel(refs{1})),cpu.AP_phase_bits,0);check(sprintf('%s R%d UE phase cost',scheme{1},numel(refs{1})),ue.AP_phase_bits,(B-1)*4);
    metric=upa.evaluate(net,sim,H,truth,ue,sample.position);assert(isfinite(metric.SE)&&metric.power_gain>=0&&metric.power_gain<=1+1e-8);
   end
  end
  for j=1:numel(designs),d=designs{j};metric=upa.evaluate(net,sim,H,truth,d,sample.position);assert(isfinite(metric.SE)&&isfinite(metric.power_gain));if ~strcmp(d.method_id,'SVD'),row=bank.budget(strcmp(bank.budget.array,net.array)&bank.budget.B==B&strcmp(bank.budget.method,d.method_id),:);check([d.method_id ' baseline feedback manifest'],d.feedback_bits-row.UE_total_bits,0);end,end
 end
end
report=struct('passed',true,'checks',numel(names),'test_names',{names},'absolute_errors',errors,'elapsed_seconds',toc(timer),'MATLAB_version',version,'note','ULA64 frozen channels, direction estimator, decoder, phase feedback and accounting checks; not 400-position performance validation');fid=fopen(fullfile(root,'verification_ula.json'),'w');fwrite(fid,jsonencode(report),'char');fclose(fid);fprintf('PASS: %d ULA checks, %.1f seconds\n',numel(names),toc(timer));
 function check(name,error,tolerance),assert(isfinite(error)&&error<=tolerance,'ula:Verification','%s: error %.4g > %.4g',name,error,tolerance);names{end+1}=name;errors(end+1)=error;end
end
