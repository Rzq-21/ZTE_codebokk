function report=verify_upa_cpu
% Verify physical conventions, frozen samples, feedback accounting and CPU boundary.
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);[net,sim,bank]=upa.load_inputs;
tests={};errors=[];timer=tic;
gold=load(fullfile(root,'shared','single_sample_alignment.mat'));gold=gold.out.UPA32x32;
sample=upa.sample_channels(net,sim,bank,1);ids=net.active_order(1:8);
for name={'LOS','LOS_NLOS_K10'}
 channel=name{1};[H,Hobs,Z,nv]=upa.observe(sample,channel,ids,30);
 if strcmp(channel,'LOS'),H0=gold.Htrue_LOS;O0=gold.Hobs_LOS_30dB;nv0=gold.noise_variance_LOS;else,H0=gold.Htrue_MIX;O0=gold.Hobs_MIX_30dB;nv0=gold.noise_variance_MIX;end
 check([channel ' truth'],max(abs(H(:)-H0(:))),1e-12);check([channel ' noisy CSI'],max(abs(Hobs(:)-O0(:))),1e-12);check([channel ' shared noise'],max(abs(Z(:)-gold.Z_standard(:))),1e-12);check([channel ' noise variance'],abs(nv-nv0),1e-16);
end
% Exact center-ray frequencies must recover the true position with either reference set.
for R=[3 5]
 if R==3,refs=net.reference3;else,refs=net.reference5;end
 freq=zeros(2,R);
 for j=1:R,b=refs(j);q=net.R(:,:,b)'*(net.pb(:,b)-sample.position)/norm(net.pb(:,b)-sample.position);freq(:,j)=-q(2:3)/2;end
 [p,~]=upa.localize(net,refs,freq);check(sprintf('Exact geometry R%d',R),norm(p-sample.position),1e-9);
end
[H,Hobs]=upa.observe(sample,'LOS',ids,30);[msg,diag]=upa.encode_directions(net,Hobs,ids,net.reference3,14);
reference=gold.HAQSE_LOS_reference_frequencies;
if isequal(size(reference),[3 2]),reference=reference.';end
check('HAQSE frozen reference',max(abs(diag.frequency_unquantized(:)-reference(:))),1e-10);
L=2^14;expected=max(0,min(L-1,round((diag.frequency_unquantized+.5)*(L-1))));check('HAQSE endpoint integer quantizer',max(abs(expected(:)-msg.HAQSE.integer_codes(:))),0);
endpoint=[-.5 .5-1e-12];k=round((endpoint+.5)*(L-1));check('Quantizer endpoints',max(abs(k-[0 L-1])),0);
[~,~,v]=upa.principal(Hobs(:,1:1024));[m,g]=upa.spectrum(net,v);f=upa.beams(net,m);check('FFT vs direct selected correlation',abs(f'*v-g(m(1)+128*m(2)+1)),1e-11);
F=upa.beams(net,[0 4;0 0]);check('Type-II orthogonal group',abs(F(:,1)'*F(:,2)),1e-12);
for Q=2:6
 [c,k]=upa.phase_weights(exp(1j*[-3 -.2 2.2].'),Q);check(sprintf('Phase quantizer Q%d',Q),max(abs(angle(c)-angle(exp(1j*2*pi*k/2^Q)))),1e-12);
 assert(k(1)==0&&all(k>=0&k<2^Q));
end
% Golden channels are not forwarded to the CPU. Reject even an accidental truth field.
bad=msg.DFT;bad.Htrue=H;rejected=false;
try,upa.cpu_from_message(net,bad,ids);catch ex,rejected=strcmp(ex.identifier,'upa:UnexpectedMessageField');end
assert(rejected,'CPU boundary failed');check('CPU rejects CSI in message',0,0);
for B=[5 8 10]
 ids=net.active_order(1:B);
 for channel={'LOS','LOS_NLOS_K10'}
  [H,Hobs]=upa.observe(sample,channel{1},ids,30);[truth.u,truth.s2]=upa.principal(H);designs=upa.design_baselines(net,Hobs,H,ids,4);messages=upa.encode_directions(net,Hobs,ids,net.reference5,14);
  for refs={net.reference3,net.reference5}
   for scheme={'HAQSE','DFT'}
    message=upa.subset_message(messages.(scheme{1}),refs{1});d=upa.cpu_from_message(net,message,ids);designs{end+1}=d;
    if strcmp(scheme{1},'DFT')
     for j=1:numel(refs{1}),slot=find(ids==refs{1}(j));check('Received DFT reference index preserved',max(abs(d.beam_indices(:,slot)-message.integer_codes(:,j))),0);end
    end
   end
  end
  for j=1:numel(designs)
   d=designs{j};metric=upa.evaluate(net,sim,H,truth,d,sample.position);assert(isfinite(metric.SE)&&isfinite(metric.power_gain));
   if strcmp(d.method_id,'SVD'),check('Ideal SVD power gain',abs(metric.power_gain-1),1e-11);else
    ref=bank.budget(strcmp(bank.budget.array,net.array)&bank.budget.B==B&strcmp(bank.budget.method,d.method_id),:);assert(height(ref)==1);check('Feedback manifest',abs(d.feedback_bits-ref.CPU_total_bits),0);
   end
  end
 end
end
% An unchanged truth channel gives exactly the same ideal precoder across observation SNRs.
ids=net.active_order(1:8);[H,O20]=upa.observe(sample,'LOS',ids,20);[~,O40]=upa.observe(sample,'LOS',ids,40);d20=upa.design_baselines(net,O20,H,ids,4);d40=upa.design_baselines(net,O40,H,ids,4);check('Ideal baseline independent of observation SNR',norm(d20{1}.w-d40{1}.w),1e-12);
report=struct('passed',true,'checks',numel(tests),'test_names',{tests},'absolute_errors',errors,'elapsed_seconds',toc(timer),'MATLAB_version',version,'note','Numerical and protocol tests; not 400-sample performance validation');
fid=fopen(fullfile(root,'verification.json'),'w');fwrite(fid,jsonencode(report),'char');fclose(fid);fprintf('PASS: %d checks, %.1f seconds\n',numel(tests),toc(timer));
 function check(name,error,tolerance)
  assert(isfinite(error)&&error<=tolerance,'upa:Verification','%s: error %.4g > %.4g',name,error,tolerance);tests{end+1}=name;errors(end+1)=error;
 end
end
