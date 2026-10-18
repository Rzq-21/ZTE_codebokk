function audit=verify_position_extension(outdir)
% Regression against archived v2 powers, plus actual integer-message boundary.
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);if ~exist(outdir,'dir'),mkdir(outdir);end
checks=0;maximum=0;beamdiff=0;positiondiff=0;allocrows=cell(0,4);
for array={'ULA','UPA'}
 ar=array{1};if strcmp(ar,'ULA'),pkg='ula';else,pkg='upa';end
 loadfun=str2func([pkg '.load_inputs']);samplefun=str2func([pkg '.sample_channels']);encodefun=str2func([pkg '.encode_directions']);localfun=str2func([pkg '.localize']);
 cpufun=str2func([pkg '.cpu_from_message']);uefun=str2func([pkg '.ue_phase_from_message']);basefun=str2func([pkg '.design_baselines']);
 [net,sim,bank]=loadfun();cfg=posfb.config(net);sample=samplefun(net,sim,bank,1);active=net.active_order(1:8);
 old=cell(1,2);for b=1:2,branches={'CPU','UE'};old{b}=readtable(fullfile(root,'results',sprintf('full400_%s_%s_20261017',ar,branches{b}),'results_trials.csv'),'VariableNamingRule','preserve');end
 for channel={'LOS','LOS_NLOS_K10'}
  ch=channel{1};[Htrue,Hobs]=upa.observe(sample,ch,active,30);[truth.u,truth.s2]=upa.principal(Htrue);[~,~,vobs]=upa.principal(Hobs);
  [messages,d]=encodefun(net,Hobs,active,net.reference5,14);bases=basefun(net,Hobs,Htrue,active,4);
  for R=[3 5]
   if R==3,refs=net.reference3;else,refs=net.reference5;end
   [~,ix]=ismember(refs,net.reference5);[pu,loc]=localfun(net,refs,d.frequency_unquantized(:,ix));
   for C=cfg.budgets
    [msg,dd]=posfb.encode(pu,refs,cfg,C);pred=posfb.decode(net,msg,active,cfg);[cpu,ue]=posfb.phases(net,pred,active,vobs,4);
    assert(isequal(cpu.beam_indices,ue.beam_indices)&&isequal(cpu.position,ue.position));
    assert(cpu.feedback_bits==C&&ue.feedback_bits==C+28&&sum(msg.coordinate_bits)==C);
    qerror=abs(pred.position(1:cfg.dim)-pu(1:cfg.dim));steps=(cfg.hi-cfg.lo)./2.^msg.coordinate_bits;
    if ~dd.overflow,assert(all(qerror<=steps/2+1e-12));end
    upa.evaluate(net,sim,Htrue,truth,cpu,sample.position);upa.evaluate(net,sim,Htrue,truth,ue,sample.position);checks=checks+5;
   end
   fake=msg;fake.true_position=sample.position;rejected=false;try,posfb.decode(net,fake,active,cfg);catch,rejected=true;end;assert(rejected);checks=checks+1;
   [sat,sd]=posfb.encode([100;100;100],refs,cfg,cfg.main_bits);assert(sd.overflow);posfb.decode(net,sat,active,cfg);checks=checks+1;
   % Continuous reference: same position must reproduce decoder-generated beams.
   pcont=posfb.prediction(net,pred.position,refs,active);assert(isequal(pcont.beam_indices,pred.beam_indices));checks=checks+1;
   % CPU and UE legacy methods must reproduce all archived sample-1 powers.
   for method={'HAQSE','DFT'}
    msg=upa.subset_message(messages.(method{1}),refs);cpu=cpufun(net,msg,active);ue=uefun(net,msg,active,Hobs,4);designs={cpu,ue};
    assert(isequal(cpu.beam_indices,ue.beam_indices));beamdiff=max(beamdiff,max(abs(cpu.beam_indices(:)-ue.beam_indices(:))));positiondiff=max(positiondiff,norm(cpu.position-ue.position));
    for b=1:2
     metric=upa.evaluate(net,sim,Htrue,truth,designs{b},sample.position);t=old{b};mask=t.sample_id==1&t.B==8&t.observation_SNR_dB==30&strcmp(t.channel_model,ch)&strcmp(t.method_id,designs{b}.method_id);
     assert(sum(mask)==1);maximum=max(maximum,abs(metric.received_power-t.received_power(mask)));checks=checks+1;
    end
   end
  end
  for j=1:numel(bases)
   metric=upa.evaluate(net,sim,Htrue,truth,bases{j},sample.position);
   for b=1:2,t=old{b};mask=t.sample_id==1&t.B==8&t.observation_SNR_dB==30&strcmp(t.channel_model,ch)&strcmp(t.method_id,bases{j}.method_id);assert(sum(mask)==1);maximum=max(maximum,abs(metric.received_power-t.received_power(mask)));checks=checks+1;end
  end
 end
 allocrows(end+1,:)={ar,cfg.main_bits,jsonencode(cfg.main_allocation),jsonencode(cfg.sensitivity_weights)}; %#ok<AGROW>
 fprintf('%s verified; main bits %s\n',ar,mat2str(cfg.main_allocation.'));
end
assert(maximum<1e-10&&beamdiff==0&&positiondiff<1e-9);
audit=struct('checks',checks,'max_legacy_power_difference',maximum,'CPU_UE_beam_index_difference',beamdiff,'CPU_UE_position_difference_m',positiondiff,'status','passed');
fid=fopen(fullfile(outdir,'preflight_audit.json'),'w');fwrite(fid,jsonencode(audit),'char');fclose(fid);
writetable(cell2table(allocrows,'VariableNames',{'array','coordinate_budget','coordinate_axis_bits','sensitivity_weights'}),fullfile(outdir,'main_allocations.csv'));
disp(audit);
end
