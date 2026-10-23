function benchmark_position_extension(out)
% Warm-cache repeated CPU wall-clock timing; no hardware/PHY latency claims.
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);rows=cell(0,8);
for array={'ULA','UPA'}
 ar=array{1};if strcmp(ar,'ULA'),pkg='ula';else,pkg='upa';end
 loadfun=str2func([pkg '.load_inputs']);samplefun=str2func([pkg '.sample_channels']);localfun=str2func([pkg '.localize']);
 oldcpu=str2func([pkg '.cpu_from_message']);encodefun=str2func([pkg '.encode_directions']);
 [net,sim,bank]=loadfun();cfg=posfb.config(net);sample=samplefun(net,sim,bank,1);active=net.active_order(1:8);
 for channel={'LOS','LOS_NLOS_K10'}
  ch=channel{1};[~,Hobs]=upa.observe(sample,ch,active,30);[~,~,vobs]=upa.principal(Hobs);messages=encodefun(net,Hobs,active,net.reference5,14);
  for R=[3 5]
   if R==3,refs=net.reference3;else,refs=net.reference5;end
   freq=estimate(net,Hobs,active,refs);pu=localfun(net,refs,freq);msg=posfb.encode(pu,refs,cfg,cfg.main_bits);pred=posfb.decode(net,msg,active,cfg);oldmsg=upa.subset_message(messages.HAQSE,refs);
   functions={@()estimate(net,Hobs,active,refs),@()localfun(net,refs,freq),@()posfb.encode(pu,refs,cfg,cfg.main_bits),...
    @()posfb.decode(net,msg,active,cfg),@()posfb.cpu_phase(net,pred,active),@()ue_all(net,Hobs,pred),@()oldcpu(net,oldmsg,active)};
   names={'UE_HAQSE_estimation','UE_ray_localization','UE_coordinate_encoding','CPU_coordinate_decode_and_beams','CPU_model_co_phase','UE_observed_target_and_co_phase','legacy_CPU_HAQSE_decode_and_co_phase'};
   for j=1:numel(functions)
    fn=functions{j};fn();times=zeros(31,1);for rep=1:31,tic;fn();times(rep)=toc*1000;end
    rows(end+1,:)={ar,ch,R,names{j},median(times),prctile(times,10),prctile(times,90),31}; %#ok<AGROW>
   end
  end
 end
end
t=cell2table(rows,'VariableNames',{'array','channel_model','R','operation','median_ms','p10_ms','p90_ms','repetitions'});writetable(t,fullfile(out,'runtime_benchmark.csv'));fprintf('Timing complete (%d operation groups)\n',height(t));
end
function f=estimate(net,Hobs,active,refs)
d=1+strcmp(net.array,'UPA32x32');f=zeros(d,numel(refs));
for j=1:numel(refs)
 slot=find(active==refs(j));Hb=Hobs(:,(slot-1)*net.Nap+(1:net.Nap));
 if d==1,f(j)=haqse_ula_reference(Hb.');else,x=MMV_CF_mQSE2_SumOfRatios(reshape(Hb.',[32,32,4]),1,1,1,3);f(:,j)=mod(x(1,:).'+.5,1)-.5;end
end
end
function d=ue_all(net,Hobs,pred)
[~,~,vobs]=upa.principal(Hobs);d=posfb.ue_phase(net,pred,vobs,4);
end
