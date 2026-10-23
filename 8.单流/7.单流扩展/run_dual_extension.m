function audit=run_dual_extension(array,varargin)
[root,~]=dmimo.setup;c=dmimo.config;
p=inputParser;addParameter(p,'Mode','pilot');addParameter(p,'OutputDir','');addParameter(p,'Samples',[]);parse(p,varargin{:});opt=p.Results;
assert(ismember(array,{'ULA','UPA'})&&ismember(opt.Mode,{'smoke','pilot','full'}));
loadfun=str2func([lower(array) '.load_inputs']);localfun=str2func([lower(array) '.localize']);[net,sim,bank]=loadfun();cfg=posfb.config(net);
if strcmp(opt.Mode,'full'),n=400;axes=c.axes;
elseif strcmp(opt.Mode,'pilot')
 n=24;axes=[8 20;8 30;8 40];rng(c.pilot_seed,'twister');bank.positions=[20+10*rand(2,n);zeros(1,n)];
 if strcmp(array,'UPA'),bank.positions(3,:)=5*rand(1,n);end
 bank.zeta=(randn(n,6)+1j*randn(n,6))/sqrt(2);bank.noise_seed=1800000+reshape(1:n*10,n,10);
else,n=2;axes=[8 30];end
if ~isempty(opt.Samples),n=opt.Samples;end
if isempty(opt.OutputDir)
 if strcmp(opt.Mode,'full'),leaf=array;else,leaf=[opt.Mode '_verified_' array];end
 out=fullfile(root,'results','dual_v1_20261022',leaf);
else,out=opt.OutputDir;end
if ~exist(out,'dir'),mkdir(out);end
contract=struct('dual',c,'position',cfg,'mode',opt.Mode,'array',array,'samples',n,'axes',axes,'rho',sim.rho_data,...
 'sample_bank',sim.sample_bank_id,'pilot_replaces_positions',strcmp(opt.Mode,'pilot'));
file=fullfile(out,'experiment_contract.json');if exist(file,'file'),prior=jsondecode(fileread(file));assert(strcmp(prior.mode,opt.Mode)&&prior.samples==n);else,dmimo.json(file,contract);end
save(fullfile(out,'public_config.mat'),'net','sim','cfg','c','-v7');
diary(fullfile(out,'run.log'));clean=onCleanup(@()diary('off'));started=tic;names=columns;
for block=1:ceil(n/20)
 first=(block-1)*20+1;last=min(n,block*20);file=fullfile(out,sprintf('chunk_%03d.csv',block));
 if exist(file,'file'),fprintf('Reuse chunk %d\n',block);continue;end
 rows=cell(0,numel(names));
 for sid=first:last
  s=dmimo.sample(net,sim,bank,sid,c);
  for sc=1:numel(c.scenarios)
   scenario=c.scenarios{sc};blocks=dmimo.channel(net,sim,s,c,scenario);
   if sc==1,settings=axes;stage='MAIN';else,settings=[8 30];stage='MISMATCH';end
   for ci=1:2
    if ci==1,channel='LOS';else,channel='MIX';end
    for setting=1:size(settings,1)
     B=settings(setting,1);snr=settings(setting,2);active=net.active_order(1:B);
     [H,Y]=dmimo.observe(blocks,s.Z,channel,active,snr);[~,e]=dmimo.modes(H,2);
     freq=dmimo.frequency(net,Y,active,net.reference5);
     [Fd,codes]=dmimo.full_beams(net,Y,B);Et=dmimo.reduce(net,H,Fd);Ey=dmimo.reduce(net,Y,Fd);
     for ns=[1 2]
      C=dmimo.coefficients(Ey,ns,B);emit('FULLDIR_UE_Q',0,ns,Et,C,B*net.index_bits,NaN,false,[],codes,[],c.qA,c.qP,stage);
      cap=dmimo.capacity(e,ns,sim.rho_data);
      rows(end+1,:)={array,stage,scenario,channel,sid,B,snr,0,ns,'CAP_RANK_LE_N',0,0,NaN,NaN,NaN,cap,cap,20*cap,20*cap,cap,NaN,1,NaN,e(2)/e(1),NaN,false,'[]','[]','[]','[]','[]','[]'}; %#ok<AGROW>
      if sc==1&&B==8&&snr==30
       emit('FULLDIR_UE_CONT',0,ns,Et,C,B*net.index_bits,NaN,false,[],codes,[],0,0,'COEFFICIENT');
       for q=[2 4],emit('FULLDIR_UE_Q',0,ns,Et,C,B*net.index_bits,NaN,false,[],codes,[],q,q+1,'COEFFICIENT');end
      end
     end
     for R=[3 5]
      if R==3,refs=net.reference3;else,refs=net.reference5;end
      [~,cols]=ismember(refs,net.reference5);pu=localfun(net,refs,freq(:,cols));
      [message,diagnostic]=posfb.encode(pu,refs,cfg,cfg.main_bits);pred=posfb.decode(net,message,active,cfg);
      err=norm(pred.position-s.position);Et=dmimo.reduce(net,H,pred.F);Ey=dmimo.reduce(net,Y,pred.F);
      Ec=dmimo.cpu_model(net,pred.position,active,pred.F,struct('E',c.E));
      if sc>1
       calibration=struct('E',c.E);
       if ismember(scenario,{'tilt20','combined'}),calibration.ue_tilt_deg=20;end
       if ismember(scenario,{'phase15','combined'}),calibration.port_phase_rad=c.port_phase_rad;end
       Ecal=dmimo.cpu_model(net,pred.position,active,pred.F,calibration);
      end
      for ns=[1 2]
       Ccpu=dmimo.coefficients(Ec,ns,B);C=dmimo.coefficients(Ey,ns,B);
       emit('POS_CPU',R,ns,Et,Ccpu,cfg.main_bits,err,diagnostic.overflow,pred.position,pred.beam_indices,message.integer_codes,0,0,stage);
       emit('POS_UE_Q',R,ns,Et,C,cfg.main_bits,err,diagnostic.overflow,pred.position,pred.beam_indices,message.integer_codes,c.qA,c.qP,stage);
       if sc>1
        Ccal=dmimo.coefficients(Ecal,ns,B);
        emit('POS_CPU_CALIBRATED',R,ns,Et,Ccal,cfg.main_bits,err,diagnostic.overflow,pred.position,pred.beam_indices,message.integer_codes,0,0,stage);
       end
       if sc==1&&B==8&&snr==30
        emit('POS_UE_CONT',R,ns,Et,C,cfg.main_bits,err,diagnostic.overflow,pred.position,pred.beam_indices,message.integer_codes,0,0,'COEFFICIENT');
        for q=[2 4],emit('POS_UE_Q',R,ns,Et,C,cfg.main_bits,err,diagnostic.overflow,pred.position,pred.beam_indices,message.integer_codes,q,q+1,'COEFFICIENT');end
       end
      end
     end
    end
   end
  end
  if mod(sid,5)==0||sid==last,fprintf('DUAL %s %s %d/%d elapsed %.1f s\n',array,opt.Mode,sid,n,toc(started));end
 end
 T=cell2table(rows,'VariableNames',names);temp=[file '.pending.csv'];writetable(T,temp);movefile(temp,file,'f');
 dmimo.json(fullfile(out,'progress.json'),struct('status','running','samples_completed',last,'samples',n,'elapsed_seconds',toc(started)));
end
T=table;files=dir(fullfile(out,'chunk_*.csv'));for k=1:numel(files),T=[T;readtable(fullfile(out,files(k).name),'TextType','string')];end %#ok<AGROW>
writetable(T,fullfile(out,'trials.csv'));
realizable=~startsWith(T.method,'CAP_');finite=realizable&~isnan(T.feedback_bits);
expected=n*2*(size(axes,1)*12+3*16+18);
assert(height(T)==expected,'dmimo:RecordCount','Unexpected row count');
keys={'array','stage','scenario','channel','sample_id','B','observation_SNR_dB','R','streams','method','qA','qP'};
assert(height(unique(T(:,keys),'rows'))==height(T));
audit=struct('status','complete','array',array,'mode',opt.Mode,'samples',n,'records',height(T),'elapsed_seconds',toc(started),...
 'max_total_power_error',max(abs(T.precoder_norm_squared(realizable)-1)),...
 'max_per_AP_power_error',max(T.per_AP_power_error(realizable)),...
 'max_joint_capacity_excess',max(T.joint_SE(realizable)-T.rank_capacity(realizable)),...
 'coordinate_overflows',sum(T.quantizer_overflow),'minimum_finite_feedback',min(T.feedback_bits(finite)));
assert(audit.max_total_power_error<1e-10&&audit.max_per_AP_power_error<1e-10&&audit.max_joint_capacity_excess<1e-7);
dmimo.json(fullfile(out,'audit.json'),audit);dmimo.json(fullfile(out,'progress.json'),struct('status','complete','samples_completed',n,'samples',n,'elapsed_seconds',toc(started)));
fprintf('COMPLETE %s %s records=%d elapsed %.1f s\n',array,opt.Mode,height(T),toc(started));

 function emit(method,R,ns,Et,C,db,err,overflow,position,beam,indices,qA,qP,stageid)
  ai=[];pi=[];cb=0;
  if qA>0
   msg=dmimo.encode_coeff(C,qA,qP);[C,cb]=dmimo.decode_coeff(msg);ai=msg.amplitude_indices;pi=msg.phase_indices;
  elseif contains(method,'CONT'),cb=NaN;end
  met=dmimo.metrics(Et,C,sim.rho_data,e,ns);
  rows(end+1,:)={array,stageid,scenario,channel,sid,B,snr,R,ns,method,qA,qP,db,cb,db+cb,met.joint,met.mmse,...
   20*met.joint,20*met.mmse,met.capacity,met.power,met.norm,met.ap_error,met.ratio,err,overflow,...
   jsonencode(position),jsonencode(beam),jsonencode(indices),jsonencode(ai),jsonencode(pi),jsonencode(met.sinr)}; %#ok<AGROW>
 end
end
function names=columns
names={'array','stage','scenario','channel','sample_id','B','observation_SNR_dB','R','streams','method','qA','qP',...
 'direction_bits','coefficient_bits','feedback_bits','joint_SE','MMSE_SE','joint_Mbps','MMSE_Mbps','rank_capacity',...
 'received_power','precoder_norm_squared','per_AP_power_error','full_sigma2_ratio','position_error_m','quantizer_overflow',...
 'decoded_position','beam_indices','coordinate_indices','coefficient_amplitude_indices','coefficient_phase_indices','stream_SINR'};
end
