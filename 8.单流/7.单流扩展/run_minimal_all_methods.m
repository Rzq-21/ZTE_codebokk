function audit=run_minimal_all_methods(mode)
% All original method families; fixed B=8, observation SNR=30dB, ranks 1/2.
if nargin<1,mode='confirm';end
[root,base]=dmimo.setup;c=dmimo.config;
assert(ismember(mode,{'smoke','confirm'}));if strcmp(mode,'smoke'),ids=[399 400];else,ids=1:32;end
out=fullfile(root,'results','minimal_all_methods_20261022',mode);if ~exist(out,'dir'),mkdir(out);end
contract=struct('version','ALL9_MINIMAL_V1','mode',mode,'sample_ids',ids,'B',8,'observation_SNR_dB',30,'ranks',[1 2],...
 'frequency_bits_per_dimension',14,'qA',3,'qP',4,'scenario','nominal','sample_selection','fixed before results; original bank subset; no tuning',...
 'TypeII','two orthogonal spatial beams combined with 3bit amplitude and 2bit phase per AP, common to two polarizations; then 2xNs coefficients',...
 'baseline_selection','leading two observed modes shared across ranks; original one-mode score generalized by summed mode and polarization energy',...
 'references','SVD equal stream power and true-channel waterfilling capacity relax per-AP constraint',...
 'dependencies',base,'dual_model',c);
dmimo.json(fullfile(out,'contract.json'),contract);rows=cell(0,numel(columns));started=tic;checks=0;maxRegression=0;
for ar={'ULA','UPA'}
 array=ar{1};loadfun=str2func([lower(array) '.load_inputs']);predictfun=str2func([lower(array) '.predict_directions']);localfun=str2func([lower(array) '.localize']);
 [net,sim,bank]=loadfun();cfg=posfb.config(net);active=net.active_order(1:8);B=8;
 for sid=ids
  s=dmimo.sample(net,sim,bank,sid,c);blocks=dmimo.channel(net,sim,s,c,'nominal');
  for channelcell={'LOS','MIX'}
   channel=channelcell{1};[H,Y]=dmimo.observe(blocks,s.Z,channel,active,30);[Vtrue,e]=dmimo.modes(H,2);
   [messages,frequency]=min9.directions(net,Y,active,net.reference5);baselines=min9.baseline_messages(net,Y,B);
   for R=[3 5]
    if R==3,refs=net.reference3;else,refs=net.reference5;end
    [~,cols]=ismember(refs,net.reference5);
    for family={'HAQSE','DFT','POS'}
     kind=family{1};
     if strcmp(kind,'POS')
      pu=localfun(net,refs,frequency(:,cols));[msg,enc]=posfb.encode(pu,refs,cfg,cfg.main_bits);pred=posfb.decode(net,msg,active,cfg);overflow=enc.overflow;
     else
      msg=messages.(kind);msg.reference_ap_ids=refs;msg.integer_codes=msg.integer_codes(:,cols);pred=predictfun(net,msg,active);overflow=false;
     end
     F=pred.F;Et=dmimo.reduce(net,H,F);Ey=dmimo.reduce(net,Y,F);Ec=dmimo.cpu_model(net,pred.position,active,F,struct('E',c.E));
     for ns=[1 2]
      for branchcell={'CPU','UE'}
       branch=branchcell{1};if strcmp(branch,'CPU'),source=Ec;else,source=Ey;end
       C=dmimo.coefficients(source,ns,B);emit([kind '_R' num2str(R)],branch,R,ns,F,Et,C,pred.direction_bits,0,msg,pred.position,overflow);
      end
     end
    end
   end
   for kindcell={'TypeI','TypeII'}
    kind=kindcell{1};msg=baselines.(kind);[F,db,ib]=min9.decode_baseline(net,msg,B);Et=dmimo.reduce(net,H,F);Ey=dmimo.reduce(net,Y,F);
    if strcmp(kind,'TypeI'),[F0,~]=dmimo.full_beams(net,Y,B);assert(norm(F-F0,'fro')<1e-12);checks=checks+1;end
    if strcmp(mode,'smoke')
     poison=msg;poison.CSI=Y;reject=false;try,min9.decode_baseline(net,poison,B);catch,reject=true;end;assert(reject);checks=checks+1;
    end
    for ns=[1 2],C=dmimo.coefficients(Ey,ns,B);emit(kind,'UE_shared',0,ns,F,Et,C,db,ib,msg,nan(3,1),false);end
   end
   for ns=[1 2]
    G=H*Vtrue(:,1:ns)/sqrt(ns);[joint,mmse,sinr]=rates(G,sim.rho_data);cap=dmimo.capacity(e,ns,sim.rho_data);assert(joint<=cap+1e-8);
    addrow('SVD_equal','ideal',0,ns,NaN,NaN,NaN,NaN,joint,mmse,cap,1,NaN,nan(3,1),false,G,sinr,struct,struct);
    addrow('CAP_rank_le_N','ideal',0,ns,NaN,NaN,NaN,NaN,cap,cap,cap,1,NaN,nan(3,1),false,zeros(size(G)),nan(size(sinr)),struct,struct);
   end
  end
  if mod(sid,8)==0||sid==ids(end),fprintf('MIN9 %s %s sid=%d elapsed %.1f s\n',mode,array,sid,toc(started));end
 end
end
T=cell2table(rows,'VariableNames',columns);assert(height(T)==numel(ids)*4*32);
assert(height(unique(T(:,{'array','channel','sample_id','method','branch','streams'}),'rows'))==height(T));
writetable(T,fullfile(out,'trials.csv'));[g,S]=findgroups(T(:,{'array','channel','method','branch','R','streams'}));
S.samples=splitapply(@numel,T.MMSE_Mbps,g);
for field={'MMSE_Mbps','joint_Mbps','feedback_bits','position_error_m','received_mode_ratio','weakest_stream_SINR'}
 name=field{1};S.(['mean_' name])=splitapply(@mean,T.(name),g);
end
S.min_received_mode_ratio=splitapply(@min,T.received_mode_ratio,g);writetable(S,fullfile(out,'summary.csv'));
realizable=~strcmp(T.branch,'ideal');dual=realizable&T.streams==2;
assert(all(T.received_mode_ratio(dual)>1e-10),'min9:RankCollapse','Effective received signal lost the second mode');
assert(max(abs(T.precoder_norm_squared(realizable)-1))<1e-10&&max(T.per_AP_power_error(realizable))<1e-10);
audit=struct('status','passed','mode',mode,'samples_per_case',numel(ids),'records',height(T),'method_families',9,...
 'extra_capacity_reference',true,'realizable_dual_records',sum(dual),'minimum_dual_received_mode_ratio',min(T.received_mode_ratio(dual)),...
 'total_power_max_error',max(abs(T.precoder_norm_squared(realizable)-1)),'per_AP_power_max_error',max(T.per_AP_power_error(realizable)),...
 'overflow_records',sum(T.quantizer_overflow),'additional_checks',checks,'max_shared_POS_regression_Mbps',maxRegression,'elapsed_seconds',toc(started));
dmimo.json(fullfile(out,'audit.json'),audit);disp(audit);

 function emit(method,branch,R,ns,F,Et,C,db,ib,msg,position,overflow)
  cm=struct;cb=0;
  if ~strcmp(branch,'CPU'),cm=dmimo.encode_coeff(C,3,4);[C,cb]=dmimo.decode_coeff(cm);end
  m=dmimo.metrics(Et,C,sim.rho_data,e,ns);G=Et*C;checks=checks+1;
  if strcmp(mode,'smoke')&&strcmp(method,'POS_R3')
   legacy=readtable(fullfile(root,'results','dual_v1_20261022',array,'trials.csv'),'TextType','string');
   meth='POS_CPU';if strcmp(branch,'UE'),meth='POS_UE_Q';end
   old=legacy(legacy.stage=="MAIN"&legacy.B==8&legacy.observation_SNR_dB==30&legacy.R==3&legacy.sample_id==sid&legacy.channel==channel&legacy.streams==ns&legacy.method==meth,:);
   assert(height(old)==1);delta=abs(20*m.mmse-old.MMSE_Mbps);assert(delta<1e-8);maxRegression=max(maxRegression,delta);
  end
  assert(max(abs(sum(abs(F).^2,1)-1))<1e-10);
  addrow(method,branch,R,ns,db,ib,cb,db+ib+cb,m.joint,m.mmse,m.capacity,m.norm,m.ap_error,position,overflow,G,m.sinr,msg,cm);
 end
 function addrow(method,branch,R,ns,db,ib,cb,bits,joint,mmse,cap,power,aperr,position,overflow,G,sinr,msg,cm)
  ev=svd(G).^2;ratio=NaN;if ns==2,ratio=ev(2)/ev(1);end
  rows(end+1,:)={array,channel,sid,method,branch,R,ns,8,30,db,ib,cb,bits,20*mmse,20*joint,20*cap,power,aperr,...
   norm(position-s.position),overflow,ratio,min(sinr),jsonencode(sinr),jsonencode(position),jsonencode(msg),jsonencode(cm)}; %#ok<AGROW>
 end
end
function names=columns
names={'array','channel','sample_id','method','branch','R','streams','B','observation_SNR_dB','direction_bits','internal_bits','coefficient_bits','feedback_bits',...
 'MMSE_Mbps','joint_Mbps','rank_capacity_Mbps','precoder_norm_squared','per_AP_power_error','position_error_m','quantizer_overflow','received_mode_ratio','weakest_stream_SINR',...
 'stream_SINR','decoded_position','geometry_message','coefficient_message'};
end
function [joint,mmse,sinr]=rates(G,rho)
Q=eye(size(G,1))+rho*G*G';R=chol((Q+Q')/2);joint=2*sum(log2(real(diag(R))));sinr=zeros(1,size(G,2));
for j=1:size(G,2),other=G;other(:,j)=[];sinr(j)=max(0,real(rho*G(:,j)'*((eye(size(G,1))+rho*other*other')\G(:,j))));end
mmse=sum(log2(1+sinr));
end
