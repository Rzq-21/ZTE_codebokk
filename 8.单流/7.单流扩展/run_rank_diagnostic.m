function run_rank_diagnostic
[root,base]=dmimo.setup;out=fullfile(root,'results','dual_v1_20261022');if ~exist(out,'dir'),mkdir(out);end
diary(fullfile(out,'rank.log'));cl=onCleanup(@()diary('off'));started=tic;
rows=cell(0,15);
for array={'ULA','UPA'}
 ar=array{1};pkg=lower(ar);loadfun=str2func([pkg '.load_inputs']);samplefun=str2func([pkg '.sample_channels']);[net,sim,bank]=loadfun();
 old=readtable(fullfile(base,'results','position_phase_v1_20261018',ar,'position_main_trials.csv'),'TextType','string');
 old=old(old.B==8&old.observation_SNR_dB==30&old.panel_branch=="CPU",:);active=net.active_order(1:8);
 for sid=1:400
  s=samplefun(net,sim,bank,sid);
  for ch={'LOS','MIX'}
   channel=ch{1};label=channel;if strcmp(channel,'MIX'),label='LOS_NLOS_K10';end
   H=reshape(s.(channel)(:,:,active),4,[]);[~,e]=dmimo.modes(H,2);
   for R=[3 5]
    t=old(old.sample_id==sid&old.channel_model==label&old.R==R,:);assert(height(t)==1);
    p=[t.estimated_x_m;t.estimated_y_m;t.estimated_z_m];pred=posfb.prediction(net,p,net.reference3,active);
    E=dmimo.reduce(net,H,pred.F,false);[~,a]=dmimo.modes(E,2);
    r1=log2(1+sim.rho_data*e(1));r2=sum(log2(1+sim.rho_data*e(1:2)/2));
    a1=log2(1+sim.rho_data*a(1));a2=sum(log2(1+sim.rho_data*a(1:2)/2));
    rows(end+1,:)={ar,channel,sid,R,e(2)/e(1),a(2)/a(1),sum(a)/sum(e),a(2)/e(2),r1,r2,dmimo.capacity(e,2,sim.rho_data),a1,a2,dmimo.capacity(a,2,sim.rho_data),sim.rho_data}; %#ok<AGROW>
   end
  end
  if mod(sid,50)==0,fprintf('RANK %s %d/400 elapsed %.1f s\n',ar,sid,toc(started));end
 end
end
names={'array','channel','sample_id','R','full_sigma2_ratio','beam_sigma2_ratio','beam_energy_retention','second_mode_energy_retention','full_R1_SE','full_R2_equal_SE','full_Rle2_capacity','beam_R1_SE','beam_R2_equal_SE','beam_Rle2_capacity','rho'};
writetable(cell2table(rows,'VariableNames',names),fullfile(out,'legacy_rank_trials.csv'));
dmimo.json(fullfile(out,'rank_audit.json'),struct('status','complete','records',size(rows,1),'samples_per_array',400,'elapsed_seconds',toc(started)));
end
