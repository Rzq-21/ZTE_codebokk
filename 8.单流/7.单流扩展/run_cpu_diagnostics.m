function run_cpu_diagnostics
% Post-main diagnosis of a residual CPU gap. Oracles are not deployable curves.
[root,~]=dmimo.setup;c=dmimo.config;out=fullfile(root,'results','dual_v1_20261022');rows=cell(0,12);maxRegression=0;started=tic;
for ar={'ULA','UPA'}
 array=ar{1};loadfun=str2func([lower(array) '.load_inputs']);[net,sim,bank]=loadfun();active=net.active_order(1:8);
 T=readtable(fullfile(out,array,'trials.csv'),'TextType','string');
 T=T(T.stage=="MAIN"&T.B==8&T.observation_SNR_dB==30&T.method=="POS_CPU"&T.R==3,:);
 for sid=1:400
  s=dmimo.sample(net,sim,bank,sid,c);blocks=dmimo.channel(net,sim,s,c,'nominal');
  for ch={'LOS','MIX'}
   channel=ch{1};H=reshape(blocks.(channel)(:,:,active),8,[]);[~,e]=dmimo.modes(H,2);
   old=T(T.sample_id==sid&T.channel==channel,:);p=jsondecode(old.decoded_position(1));p=p(:);
   pred=posfb.prediction(net,p,net.reference3,active);Et=dmimo.reduce(net,H,pred.F);
   E_est=dmimo.cpu_model(net,p,active,pred.F,struct('E',c.E));
   E_true_p=dmimo.cpu_model(net,s.position,active,pred.F,struct('E',c.E));
   for ns=[1 2]
    sources={E_est,E_true_p,Et};labels={'CPU_estimated_p','ORACLE_true_p_same_beams','ORACLE_true_effective_H_same_beams'};
    for k=1:3
     C=dmimo.coefficients(sources{k},ns,8);m=dmimo.metrics(Et,C,sim.rho_data,e,ns);
     rows(end+1,:)={array,channel,sid,3,ns,labels{k},20*m.mmse,20*m.joint,m.norm,m.ap_error,norm(p-s.position),NaN}; %#ok<AGROW>
     if k==1
      delta=abs(20*m.mmse-old.MMSE_Mbps(old.streams==ns));assert(delta<1e-8);maxRegression=max(maxRegression,delta);
     end
    end
   end
  end
 end
 fprintf('CPU diagnosis %s complete elapsed %.1f s\n',array,toc(started));
end
names={'array','channel','sample_id','R','streams','variant','MMSE_Mbps','joint_Mbps','precoder_norm_squared','per_AP_power_error','estimated_position_error_m','oracle_finite_feedback_bits'};
writetable(cell2table(rows,'VariableNames',names),fullfile(out,'cpu_oracle_diagnostics.csv'));
dmimo.json(fullfile(out,'cpu_diagnostic_audit.json'),struct('status','passed','records',size(rows,1),'maximum_main_regression_Mbps',maxRegression,'elapsed_seconds',toc(started)));
end
