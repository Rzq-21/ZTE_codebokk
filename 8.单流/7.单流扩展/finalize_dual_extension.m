function finalize_dual_extension
[root,~]=dmimo.setup;out=fullfile(root,'results','dual_v1_20261022');
for ar={'ULA','UPA'}
 array=ar{1};path=fullfile(out,array);T=readtable(fullfile(path,'trials.csv'),'TextType','string');
 keys={'array','stage','scenario','channel','B','observation_SNR_dB','R','streams','method','qA','qP'};
 [G,S]=findgroups(T(:,keys));S.samples=splitapply(@numel,T.MMSE_SE,G);
 for col={'feedback_bits','MMSE_Mbps','joint_Mbps','rank_capacity','position_error_m','full_sigma2_ratio'}
  name=col{1};S.(['mean_' name])=splitapply(@mean,T.(name),G);
 end
 S.MMSE_std_Mbps=splitapply(@std,T.MMSE_Mbps,G);S.MMSE_CI95_halfwidth_Mbps=1.96*S.MMSE_std_Mbps./sqrt(S.samples);
 assert(all(S.samples==400));writetable(S,fullfile(path,'summary.csv'));
 selected=S(S.stage=="MAIN"&S.B==8&S.observation_SNR_dB==30,:);writetable(selected,fullfile(path,'selected_B8_OBS30.csv'));
 methods=["POS_CPU" "POS_CPU" "POS_UE_Q" "POS_UE_Q" "FULLDIR_UE_Q" "CAP_RANK_LE_N"];
 refs=[3 5 3 5 0 0];labels={'POS CPU R3','POS CPU R5','POS UE R3','POS UE R5','Full direction UE','Capacity rank <= N'};
 colors=[.1 .3 .65;.35 .55 .85;.05 .55 .35;.35 .75 .55;.9 .4 .1;.35 .35 .35];
 for axisname={'SNR','B'}
  figure('Visible','off','Position',[100 100 1250 760],'Color','w');tl=tiledlayout(2,2,'TileSpacing','compact');
  for ci=1:2
   for ns=1:2
    nexttile;hold on;channel="LOS";if ci==2,channel="MIX";end
    for m=1:numel(methods)
     a=S(S.stage=="MAIN"&S.channel==channel&S.streams==ns&S.method==methods(m)&S.R==refs(m),:);
     if strcmp(axisname{1},'SNR'),a=a(a.B==8,:);x=a.observation_SNR_dB;else,a=a(a.observation_SNR_dB==30,:);x=a.B;end
     [x,ix]=sort(x);a=a(ix,:);style='-o';if m==6,style='--';end
     plot(x,a.mean_MMSE_Mbps,style,'Color',colors(m,:),'LineWidth',1.3,'MarkerSize',4);
    end
    grid on;ylabel('MMSE throughput (Mbps)');if strcmp(axisname{1},'SNR'),xlabel('CSI observation SNR (dB), B=8');else,xlabel('Cooperating AP count, observation SNR=30dB');end
    title(sprintf('%s | %s | %d stream(s)',array,channel,ns));if ci==1&&ns==1,legend(labels,'Location','best','FontSize',8);end
   end
  end
  title(tl,'Dual polarization: same hardware and total transmit power');saveplot(path,['rate_' axisname{1}]);
 end
 figure('Visible','off','Position',[100 100 1250 760],'Color','w');tl=tiledlayout(2,2,'TileSpacing','compact');
 for ci=1:2
  for ns=1:2
   channel="LOS";if ci==2,channel="MIX";end
   nexttile;hold on;v=zeros(4,3);scenarios=["nominal" "tilt20" "phase15" "combined"];
   for sc=1:4
    a=S(S.scenario==scenarios(sc)&S.channel==channel&S.streams==ns&S.R==3&S.B==8&S.observation_SNR_dB==30,:);
    v(sc,1)=a.mean_MMSE_Mbps(a.method=="POS_CPU");v(sc,2)=a.mean_MMSE_Mbps(a.method=="POS_UE_Q"&a.qA==3);
    if sc==1,v(sc,3)=v(sc,1);else,v(sc,3)=a.mean_MMSE_Mbps(a.method=="POS_CPU_CALIBRATED");end
   end
   bar(v);xticks(1:4);xticklabels({'Nominal','UE tilt 20deg','Port phase 15deg','Combined'});xtickangle(12);grid on;ylabel('MMSE throughput (Mbps)');title(sprintf('%s | %s | %d stream(s)',array,channel,ns));
   if ci==1&&ns==1,legend({'CPU assumes nominal','UE observed coefficients','CPU known calibration*'},'Location','best','FontSize',8);end
  end
 end
 title(tl,'R3, B=8, observation SNR=30dB; *extra calibration information, unpriced');saveplot(path,'polarization_mismatch');
 figure('Visible','off','Position',[100 100 1250 760],'Color','w');tl=tiledlayout(2,2,'TileSpacing','compact');
 for ci=1:2
  for ns=1:2
   channel="LOS";if ci==2,channel="MIX";end
   nexttile;hold on;
   for m=[3 4 5]
    a=S(S.scenario=="nominal"&S.channel==channel&S.streams==ns&S.method==methods(m)&S.R==refs(m)&S.B==8&S.observation_SNR_dB==30,:);
    [~,ix]=sort(a.mean_feedback_bits);a=a(ix,:);plot(a.mean_feedback_bits,a.mean_MMSE_Mbps,'-o','Color',colors(m,:),'LineWidth',1.4);
   end
   a=S(S.stage=="MAIN"&S.channel==channel&S.streams==ns&S.method=="POS_CPU"&S.R==3&S.B==8&S.observation_SNR_dB==30,:);
   plot(a.mean_feedback_bits,a.mean_MMSE_Mbps,'d','Color',colors(1,:),'MarkerSize',8,'LineWidth',1.4);
   grid on;ylabel('MMSE throughput (Mbps)');xlabel('Total feedback bits (coordinates/indices + coefficients)');title(sprintf('%s | %s | %d stream(s)',array,channel,ns));
   if ci==1&&ns==1,legend({'POS UE R3','POS UE R5','Full direction UE','POS CPU R3'},'Location','best','FontSize',8);end
  end
 end
 title(tl,'B=8, observation SNR=30dB; qA/qP = 2/3, 3/4, 4/5');saveplot(path,'feedback_rate');
 fprintf('Finalized %s: %d groups, 4 PNG + FIG pairs\n',array,height(S));
end
end
function saveplot(path,name)
exportgraphics(gcf,fullfile(path,[name '.png']),'Resolution',150);savefig(gcf,fullfile(path,[name '.fig']));close(gcf);
end
