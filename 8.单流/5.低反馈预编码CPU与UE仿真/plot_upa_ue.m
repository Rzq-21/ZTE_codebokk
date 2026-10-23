function plot_upa_ue(outdir)
% UE-branch mean figures: throughput vs observation SNR, AP count, BF gain,
% and total uplink feedback. No confidence bands or sample filtering.
t=readtable(fullfile(outdir,'results_summary.csv'),'VariableNamingRule','preserve');
b=readtable(fullfile(outdir,'feedback_budget.csv'),'VariableNamingRule','preserve');
methods={'SVD','TypeII','TypeI','HAQSE_R3','DFT_R3','HAQSE_R5','DFT_R5'};
labels={'Ideal SVD','Type-II','Type-I','HAQSE R=3 (UE phase)','DFT R=3 (UE phase)','HAQSE R=5 (UE phase)','DFT R=5 (UE phase)'};
colors=[0 0 0;.15 .40 .80;0 .65 .70;.85 .20 .15;.65 .20 .70;.95 .55 .10;.15 .60 .30];
markers={'*','o','d','s','^','+','x'};channels={'LOS','LOS_NLOS_K10'};channel_titles={'LOS','LOS + NLOS, K=10 dB'};
settings={'E1_UE_throughput_SNR','E2_UE_throughput_AP','E1_UE_BF_gain_SNR'};
for kind=1:3
 f=figure('Visible','off','Color','w','Position',[100 100 1200 470]);cleanup=onCleanup(@()close(f));tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 for c=1:2
  nexttile;hold on;
  for j=1:numel(methods)
   mask=strcmp(t.channel_model,channels{c})&strcmp(t.method_id,methods{j});
   if kind==2,mask=mask&t.observation_SNR_dB==30;x=t.B(mask);else,mask=mask&t.B==8;x=t.observation_SNR_dB(mask);end
   if kind==3,y=t.bf_gain_report(mask);else,y=t.mean_throughput_Mbps(mask);end
   [x,order]=sort(x);plot(x,y(order),'-','Color',colors(j,:),'Marker',markers{j},'LineWidth',1.5,'MarkerSize',6,'DisplayName',labels{j});
  end
  grid on;box on;title(channel_titles{c});
  if kind==2,xlabel('Number of cooperating APs');xticks(5:10);else,xlabel('Observation SNR (dB)');xticks(20:5:40);end
  if kind==3,ylabel('Normalized BF gain: sqrt(mean(G))');ylim([0 1.03]);else,ylabel('Mean throughput (Mbps)');end
  if c==2,legend('Location','best','FontSize',8);end
 end
 title(tl,sprintf('UPA 32x32 | UE co-phase feedback | %d shared positions per point',t.samples(1)));
 exportgraphics(f,fullfile(outdir,[settings{kind} '.png']),'Resolution',180);savefig(f,fullfile(outdir,[settings{kind} '.fig']));clear cleanup;
end
f=figure('Visible','off','Color','w','Position',[100 100 750 470]);cleanup=onCleanup(@()close(f));hold on;
for j=2:numel(methods)
 mask=strcmp(b.method,methods{j});[x,order]=sort(b.B(mask));y=b.UE_total_bits(mask);plot(x,y(order),'-','Color',colors(j,:),'Marker',markers{j},'LineWidth',1.5,'DisplayName',labels{j});
end
grid on;box on;xticks(5:10);xlabel('Number of cooperating APs');ylabel('UE feedback bits per update');title('UPA 32x32 | UE branch total feedback budget');legend('Location','best','FontSize',9);
exportgraphics(f,fullfile(outdir,'E3_UE_feedback_AP.png'),'Resolution',180);savefig(f,fullfile(outdir,'E3_UE_feedback_AP.fig'));
end
