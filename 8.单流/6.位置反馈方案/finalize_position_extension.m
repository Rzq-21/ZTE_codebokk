function finalize_position_extension(array,out)
% Merge unchanged archived seven methods with the two new coordinate methods.
root=fileparts(mfilename('fullpath'));addpath(root);maxNumCompThreads(1);
main=readtable(fullfile(out,'position_main_trials.csv'),'VariableNamingRule','preserve');
supp=readtable(fullfile(out,'supplement_trials.csv'),'VariableNamingRule','preserve');
for branch={'CPU','UE'}
 br=branch{1};destination=fullfile(out,br);if ~exist(destination,'dir'),mkdir(destination);end
 oldpath=fullfile(root,'results',sprintf('full400_%s_%s_20261017',array,br),'results_trials.csv');
 old=readtable(oldpath,'VariableNamingRule','preserve');names=posfb.columns;
 for j=numel(upa.schema)+1:numel(names)
  name=names{j};if iscell(main.(name)),old.(name)=repmat({'legacy_v2'},height(old),1);else,old.(name)=nan(height(old),1);end
 end
 old.experiment_id=repmat({'MAIN'},height(old),1);position=main(strcmp(main.panel_branch,br),:);
 trials=[old;position];assert(height(trials)==72000&&height(unique(trials(:,{'case_id','sample_id'}),'rows'))==72000);
 summary=posfb.summarize(trials);assert(height(summary)==180&&all(summary.samples==400));
 writetable(trials,fullfile(destination,'results_trials.csv'));writetable(summary,fullfile(destination,'results_summary.csv'));
 budgets=cell(0,4);methods={'TypeII','TypeI','HAQSE_R3','DFT_R3','HAQSE_R5','DFT_R5','POS_R3','POS_R5'};
 for B=5:10,for j=1:numel(methods),t=summary(summary.B==B&strcmp(summary.method_id,methods{j}),:);assert(~isempty(t)&&numel(unique(t.feedback_bits))==1);budgets(end+1,:)={methods{j},B,t.feedback_bits(1),br};end,end %#ok<AGROW>
 budget=cell2table(budgets,'VariableNames',{'method','B','total_bits','branch'});writetable(budget,fullfile(destination,'feedback_budget.csv'));
 paired=pair_comparison(trials);writetable(paired,fullfile(destination,'paired_vs_TypeI.csv'));
 baseline=strcmp(trials.panel_branch,'SHARED');metadata=struct('array',array,'panel_branch',br,'cases',height(summary),'records',height(trials),...
  'old_records_reused',height(old),'new_position_records',height(position),'baseline_records',sum(baseline),'source_v2_trials',oldpath,...
  'coordinate_allocation','AP1 path-difference sensitivity; B8 geometry; frozen over B and both branches','main_method_count',9,'status','complete');
 writejson(fullfile(destination,'delivery_audit.json'),metadata);
 plot_main(summary,budget,destination,array,br);
end
ss=posfb.summarize(supp);plot_budget(ss,out,array);
threshold=thresholds(ss);writetable(threshold,fullfile(out,'one_percent_thresholds.csv'));
writejson(fullfile(out,'finalization_complete.json'),struct('array',array,'status','complete','PNG_figures',12,'MATLAB_FIG_figures',12));
fprintf('FINALIZED %s: CPU and UE each 72000 records, main figures and budget figures complete\n',array);
end
function t=pair_comparison(trials)
rows=cell(0,12);ids=unique(trials.case_id,'stable');
for j=1:numel(ids)
 t=sortrows(trials(strcmp(trials.case_id,ids{j}),:),'sample_id');if ~all(isfinite(t.feedback_bits)),continue;end
 base=sortrows(trials(trials.B==t.B(1)&trials.observation_SNR_dB==t.observation_SNR_dB(1)&strcmp(trials.channel_model,t.channel_model{1})&strcmp(trials.method_id,'TypeI'),:),'sample_id');
 assert(isequal(t.sample_id,base.sample_id));delta=t.throughput_Mbps-base.throughput_Mbps;m=mean(delta);ci=1.96*std(delta)/sqrt(height(t));
 rows(end+1,:)={t.case_id{1},t.channel_model{1},t.method_id{1},t.B(1),t.observation_SNR_dB(1),t.feedback_bits(1),...
  100*(1-t.feedback_bits(1)/base.feedback_bits(1)),100*(1-mean(t.throughput_Mbps)/mean(base.throughput_Mbps)),m,m-ci,m+ci,prctile(delta,5)}; %#ok<AGROW>
end
t=cell2table(rows,'VariableNames',{'case_id','channel_model','method_id','B','observation_SNR_dB','feedback_bits','feedback_saving_pct',...
 'throughput_loss_pct','paired_delta_Mbps','paired_delta_CI95_low_Mbps','paired_delta_CI95_high_Mbps','paired_delta_p05_Mbps'});
end
function result=thresholds(t)
rows=cell(0,9);
for ch={'LOS','LOS_NLOS_K10'},for br={'CPU','UE'},for R=[3 5]
 mask=strcmp(t.channel_model,ch{1})&strcmp(t.panel_branch,br{1})&strcmp(t.experiment_id,'BUDGET');
 cont=t(mask&strcmp(t.method_id,sprintf('POS_CONT_R%d',R)),:);assert(height(cont)==1);
 q=t(mask&strcmp(t.method_id,sprintf('POS_R%d',R)),:);loss=100*(1-q.mean_throughput_Mbps/cont.mean_throughput_Mbps);valid=find(loss<=1);
 if isempty(valid),C=NaN;total=NaN;actual=NaN;else,[C,ix]=min(q.coordinate_budget(valid));k=valid(ix);total=q.feedback_bits(k);actual=loss(k);end
 rows(end+1,:)={ch{1},br{1},R,C,total,actual,cont.mean_throughput_Mbps,min(q.coordinate_budget),max(q.coordinate_budget)}; %#ok<AGROW>
end,end,end
result=cell2table(rows,'VariableNames',{'channel_model','branch','R','minimum_tested_coordinate_bits','total_feedback_bits','loss_to_continuous_pct',...
 'continuous_throughput_Mbps','minimum_scanned_bits','maximum_scanned_bits'});
end
function plot_main(t,b,out,array,branch)
methods={'SVD','TypeII','TypeI','HAQSE_R3','DFT_R3','HAQSE_R5','DFT_R5','POS_R3','POS_R5'};
labels={'Ideal SVD','Type-II','Type-I','HAQSE R=3','DFT R=3','HAQSE R=5','DFT R=5','Position R=3','Position R=5'};
if strcmp(branch,'CPU'),labels{2}='Type-II (UE phase)';labels{3}='Type-I (UE phase)';end
colors=[0 0 0;.50 .50 .50;0 .55 .55;.85 .25 .15;.60 .25 .70;.95 .60 .15;.25 .60 .30;.10 .30 .85;.10 .65 .95];
markers={'*','o','d','s','^','+','x','p','h'};styles={':','--','--','-','-','-','-','-','-'};
channels={'LOS','LOS_NLOS_K10'};titles={'LOS','LOS + NLOS (K = 10 dB)'};
files={'E1_throughput_SNR','E2_throughput_AP','E1_BF_gain_SNR'};
for kind=1:3
 f=figure('Visible','off','Color','w','Position',[50 50 1400 650]);tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
 for c=1:2
  ax=nexttile;hold(ax,'on');
  for j=1:numel(methods)
   mask=strcmp(t.channel_model,channels{c})&strcmp(t.method_id,methods{j});
   if kind==2,mask=mask&t.observation_SNR_dB==30;x=t.B(mask);else,mask=mask&t.B==8;x=t.observation_SNR_dB(mask);end
   if kind==3,y=t.bf_gain_report(mask);else,y=t.mean_throughput_Mbps(mask);end
   [x,ix]=sort(x);plot(ax,x,y(ix),'Color',colors(j,:),'LineStyle',styles{j},'Marker',markers{j},'LineWidth',1.6,'MarkerSize',6,'DisplayName',labels{j});
  end
  set(ax,'FontSize',11,'LineWidth',.8);grid(ax,'on');box(ax,'on');title(ax,titles{c});
  if kind==2,xlabel(ax,'Number of cooperating APs');xticks(ax,5:10);else,xlabel(ax,'Observation SNR (dB)');xticks(ax,20:5:40);end
  if kind==3,ylabel(ax,'Normalized BF gain: sqrt(mean(G))');ylim(ax,[0 1.03]);else,ylabel(ax,'Mean throughput (Mbps)');end
 end
 lg=legend('NumColumns',3,'FontSize',10);lg.Layout.Tile='south';
 title(tl,sprintf('%s | %s co-phase | 400 paired positions | fixed total power',array,branch),'FontSize',14);
 exportgraphics(f,fullfile(out,[files{kind} '.png']),'Resolution',180);savefig(f,fullfile(out,[files{kind} '.fig']));close(f);
end
f=figure('Visible','off','Color','w','Position',[50 50 1000 650]);hold on;
for j=2:numel(methods)
 mask=strcmp(b.method,methods{j});[x,ix]=sort(b.B(mask));y=b.total_bits(mask);
 plot(x,y(ix),'Color',colors(j,:),'LineStyle',styles{j},'Marker',markers{j},'LineWidth',1.6,'MarkerSize',6,'DisplayName',labels{j});
end
grid on;box on;xticks(5:10);xlabel('Number of cooperating APs');ylabel('Feedback bits per update');set(gca,'FontSize',11);
title(sprintf('%s | %s branch | fixed coordinate budget',array,branch));legend('Location','northwest','NumColumns',2);
exportgraphics(f,fullfile(out,'E3_feedback_AP.png'),'Resolution',180);savefig(f,fullfile(out,'E3_feedback_AP.fig'));close(f);
end
function plot_budget(t,out,array)
channels={'LOS','LOS_NLOS_K10'};titles={'LOS','LOS + NLOS (K = 10 dB)'};
methods={'POS_R3','POS_R5','HAQSE_R3','HAQSE_R5'};labels={'Position R=3','Position R=5','HAQSE freq. R=3','HAQSE freq. R=5'};
colors=[.10 .30 .85;.10 .65 .95;.85 .25 .15;.95 .60 .15];markers={'p','h','s','+'};
for branch={'CPU','UE'}
 br=branch{1};
 for kind=1:2
  f=figure('Visible','off','Color','w','Position',[50 50 1400 650]);tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
  for c=1:2
   ax=nexttile;hold(ax,'on');
   for j=1:4
    mask=strcmp(t.channel_model,channels{c})&strcmp(t.panel_branch,br)&strcmp(t.method_id,methods{j});
    if j<=2,mask=mask&strcmp(t.experiment_id,'BUDGET');else,mask=mask&strcmp(t.experiment_id,'FREQUENCY');end
    q=sortrows(t(mask,:),'feedback_bits');if kind==1,y=q.mean_throughput_Mbps;else,y=q.mean_power_gain;end
    plot(ax,q.feedback_bits,y,'-','Color',colors(j,:),'Marker',markers{j},'LineWidth',1.6,'MarkerSize',6,'DisplayName',labels{j});
   end
   for R=[3 5]
    q=t(strcmp(t.channel_model,channels{c})&strcmp(t.panel_branch,br)&strcmp(t.method_id,sprintf('POS_CONT_R%d',R)),:);
    if kind==1,y=q.mean_throughput_Mbps;else,y=q.mean_power_gain;end
    col=colors(1+(R==5),:);yline(ax,y,':','Color',col,'LineWidth',1.2,'DisplayName',sprintf('Continuous position R=%d',R));
   end
   old=readtable(fullfile(out,br,'results_summary.csv'),'VariableNamingRule','preserve');
   base=old(old.B==8&old.observation_SNR_dB==30&strcmp(old.channel_model,channels{c})&strcmp(old.method_id,'TypeI'),:);
   if kind==1,y=base.mean_throughput_Mbps;else,y=base.mean_power_gain;end
   plot(ax,base.feedback_bits,y,'kd','MarkerFaceColor','k','MarkerSize',7,'DisplayName','Type-I (UE phase)');
   grid(ax,'on');box(ax,'on');title(ax,titles{c});xlabel(ax,'Total feedback bits per update');set(ax,'FontSize',11);
   if kind==1,ylabel(ax,'Mean throughput (Mbps)');else,ylabel(ax,'Mean normalized power gain: mean(G)');ylim(ax,[0 1.03]);end
  end
  lg=legend('NumColumns',3,'FontSize',10);lg.Layout.Tile='south';title(tl,sprintf('%s | %s co-phase | B=8 | observation SNR=30 dB',array,br),'FontSize',14);
  if kind==1,name=['E4_' br '_throughput_bits'];else,name=['E4_' br '_power_gain_bits'];end
  exportgraphics(f,fullfile(out,[name '.png']),'Resolution',180);savefig(f,fullfile(out,[name '.fig']));close(f);
 end
end
end
function writejson(path,value)
fid=fopen(path,'w');assert(fid>0);clean=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value),'char');
end
