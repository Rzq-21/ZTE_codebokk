function result=run_ula_sim(varargin)
% Official ULA64 CPU/UE simulation. Example: run_ula_sim('Mode','full','PanelBranch','CPU').
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);
parser=inputParser;addParameter(parser,'Mode','quick');addParameter(parser,'Samples',[]);addParameter(parser,'OutputDir','');addParameter(parser,'MakePlots',true);addParameter(parser,'PanelBranch','CPU');parse(parser,varargin{:});opts=parser.Results;
assert(ismember(opts.Mode,{'smoke','quick','full'}));assert(ismember(opts.PanelBranch,{'CPU','UE'}));[net,sim,bank]=ula.load_inputs;
if strcmp(opts.Mode,'full'),n=400;elseif strcmp(opts.Mode,'quick'),n=20;else,n=3;end
if ~isempty(opts.Samples),n=opts.Samples;end;validateattributes(n,{'numeric'},{'integer','scalar','>=',1,'<=',400});
if strcmp(opts.Mode,'smoke'),snrs=[20 30 40];counts=[5 8 10];else,snrs=sim.observation_SNR_dB(:).';counts=sim.AP_total_counts(:).';end
axes=unique([repmat(8,numel(snrs),1),snrs(:);counts(:),repmat(30,numel(counts),1)],'rows');channels={'LOS','LOS_NLOS_K10'};expected_rows=n*size(axes,1)*2*7;rows=cell(expected_rows,numel(upa.schema));row_number=0;
if isempty(opts.OutputDir),stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));outdir=fullfile(root,'results',sprintf('%s_%s_%d_%s',lower(opts.PanelBranch),opts.Mode,n,stamp));else,outdir=opts.OutputDir;end
if exist(fullfile(outdir,'results_trials.csv'),'file'),error('ula:ExistingOutput','Choose a new OutputDir to preserve existing results');end
if ~exist(outdir,'dir'),mkdir(outdir);end;diary(fullfile(outdir,'run.log'));cleanup=onCleanup(@()diary('off'));started=tic;
fprintf('ULA64 %s v2: mode=%s samples=%d settings=%d channels=2 methods=7\n',opts.PanelBranch,opts.Mode,n,size(axes,1));
for sample_id=1:n
 sample=ula.sample_channels(net,sim,bank,sample_id);
 for channel_index=1:2
  channel=channels{channel_index};
  for setting_index=1:size(axes,1)
   B=axes(setting_index,1);snr=axes(setting_index,2);active=net.active_order(1:B);[Htrue,Hobs,~,noise_variance]=ula.observe(sample,channel,active,snr);[truth.u,truth.s2]=upa.principal(Htrue);
   designs=ula.design_baselines(net,Hobs,Htrue,active,sim.UE_phase_bits);[messages,direction_diagnostic]=ula.encode_directions(net,Hobs,active,net.reference5,sim.frequency_bits_per_dimension);
   for refs_count=[3 5]
    if refs_count==3,refs=net.reference3;else,refs=net.reference5;end
    for scheme={'HAQSE','DFT'}
     msg=upa.subset_message(messages.(scheme{1}),refs);
     if strcmp(opts.PanelBranch,'CPU'),designs{end+1}=ula.cpu_from_message(net,msg,active);else,designs{end+1}=ula.ue_phase_from_message(net,msg,active,Hobs,sim.UE_phase_bits);end
    end
   end
   for method_index=1:7
    design=designs{method_index};metric=upa.evaluate(net,sim,Htrue,truth,design,sample.position);row_number=row_number+1;rows(row_number,:)=ula.result_row(net,sim,design,metric,channel,active,snr,sample_id,sample.position);
   end
   if sample_id==1&&B==8&&snr==30
    debug=struct('sample_id',sample_id,'channel',channel,'active_ids',active,'Htrue',Htrue,'Hobs',Hobs,'noise_variance',noise_variance,'messages',messages,'direction_diagnostic',direction_diagnostic,'designs',{designs});save(fullfile(outdir,['debug_sample1_' channel '.mat']),'debug','-v7');
   end
  end
 end
 fprintf('Completed sample %d/%d, elapsed %.1fs\n',sample_id,n,toc(started));
 if mod(sample_id,10)==0||sample_id==n,completed_samples=sample_id;save(fullfile(outdir,'checkpoint.mat'),'rows','row_number','completed_samples','opts','-v7');end
end
assert(row_number==expected_rows);trials=cell2table(rows,'VariableNames',upa.schema);assert(isequal(height(unique(trials(:,{'case_id','sample_id'}),'rows')),height(trials)));summary=upa.summarize(trials);assert(height(summary)==size(axes,1)*2*7);
if strcmp(opts.PanelBranch,'CPU'),branch_mask=~strcmp(bank.manifest.panel_branch,'UE');else,branch_mask=strcmp(bank.manifest.panel_branch,'UE')|strcmp(bank.manifest.panel_branch,'SHARED');end
expected=bank.manifest(strcmp(bank.manifest.array,'ULA64')&branch_mask,:);assert(all(ismember(summary.case_id,expected.case_id)),'ula:Manifest','Simulation cases differ from frozen manifest');
for j=1:height(summary),ref=expected(strcmp(expected.case_id,summary.case_id{j}),:);assert(height(ref)==1&&isequaln(summary.feedback_bits(j),ref.feedback_bits));end
writetable(trials,fullfile(outdir,'results_trials.csv'));writetable(summary,fullfile(outdir,'results_summary.csv'));budgets=bank.budget(strcmp(bank.budget.array,'ULA64'),:);writetable(budgets,fullfile(outdir,'feedback_budget.csv'));
metadata=struct('mode',opts.Mode,'panel_branch',opts.PanelBranch,'samples_per_case',n,'full_400_samples',n==400,'case_count',height(summary),'trial_count',height(trials),'elapsed_seconds',toc(started),'protocol_version','v2.0','geometry_id',net.geometry_id,'sample_bank_id',sim.sample_bank_id,'low_feedback_phase_source',ternary(strcmp(opts.PanelBranch,'CPU'),'CPU geometric channel reconstruction','UE noisy measured CSI with quantized phase feedback'),'complete_main_grid',size(axes,1)==10,'true_channel','exact spherical LOS or LOS+NLOS K10','observation_noise','global active H mean power','DFT_codebook_size',net.M,'DFT_oversampling',net.M/net.Nap,'frequency_bits',sim.frequency_bits_per_dimension,'UE_baseline_phase_bits',sim.UE_phase_bits,'reference_data_SNR_dB',sim.reference_data_SNR_dB,'reference_power',sim.reference_power,'bandwidth_Hz',sim.bandwidth_Hz);
fid=fopen(fullfile(outdir,'run_metadata.json'),'w');fwrite(fid,jsonencode(metadata),'char');fclose(fid);copyfile(fullfile(bank.shared_dir,'source_hashes.json'),fullfile(outdir,'input_hashes.json'));save(fullfile(outdir,'results.mat'),'trials','summary','metadata','net','sim','-v7');
if opts.MakePlots,plot_ula_sim(outdir,opts.PanelBranch);end
fprintf('COMPLETE %d cases %d rows: %s\n',height(summary),height(trials),outdir);result=struct('output_dir',outdir,'summary',summary,'metadata',metadata);
end
function out=ternary(condition,a,b),if condition,out=a;else,out=b;end,end
