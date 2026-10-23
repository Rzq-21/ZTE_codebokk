function result=run_position_extension(array,varargin)
% Add POS to frozen v2 observations. CPU/UE share every coordinate and beam.
% run_position_extension('ULA','Mode','full','OutputDir',new_directory)
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));maxNumCompThreads(1);
p=inputParser;addParameter(p,'Mode','smoke');addParameter(p,'OutputDir','');addParameter(p,'Samples',[]);parse(p,varargin{:});opt=p.Results;
assert(ismember(array,{'ULA','UPA'})&&ismember(opt.Mode,{'smoke','pilot','full'}));
if strcmp(array,'ULA'),pkg='ula';else,pkg='upa';end
loadfun=str2func([pkg '.load_inputs']);samplefun=str2func([pkg '.sample_channels']);encodefun=str2func([pkg '.encode_directions']);localfun=str2func([pkg '.localize']);predictfun=str2func([pkg '.predict_directions']);
[net,sim,bank]=loadfun();cfg=posfb.config(net);channels={'LOS','LOS_NLOS_K10'};
if strcmp(opt.Mode,'full'),n=400;axes=unique([repmat(8,5,1),(20:5:40)';(5:10)',repmat(30,6,1)],'rows');
elseif strcmp(opt.Mode,'pilot')
 n=24;axes=[8 20;8 30;8 40];original_xy=bank.positions(1:2,:).';rng(20261018,'twister');bank.positions=[20+10*rand(2,n);zeros(1,n)];
 if strcmp(array,'UPA'),bank.positions(3,:)=5*rand(1,n);end
 bank.zeta=(randn(n,6)+1j*randn(n,6))/sqrt(2);bank.noise_seed=1700000+reshape(1:n*10,n,10);sim.sample_bank_id='POS_INDEPENDENT_PILOT_20261018';
 assert(isempty(intersect(bank.positions(1:2,:).',original_xy,'rows')));
else,n=3;axes=[5 30;8 20;8 30;8 40;10 30];end
if ~isempty(opt.Samples),n=opt.Samples;end
if isempty(opt.OutputDir),out=fullfile(root,'results',['position_phase_v1_' opt.Mode '_' array]);else,out=opt.OutputDir;end
if ~exist(out,'dir'),mkdir(out);end
cfg.mode=opt.Mode;cfg.samples=n;cfg.settings=axes;cfg.sample_bank_id=sim.sample_bank_id;
configfile=fullfile(out,'position_config.json');
if exist(configfile,'file'),prior=jsondecode(fileread(configfile));assert(isequal(prior.main_allocation(:),cfg.main_allocation(:))&&prior.samples==n&&strcmp(prior.sample_bank_id,cfg.sample_bank_id),'posfb:ConfigMismatch','Output config differs');else,writejson(configfile,cfg);end
save(fullfile(out,'public_config.mat'),'net','sim','cfg','-v7');
allocation_rows=cell(numel(cfg.budgets),5);
for j=1:numel(cfg.budgets),bb=posfb.allocate(cfg,cfg.budgets(j));allocation_rows(j,:)={cfg.budgets(j),jsonencode(bb),jsonencode((cfg.hi-cfg.lo)./2.^bb),cfg.id,cfg.domain_id};end
writetable(cell2table(allocation_rows,'VariableNames',{'coordinate_budget','coordinate_axis_bits','step_m','allocation_id','quantizer_domain_id'}),fullfile(out,'coordinate_allocations.csv'));
diary(fullfile(out,'run.log'));cleanup=onCleanup(@()diary('off'));started=tic;
blocksize=20;blockcount=ceil(n/blocksize);totalmain=0;totalsupp=0;
for block=1:blockcount
 first=(block-1)*blocksize+1;last=min(n,block*blocksize);mainfile=fullfile(out,sprintf('main_chunk_%03d.csv',block));suppfile=fullfile(out,sprintf('supplement_chunk_%03d.csv',block));
 if exist(mainfile,'file')&&exist(suppfile,'file'),fprintf('Reuse committed block %d\n',block);continue;end
 mainrows=cell(0,numel(posfb.columns));supprows=mainrows;
 for sid=first:last
  sample=samplefun(net,sim,bank,sid);
  for ci=1:2
   channel=channels{ci};
   for setting=1:size(axes,1)
    B=axes(setting,1);snr=axes(setting,2);active=net.active_order(1:B);
    [Htrue,Hobs]=upa.observe(sample,channel,active,snr);[truth.u,truth.s2]=upa.principal(Htrue);[~,~,vobs]=upa.principal(Hobs);
    [messages,diagnostic]=encodefun(net,Hobs,active,net.reference5,sim.frequency_bits_per_dimension);freq5=diagnostic.frequency_unquantized;
    for R=[3 5]
     if R==3,refs=net.reference3;else,refs=net.reference5;end
     [~,columns]=ismember(refs,net.reference5);frequency=freq5(:,columns);[pu,loc]=localfun(net,refs,frequency);
     [msg,encdiag]=posfb.encode(pu,refs,cfg,cfg.main_bits);prediction=posfb.decode(net,msg,active,cfg);
     extra=extras('MAIN',cfg.main_bits,msg.coordinate_bits,'POS','ray_unweighted',pu,loc,encdiag.overflow,cfg);
     pair=pairrows(prediction,extra);mainrows=[mainrows;pair]; %#ok<AGROW>
     % Low/high budget and continuous limits are also checked at pilot SNRs.
     dosupp=(B==8&&snr==30)||strcmp(opt.Mode,'pilot');
     if dosupp
      if strcmp(opt.Mode,'pilot'),budgets=unique([cfg.main_bits,max(cfg.budgets)]);else,budgets=cfg.budgets;end
      for C=budgets
       [msg,encdiag]=posfb.encode(pu,refs,cfg,C);prediction=posfb.decode(net,msg,active,cfg);
       extra=extras('BUDGET',C,msg.coordinate_bits,'POS','ray_unweighted',pu,loc,encdiag.overflow,cfg);supprows=[supprows;pairrows(prediction,extra)]; %#ok<AGROW>
      end
      prediction=posfb.prediction(net,pu,refs,active);prediction.method_id=sprintf('POS_CONT_R%d',R);
      extra=extras('BUDGET',NaN,[],'POS','ray_unweighted',pu,loc,false,cfg);extra.allocation='continuous_no_finite_message';supprows=[supprows;pairrows(prediction,extra)]; %#ok<AGROW>
      if ~strcmp(opt.Mode,'pilot')
       for q=4:16
        msg=messages.HAQSE;msg.reference_ap_ids=refs;msg.q_frequency=q;L=2^q;
        if strcmp(array,'ULA'),msg.integer_codes=mod(round(L*frequency),L);else,msg.integer_codes=max(0,min(L-1,round((frequency+.5)*(L-1))));end
        prediction=predictfun(net,msg,active);extra=extras('FREQUENCY',q,[],'FREQ','ray_unweighted',pu,loc,false,cfg);extra.allocation='original_frequency_grid';
        supprows=[supprows;pairrows(prediction,extra)]; %#ok<AGROW>
       end
       pw=posfb.weighted_position(net,pu,frequency,refs,Hobs,active);
       for C=[cfg.main_bits NaN]
        if isnan(C),prediction=posfb.prediction(net,pw,refs,active);prediction.method_id=sprintf('POS_WLS_CONT_R%d',R);bits=[];overflow=false;
        else,[msg,encdiag]=posfb.encode(pw,refs,cfg,C);prediction=posfb.decode(net,msg,active,cfg);prediction.method_id=sprintf('POS_WLS_R%d',R);bits=msg.coordinate_bits;overflow=encdiag.overflow;end
        extra=extras('WEIGHTED',C,bits,'POS','energy_weighted_frequency',pw,loc,overflow,cfg);
        if isnan(C),extra.allocation='continuous_no_finite_message';end
        supprows=[supprows;pairrows(prediction,extra)]; %#ok<AGROW>
       end
      end
     end
    end
   end
  end
  if mod(sid,5)==0||sid==last,fprintf('%s %s sample %d/%d elapsed %.1f s\n',array,opt.Mode,sid,n,toc(started));end
 end
 writechunk(mainrows,mainfile);writechunk(supprows,suppfile);totalmain=totalmain+size(mainrows,1);totalsupp=totalsupp+size(supprows,1);
 writejson(fullfile(out,'progress.json'),struct('completed_samples',last,'samples',n,'elapsed_seconds',toc(started),'status','running'));
end
main=readchunks(out,'main_chunk_*.csv');supp=readchunks(out,'supplement_chunk_*.csv');
assert(height(main)==n*size(axes,1)*8);assert(height(unique(main(:,{'case_id','sample_id'}),'rows'))==height(main));
assert(height(unique(supp(:,{'case_id','sample_id'}),'rows'))==height(supp));
writetable(main,fullfile(out,'position_main_trials.csv'));writetable(posfb.summarize(main),fullfile(out,'position_main_summary.csv'));
writetable(supp,fullfile(out,'supplement_trials.csv'));writetable(posfb.summarize(supp),fullfile(out,'supplement_summary.csv'));
audit=struct('array',array,'mode',opt.Mode,'samples',n,'settings',size(axes,1),'new_main_records',height(main),'supplement_records',height(supp),...
 'maximum_precoder_norm_error',max(abs(main.precoder_norm_squared-1)),'maximum_power_gain',max([main.power_gain;supp.power_gain]),...
 'main_overflow_count',sum(main.quantizer_overflow),'supplement_overflow_count',sum(supp.quantizer_overflow),...
 'main_axis_bits',cfg.main_allocation,'cpu_ue_coordinate_and_beam_equality','asserted per prediction',...
 'elapsed_seconds',toc(started),'status','complete');
assert(audit.maximum_precoder_norm_error<1e-10&&audit.maximum_power_gain<=1+1e-8);writejson(fullfile(out,'audit.json'),audit);
writejson(fullfile(out,'progress.json'),struct('completed_samples',n,'samples',n,'elapsed_seconds',toc(started),'status','complete'));
fprintf('COMPLETE %s main=%d supplement=%d elapsed %.1f s\n',array,height(main),height(supp),toc(started));result=audit;

 function rr=pairrows(pred,ex)
  [cpu,ue]=posfb.phases(net,pred,active,vobs,sim.UE_phase_bits);designs={cpu,ue};rr=cell(2,numel(posfb.columns));
  for branch=1:2
   dd=designs{branch};metric=upa.evaluate(net,sim,Htrue,truth,dd,sample.position);
   rr(branch,:)=posfb.result_row(net,sim,dd,metric,channel,active,snr,sid,sample.position,ex);
  end
 end
end
function ex=extras(experiment,C,bits,kind,localizer,pu,loc,overflow,cfg)
ex=struct('experiment',experiment,'budget',C,'bits',bits,'kind',kind,'localizer',localizer,'pu',pu,'condition',loc.condition,'overflow',overflow,'allocation',cfg.id,'domain',cfg.domain_id);
end
function writechunk(rows,path)
t=cell2table(rows,'VariableNames',posfb.columns);temporary=[path '.pending.csv'];writetable(t,temporary);movefile(temporary,path,'f');
end
function t=readchunks(out,pattern)
files=dir(fullfile(out,pattern));t=table;for i=1:numel(files),part=readtable(fullfile(out,files(i).name),'VariableNamingRule','preserve','TextType','char');t=[t;part];end
end
function writejson(path,value)
fid=fopen(path,'w');assert(fid>0);clean=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value),'char');
end
