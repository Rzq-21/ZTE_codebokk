function audit = run_single_stream_task2(varargin)
%RUN_SINGLE_STREAM_TASK2 单流、400 MHz 任务二仿真入口。
%   默认 smoke 模式只使用少量固定样本检查数据链路与全部出图。
%   完整仿真示例：
%   run_single_stream_task2('Mode','full','OutputDir','results/single_stream_full');
%
%   CPU 侧只由固定 R=3 条方向/位置消息重构波束和共相位；UE 侧
%   使用观测 CSI 计算共相位。Type-I/II 是全反馈（UE-CSI）基线；
%   Docomo 方案仅在 UE 侧吞吐量和反馈开销图中比较。

root = fileparts(mfilename('fullpath'));
dep = fullfile(root,'6.位置反馈方案');
addpath(dep,fullfile(dep,'lib'));
parser = inputParser;
addParameter(parser,'Mode','smoke',@(x) ismember(x,{'smoke','full'}));
addParameter(parser,'OutputDir','',@(x) ischar(x) || isstring(x));
addParameter(parser,'Samples',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
addParameter(parser,'MakeFigures',true,@islogical);
parse(parser,varargin{:});
opt = parser.Results;

[net0,sim,bank] = ula.load_inputs();
sim.bandwidth_Hz = 400e6;
sim.rho_data = 10^(sim.reference_data_SNR_dB/10)/sim.reference_power;
sim.task_id = 'single_stream_400MHz_task2';
if isempty(opt.OutputDir)
    out = fullfile(root,'results',['single_stream_400MHz_' char(opt.Mode)]);
else
    out = char(opt.OutputDir);
end
if ~exist(out,'dir'), mkdir(out); end

if isempty(opt.Samples)
    if strcmp(opt.Mode,'full'), sampleIds = 1:400; else, sampleIds = 1:6; end
else
    sampleIds = 1:min(400,round(opt.Samples));
end
settings = make_settings(opt.Mode);
rows = cell(0,19);
started = tic;

for si = 1:numel(settings)
    q = settings(si);
    net = configure_ula(net0,q.Nant);
    ids = net.active_order(1:q.B);
    refs = net.reference3;
    for sid = sampleIds
        sample = make_sample(net,sim,bank,sid,q.KdB);
        [Htrue,Hobs] = observe_sample(sample,ids,q.SNRdB);
        [utrue,s2,vtrue] = upa.principal(Htrue);
        truth = struct('u',utrue,'s2',s2,'v',vtrue);
        [type1,type2,docomo,idealDesign] = full_feedback_design(net,Hobs,truth,q.B,sim.UE_phase_bits);
        [haqse,dftindex,position] = partial_designs(net,Hobs,ids,refs,sim);
        % Docomo 的码字选择以 UE 观测 CSI 为依据，按前次讨论不并入
        % CPU 几何反馈链路；Type-I/II 则作为完整 UE-CSI 反馈基线同时列出。
        designs = {idealDesign,type1,type2,haqse.cpu,dftindex.cpu,position.cpu, ...
                   idealDesign,type1,type2,docomo,haqse.ue,dftindex.ue,position.ue};
        branches = [repmat({'CPU'},1,6),repmat({'UE'},1,7)];
        for di = 1:numel(designs)
            design = designs{di};
            metric = evaluate_design(Htrue,truth,design,sim,sample.position);
            rows(end+1,:) = {branches{di},q.curve,q.xvalue,design.method_id,design.feedback_scope, ...
                design.cpu_estimation_dimension,q.B,q.Nant,q.SNRdB,q.KdB,sid, ...
                metric.position_error_m,metric.throughput_Mbps,design.feedback_bits, ...
                design.coordinate_bits_x,design.coordinate_bits_y,sample.position(1),sample.position(2), ...
                design.estimated_position(1)}; %#ok<AGROW>
        end
    end
    fprintf('单流 400 MHz: %s %d/%d，%.1f 秒\n',q.curve,si,numel(settings),toc(started));
end

names = {'branch','curve','x_value','method','feedback_scope','cpu_estimation_dimension', ...
    'B','antenna_count','observation_SNR_dB','Rician_K_dB','sample_id', ...
    'position_error_m','throughput_Mbps','feedback_bits','coordinate_bits_x','coordinate_bits_y', ...
    'true_x_m','true_y_m','estimated_x_m'};
trials = cell2table(rows,'VariableNames',names);
writetable(trials,fullfile(out,'single_stream_trials.csv'));
summary = summarize_trials(trials);
writetable(summary,fullfile(out,'single_stream_summary.csv'));

contract = struct('task','单流','bandwidth_Hz',sim.bandwidth_Hz,'streams',1, ...
    'samples',numel(sampleIds),'reference_AP_count',3,'partial_feedback_methods',{{'HAQSE','DFT-index','位置反馈'}}, ...
    'full_feedback_baselines',{{'Type-I','Type-II','Docomo 近场 LoS 码本（仅 UE）','完美 CSI SVD 上界'}}, ...
    'curves',{{'定位误差-SNR','定位误差-天线数','定位误差-AP数','吞吐量-SNR','吞吐量-天线数','吞吐量-莱斯因子','比特数-AP数'}}, ...
    'CPU_note','部分反馈方案的 CPU 信道估计维度固定为 R=3，不随协作 AP 数 B 增加。', ...
    'UE_note','UE 使用本地观测 CSI 计算 AP 间共相位；IDEAL_SVD 使用完美 CSI，仅作性能上界。', ...
    'docomo_note','Docomo 采用本地 4 倍过采样 DFT、QPSK beta、AP 逻辑栅格 DFT 和三轮交替更新。对 B=5:10 以 QPSK 坐标更新代替原文仅适用于小 B 的全枚举。', ...
    'status','complete','elapsed_seconds',toc(started));
write_json(fullfile(out,'experiment_contract.json'),contract);
if opt.MakeFigures
    plot_single_stream_task2(out);
end
audit = struct('status','complete','output_directory',out,'records',height(trials), ...
    'samples',numel(sampleIds),'settings',numel(settings),'bandwidth_Hz',sim.bandwidth_Hz, ...
    'elapsed_seconds',toc(started));
write_json(fullfile(out,'audit.json'),audit);
end

function settings = make_settings(mode)
if strcmp(mode,'full')
    snrs = 20:5:40; antennas = [16 32 64]; aps = 5:10; ks = 0:5:20;
else
    snrs = [20 30 40]; antennas = [16 32 64]; aps = [5 8 10]; ks = [0 10 20];
end
settings = struct('curve',{},'xvalue',{},'Nant',{},'B',{},'SNRdB',{},'KdB',{});
for x = snrs, settings(end+1) = item('SNR',x,64,8,x,10); end %#ok<AGROW>
for x = antennas, settings(end+1) = item('天线数',x,x,8,30,10); end %#ok<AGROW>
for x = aps, settings(end+1) = item('AP数',x,64,x,30,10); end %#ok<AGROW>
for x = ks, settings(end+1) = item('莱斯因子',x,64,8,30,x); end %#ok<AGROW>
end

function q = item(curve,x,Nant,B,SNRdB,KdB)
q = struct('curve',curve,'xvalue',x,'Nant',Nant,'B',B,'SNRdB',SNRdB,'KdB',KdB);
end

function net = configure_ula(net0,N)
net = net0;
net.Nap = N; net.M = 8*N; net.index_bits = round(log2(net.M));
net.local = [zeros(1,N);((0:N-1)-(N-1)/2)*net.d;zeros(1,N)];
net.ap_elements = zeros(3,N,10);
for b = 1:10
    net.ap_elements(:,:,b) = net.pb(:,b) + net.R(:,:,b)*net.local;
end
end

function sample = make_sample(net,sim,bank,sid,KdB)
p = bank.positions(:,sid); scatter = sim.scatterer_positions_m.';
sample.position = p; sample.H = zeros(4,net.Nap,10); sample.Z = sample.H;
K = 10^(KdB/10);
for b = 1:10
    ap = net.ap_elements(:,:,b); los = zeros(4,net.Nap); nlos = los;
    for u = 1:4
        user = p + net.ue_offsets(:,u); r = vecnorm(ap-user);
        los(u,:) = exp(-1j*net.k*r)./r;
        for ell = 1:6
            L = vecnorm(ap-scatter(:,ell)) + norm(scatter(:,ell)-user);
            nlos(u,:) = nlos(u,:) + bank.zeta(sid,ell)*exp(-1j*net.k*L)./L;
        end
    end
    nlos = nlos*norm(los,'fro')/norm(nlos,'fro');
    sample.H(:,:,b) = sqrt(K/(K+1))*los + sqrt(1/(K+1))*nlos;
    rng(bank.noise_seed(sid,b),'twister');
    sample.Z(:,:,b) = (randn(4,net.Nap)+1j*randn(4,net.Nap))/sqrt(2);
end
end

function [Htrue,Hobs] = observe_sample(sample,ids,snr)
Htrue = reshape(sample.H(:,:,ids),4,[]); Z = reshape(sample.Z(:,:,ids),4,[]);
noiseVariance = mean(abs(Htrue(:)).^2)/10^(snr/10);
Hobs = Htrue + sqrt(noiseVariance)*Z;
end

function [type1,type2,docomo,ideal] = full_feedback_design(net,Hobs,truth,B,Q)
% Type-I/II 与原仓库的 ULA 基线相同：本地 DFT 波束分别为 1/2 条，
% 后续使用量化 AP 间相位。全部索引和量化系数都计入 feedback_bits。
F1 = zeros(net.Nap,B);
F2 = zeros(net.Nap,B);
codebook = local_dft_codebook(net,net.M);
for slot = 1:B
    ix = (slot-1)*net.Nap + (1:net.Nap);
    Hb = Hobs(:,ix);
    [~,order] = sort(sum(abs(Hb*codebook).^2,1),'descend');
    F1(:,slot) = codebook(:,order(1));
    [~,~,v] = upa.principal(Hb);
    a1 = F1(:,slot)'*v; a2 = codebook(:,order(2))'*v;
    amp = round(7*min(1,abs(a2)/max(abs(a1),realmin)))/7;
    phase = 2*pi*mod(round(angle(a2*conj(a1))*4/(2*pi)),4)/4;
    f = F1(:,slot) + amp*exp(1j*phase)*codebook(:,order(2));
    F2(:,slot) = f/norm(f);
end
c1 = greedy_qpsk_phases(Hobs,F1,Q);
c2 = greedy_qpsk_phases(Hobs,F2,Q);
type1 = design('TYPE_I','全反馈（UE-CSI）',B*net.index_bits+(B-1)*Q,F1,c1,nan(3,1),B,NaN,NaN);
type2 = design('TYPE_II','全反馈（UE-CSI）',2*B*net.index_bits+5*B+(B-1)*Q,F2,c2,nan(3,1),B,NaN,NaN);
docomo = docomo_design(net,Hobs,B);
ideal = design('IDEAL_SVD','全反馈',NaN,[],[],nan(3,1),-1,NaN,NaN);
ideal.w = truth.v; ideal.feedback_scope = '全反馈';
end

function [haqse,dftindex,pos] = partial_designs(net,Hobs,ids,refs,sim)
R = numel(refs); fq = zeros(1,R); dft = zeros(1,R);
for j = 1:R
    slot = find(ids==refs(j),1); ix = (slot-1)*net.Nap+(1:net.Nap); Hb = Hobs(:,ix);
    raw = haqse_ula_reference(Hb.'); L = 2^sim.frequency_bits_per_dimension;
    code = mod(round(L*raw),L); fq(j) = mod(code/L+.5,1)-.5;
    [dft(j),~] = dft_beam(net,Hb);
end
df = mod(-dft/net.M+.5,1)-.5;
[pH,condH] = localize_ula(net,refs,fq); [pD,condD] = localize_ula(net,refs,df);
Fh = predict_beams(net,ids,refs,pH,mod(round(-fq*net.M),net.M));
Fd = predict_beams(net,ids,refs,pD,dft);
[pP,bits] = encode_position(net,pH,ids,24);
Fp = predict_beams(net,ids,[],pP,[]);
haqse = paired_design(net,Hobs,ids,Fh,pH,condH,'HAQSE_R3',R*sim.frequency_bits_per_dimension,R,sim.UE_phase_bits,NaN,NaN);
dftindex = paired_design(net,Hobs,ids,Fd,pD,condD,'DFT_INDEX_R3',R*net.index_bits,R,sim.UE_phase_bits,NaN,NaN);
pos = paired_design(net,Hobs,ids,Fp,pP,NaN,'POS_R3',sum(bits),2,sim.UE_phase_bits,bits(1),bits(2));
end

function pair = paired_design(net,Hobs,ids,F,p,condition,method,bits,dimension,Q,bx,by)
B = numel(ids); Hmodel = geometry_channel(net,p,ids); [~,~,vmodel] = upa.principal(Hmodel);
[~,~,vobs] = upa.principal(Hobs);
ccpu = co_phases(F,vmodel,net.Nap,Inf); cue = co_phases(F,vobs,net.Nap,Q);
pair.cpu = design(method,'部分反馈',bits,F,ccpu,p,dimension,bx,by);
pair.cpu.localization_condition = condition;
pair.ue = design(method,'部分反馈',bits+(B-1)*Q,F,cue,p,dimension,bx,by);
pair.ue.localization_condition = condition;
end

function d = design(method,scope,bits,F,c,p,dimension,bx,by)
d = struct('method_id',method,'feedback_scope',scope,'feedback_bits',bits, ...
    'F',F,'c',c,'w',[],'estimated_position',p,'cpu_estimation_dimension',dimension, ...
    'coordinate_bits_x',bx,'coordinate_bits_y',by,'localization_condition',NaN);
if ~isempty(F), d.w = reshape(F.*reshape(c,1,[]),[],1)/sqrt(size(F,2)); end
end

function [code,f] = dft_beam(net,Hb)
[~,~,v] = upa.principal(Hb); spectrum = fft(v,net.M)/sqrt(net.Nap);
[~,index] = max(abs(spectrum).^2); code = index-1;
n = (0:net.Nap-1).'; f = exp(1j*2*pi*n*code/net.M)/sqrt(net.Nap);
end

function C = local_dft_codebook(net,M)
% 所有本地 DFT 码字均为单位范数列向量；M/N 是过采样倍数。
n = (0:net.Nap-1).';
C = exp(1j*2*pi*n*(0:M-1)/M)/sqrt(net.Nap);
end

function c = greedy_qpsk_phases(H,F,Q)
% 与 Type-I 参考实现一致：固定第一个 AP，按 AP 次序贪婪地选量化相位。
B = size(F,2); Heff = zeros(size(H,1),B);
for slot = 1:B
    ix = (slot-1)*size(F,1)+(1:size(F,1));
    Heff(:,slot) = H(:,ix)*F(:,slot);
end
L = 2^Q; alphabet = exp(1j*2*pi*(0:L-1)/L); c = ones(B,1);
acc = Heff(:,1);
for slot = 2:B
    trials = acc + Heff(:,slot)*alphabet;
    [~,best] = max(sum(abs(trials).^2,1));
    c(slot) = alphabet(best); acc = trials(:,best);
end
end

function d = docomo_design(net,Hobs,B)
% Docomo Near-Field LoS 方案的 ULA 适配：
% 本地 4 倍过采样 DFT、QPSK beta、AP 逻辑栅格 DFT，以及 beta/v/f 交替更新。
% 原 simulate_proposed 对 beta 全枚举 4^(B-1)；此处 B=5:10，采用坐标更新。
Mlocal = 4*net.Nap; C = local_dft_codebook(net,Mlocal);
G = cell(B,1); index = ones(B,1);
for slot = 1:B
    ix = (slot-1)*net.Nap+(1:net.Nap); G{slot} = Hobs(:,ix)*C;
    [~,index(slot)] = max(sum(abs(G{slot}).^2,1));
end
V = ap_grid_dft_basis(B); beta = ones(B,1);
[~,best] = max(arrayfun(@(k) docomo_power(G,index,beta,V(:,k)),1:B));
v = V(:,best);
for outer = 1:3
    beta = docomo_beta_update(G,index,beta,v);
    [~,best] = max(arrayfun(@(k) docomo_power(G,index,beta,V(:,k)),1:B));
    v = V(:,best);
    total = docomo_effective_channel(G,index,beta,v);
    for slot = 1:B
        old = beta(slot)*v(slot)*G{slot}(:,index(slot));
        residual = total-old;
        trials = residual + beta(slot)*v(slot)*G{slot};
        [~,index(slot)] = max(sum(abs(trials).^2,1));
        total = residual + beta(slot)*v(slot)*G{slot}(:,index(slot));
    end
end
F = zeros(net.Nap,B);
for slot = 1:B, F(:,slot) = C(:,index(slot)); end
% design() 以 sqrt(B) 归一化，因此将 AP 栅格 DFT 的 1/sqrt(B) 还原到 c。
c = sqrt(B)*beta.*v;
bits = B*ceil(log2(Mlocal)) + 2*(B-1) + ceil(log2(B));
d = design('DOCOMO_NFLOS','全反馈（UE-CSI）',bits,F,c,nan(3,1),B,NaN,NaN);
end

function V = ap_grid_dft_basis(B)
% 当前活动 AP 的个数可变且不总是矩形网格，按 active_order 建立一维逻辑栅格。
n = (0:B-1).'; V = exp(1j*2*pi*n*(0:B-1)/B)/sqrt(B);
end

function power = docomo_power(G,index,beta,v)
h = docomo_effective_channel(G,index,beta,v); power = sum(abs(h).^2);
end

function h = docomo_effective_channel(G,index,beta,v)
h = zeros(size(G{1},1),1);
for slot = 1:numel(G)
    h = h + beta(slot)*v(slot)*G{slot}(:,index(slot));
end
end

function beta = docomo_beta_update(G,index,beta,v)
% beta(1)=1；其余 AP 依序在 QPSK 上做坐标最大化。
h = docomo_effective_channel(G,index,beta,v);
alphabet = exp(1j*pi/2*(0:3));
for slot = 2:numel(G)
    residual = h-beta(slot)*v(slot)*G{slot}(:,index(slot));
    trials = residual + v(slot)*G{slot}(:,index(slot))*alphabet;
    [~,best] = max(sum(abs(trials).^2,1));
    beta(slot) = alphabet(best); h = trials(:,best);
end
end

function F = predict_beams(net,ids,refs,p,referenceCodes)
B = numel(ids); F = zeros(net.Nap,B); n = (0:net.Nap-1).';
for slot = 1:B
    b = ids(slot); q = net.R(:,:,b)'*((net.pb(:,b)-p)/norm(net.pb(:,b)-p));
    m = mod(round(net.M*q(2)/2),net.M);
    if ~isempty(refs)
        ref = find(refs==b,1); if ~isempty(ref), m = referenceCodes(ref); end
    end
    F(:,slot) = exp(1j*2*pi*n*m/net.M)/sqrt(net.Nap);
end
end

function [p,condition] = localize_ula(net,refs,frequency)
A = zeros(2); rhs = zeros(2,1);
for j = 1:numel(refs)
    qy = -2*frequency(j); qx = -sqrt(max(0,1-qy*qy)); q = [qx;qy;0];
    g = -net.R(:,:,refs(j))*q; g = g(1:2)/norm(g(1:2)); P = eye(2)-g*g.';
    A = A+P; rhs = rhs+P*net.pb(1:2,refs(j));
end
condition = cond(A); p = [A\rhs;0];
end

function H = geometry_channel(net,p,ids)
H = zeros(4,numel(ids)*net.Nap);
for slot = 1:numel(ids)
    b = ids(slot); delta = net.pb(:,b)-p; r = norm(delta); q = delta/r; ql = net.R(:,:,b)'*q;
    ix = (slot-1)*net.Nap+(1:net.Nap);
    H(:,ix) = exp(-1j*net.k*r)/r*exp(1j*net.k*net.ue_offsets.'*q)*exp(-1j*net.k*ql.'*net.local);
end
end

function c = co_phases(F,v,N,Q)
B = size(F,2); a = zeros(B,1);
for slot = 1:B, ix = (slot-1)*N+(1:N); a(slot) = F(:,slot)'*v(ix); end
theta = angle(exp(1j*(angle(a)-angle(a(1)))));
if isinf(Q), c = exp(1j*theta); else, L = 2^Q; c = exp(1j*2*pi*mod(round(theta*L/(2*pi)),L)/L); end
end

function [pdecoded,bits] = encode_position(net,p,ids,totalBits)
bits = coordinate_bits(net,p,ids,totalBits); lo = [19;19]; hi = [31;31]; count = 2.^bits;
index = max(0,min(count-1,floor((p(1:2)-lo)./(hi-lo).*count)));
pdecoded = [lo+(index+.5).*(hi-lo)./count;0];
end

function bits = coordinate_bits(net,p,ids,totalBits)
lo = [19;19]; hi = [31;31]; weight = zeros(2,1);
for b = ids
    g = (net.pb(:,b)-p)/norm(net.pb(:,b)-p); weight = weight + g(1:2).^2;
end
best = Inf; bits = [4;totalBits-4];
for bx = 4:(totalBits-4)
    by = totalBits-bx; cost = sum(weight.*((hi-lo)./[2^bx;2^by]).^2);
    if cost < best, best = cost; bits = [bx;by]; end
end
end

function metric = evaluate_design(H,truth,d,sim,truePosition)
metric.position_error_m = NaN;
if all(isfinite(d.estimated_position)), metric.position_error_m = norm(d.estimated_position-truePosition); end
power = norm(H*d.w)^2; metric.throughput_Mbps = sim.bandwidth_Hz/1e6*log2(1+sim.rho_data*power);
assert(abs(norm(d.w)^2-1)<1e-9);
assert(power <= truth.s2*(1+1e-8));
end

function summary = summarize_trials(t)
keys = {'branch','curve','x_value','method','feedback_scope','cpu_estimation_dimension','B','antenna_count','observation_SNR_dB','Rician_K_dB'};
[g,summary] = findgroups(t(:,keys));
summary.samples = splitapply(@numel,t.sample_id,g);
summary.mean_position_error_m = splitapply(@(x) mean(x,'omitnan'),t.position_error_m,g);
summary.mean_throughput_Mbps = splitapply(@mean,t.throughput_Mbps,g);
summary.mean_feedback_bits = splitapply(@(x) mean(x,'omitnan'),t.feedback_bits,g);
end

function write_json(path,value)
fid = fopen(path,'w'); assert(fid>0);
fwrite(fid,jsonencode(value),'char'); fclose(fid);
end
