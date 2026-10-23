function plot_single_stream_task2(resultDir)
%PLOT_SINGLE_STREAM_TASK2 集中生成任务二的所有图；不设置任何图片标题。
if nargin<1 || isempty(resultDir), error('请提供 run_single_stream_task2 的结果目录。'); end
T = readtable(fullfile(resultDir,'single_stream_trials.csv'),'TextType','string');
for branch = ["CPU" "UE"]
    out = fullfile(resultDir,char(branch)); if ~exist(out,'dir'), mkdir(out); end
    tb = T(T.branch==branch,:);
    make_position_plot(tb,branch,'SNR','观测 SNR（dB）','定位误差（m）',out);
    make_position_plot(tb,branch,'天线数','每 AP 天线数','定位误差（m）',out);
    make_position_plot(tb,branch,'AP数','协作 AP 数','定位误差（m）',out);
    make_rate_plot(tb,branch,'SNR','观测 SNR（dB）',out);
    make_rate_plot(tb,branch,'天线数','每 AP 天线数',out);
    make_rate_plot(tb,branch,'莱斯因子','莱斯因子 K（dB）',out);
    make_bits_plot(tb,branch,out);
end
end

function make_position_plot(T,branch,curve,xlabelText,ylabelText,out)
T = T(T.curve==curve & isfinite(T.position_error_m),:); f = new_figure(); ax = axes(f); hold(ax,'on');
for method = ["HAQSE_R3" "DFT_INDEX_R3" "POS_R3"]
    q = aggregate(T(T.method==method,:),{'x_value'},'position_error_m');
    plot(ax,q.x_value,q.mean_value,'-o','Color',method_color(method),'LineWidth',2.0,'MarkerSize',7,'DisplayName',method_label(method));
end
format_axes(ax,xlabelText,ylabelText); legend(ax,'Location','best','FontSize',13);
save_figure(f,out,[char(branch) '_定位误差_' char(curve)]);
end

function make_rate_plot(T,branch,curve,xlabelText,out)
T = T(T.curve==curve,:); f = new_figure(); ax = axes(f); hold(ax,'on');
if branch=="CPU"
    methods = ["IDEAL_SVD" "TYPE_I" "TYPE_II" "HAQSE_R3" "DFT_INDEX_R3" "POS_R3"];
else
    methods = ["IDEAL_SVD" "TYPE_I" "TYPE_II" "DOCOMO_NFLOS" "HAQSE_R3" "DFT_INDEX_R3" "POS_R3"];
end
for method = methods
    q = aggregate(T(T.method==method,:),{'x_value'},'throughput_Mbps');
    plot(ax,q.x_value,q.mean_value,'-o','Color',method_color(method),'LineWidth',2.0,'MarkerSize',7,'DisplayName',method_label(method));
end
format_axes(ax,xlabelText,'吞吐量（Mbps）'); legend(ax,'Location','best','FontSize',12);
save_figure(f,out,[char(branch) '_吞吐量_' char(curve)]);
end

function make_bits_plot(T,branch,out)
T = T(T.curve=="AP数" & T.method~="IDEAL_SVD",:); f = new_figure(); ax = axes(f); hold(ax,'on');
if branch=="CPU"
    methods = ["TYPE_I" "TYPE_II" "HAQSE_R3" "DFT_INDEX_R3" "POS_R3"];
else
    methods = ["TYPE_I" "TYPE_II" "DOCOMO_NFLOS" "HAQSE_R3" "DFT_INDEX_R3" "POS_R3"];
end
for method = methods
    q = aggregate(T(T.method==method,:),{'x_value'},'feedback_bits');
    plot(ax,q.x_value,q.mean_value,'-o','Color',method_color(method),'LineWidth',2.0,'MarkerSize',7,'DisplayName',method_label(method));
end
format_axes(ax,'协作 AP 数','每次更新反馈比特数'); legend(ax,'Location','best','FontSize',12);
save_figure(f,out,[char(branch) '_比特数_AP数']);
end

function q = aggregate(T,keys,field)
[g,q] = findgroups(T(:,keys)); q.mean_value = splitapply(@mean,T.(field),g); q = sortrows(q,keys);
end

function c = method_color(method)
switch method
    case "IDEAL_SVD", c = [0 0 0];
    case "TYPE_I", c = [.80 .24 .14];
    case "TYPE_II", c = [.91 .55 .10];
    case "DOCOMO_NFLOS", c = [.12 .62 .63];
    case "HAQSE_R3", c = [.05 .38 .75];
    case "DFT_INDEX_R3", c = [.10 .60 .28];
    case "POS_R3", c = [.60 .25 .72];
end
end

function label = method_label(method)
switch method
    case "IDEAL_SVD", label = '全反馈-完美 CSI（上界）';
    case "TYPE_I", label = '全反馈-Type-I（UE-CSI）';
    case "TYPE_II", label = '全反馈-Type-II（UE-CSI）';
    case "DOCOMO_NFLOS", label = 'Docomo 近场 LoS（UE-CSI）';
    case "HAQSE_R3", label = '部分反馈-HAQSE（R=3）';
    case "DFT_INDEX_R3", label = '部分反馈-DFT-index（R=3）';
    case "POS_R3", label = '部分反馈-位置反馈（R=3）';
end
end

function f = new_figure
f = figure('Visible','off','Color','w','Position',[100 100 1060 700]);
end

function format_axes(ax,xlabelText,ylabelText)
xlabel(ax,xlabelText); ylabel(ax,ylabelText); grid(ax,'on'); box(ax,'on');
set(ax,'FontSize',16,'LineWidth',1.1);
end

function save_figure(f,out,name)
exportgraphics(f,fullfile(out,[name '.png']),'Resolution',220);
savefig(f,fullfile(out,[name '.fig'])); close(f);
end
