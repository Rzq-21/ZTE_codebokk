function report=audit_ula_results
root=fileparts(mfilename('fullpath'));cpu=fullfile(root,'results','full400_ULA_CPU_20261017');ue=fullfile(root,'results','full400_ULA_UE_20261017');T1=readtable(fullfile(cpu,'results_trials.csv'),'VariableNamingRule','preserve');T2=readtable(fullfile(ue,'results_trials.csv'),'VariableNamingRule','preserve');S1=readtable(fullfile(cpu,'results_summary.csv'),'VariableNamingRule','preserve');S2=readtable(fullfile(ue,'results_summary.csv'),'VariableNamingRule','preserve');
assert(height(T1)==56000&&height(T2)==56000);assert(height(S1)==140&&height(S2)==140);assert(all(T1.sample_id>=1&T1.sample_id<=400)&&all(T2.sample_id>=1&T2.sample_id<=400));assert(all(abs(T1.precoder_norm_squared-1)<1e-10)&&all(abs(T2.precoder_norm_squared-1)<1e-10));assert(all(isfinite(T1.throughput_Mbps))&&all(isfinite(T2.throughput_Mbps)));
shared=ismember(T1.method_id,{'SVD','TypeI','TypeII'});key1=make_key(T1);key2=make_key(T2);[found,idx2]=ismember(key1,key2);assert(all(found));
max_shared=max(abs(T1.received_power(shared)-T2.received_power(idx2(shared))));
assert(max_shared<1e-13,'Shared CPU/UE baseline mismatch');low=~shared;ri=idx2(low);
max_direction=any(~strcmp(T1.direction_integer_codes(low),T2.direction_integer_codes(ri)));max_beam=any(~strcmp(T1.actual_beam_indices(low),T2.actual_beam_indices(ri)));
dx=T1.estimated_x_m(low)-T2.estimated_x_m(ri);dy=T1.estimated_y_m(low)-T2.estimated_y_m(ri);dz=T1.estimated_z_m(low)-T2.estimated_z_m(ri);max_position=max(sqrt(dx.^2+dy.^2+dz.^2));qphase=(T1.B(low)-1)*4;
assert(all(T1.feedback_bits(low)==T1.direction_bits(low)));assert(all(T2.feedback_bits(ri)==T2.direction_bits(ri)+qphase));assert(all(T2.AP_phase_bits(ri)==qphase));
assert(max_direction==0&&max_beam==0&&max_position<1e-12,'CPU/UE low-feedback direction or location mismatch');
assert(all(S1.samples==400)&&all(S2.samples==400));same_cost=(S1.feedback_bits==S2.feedback_bits)|(isnan(S1.feedback_bits)&isnan(S2.feedback_bits));assert(all(same_cost|contains(S1.method_id,'HAQSE_R')|contains(S1.method_id,'DFT_R'))); % costs differ only on low-feedback UE phase bits
channels={'LOS','LOS_NLOS_K10'};methods={'SVD','TypeII','TypeI','HAQSE_R3','DFT_R3','HAQSE_R5','DFT_R5'};out=cell(numel(channels)*numel(methods),16);n=0;
for c=1:numel(channels),for m=1:numel(methods)
 i=find(strcmp(S1.channel_model,channels{c})&S1.B==8&S1.observation_SNR_dB==30&strcmp(S1.method_id,methods{m}));j=find(strcmp(S2.channel_model,channels{c})&S2.B==8&S2.observation_SNR_dB==30&strcmp(S2.method_id,methods{m}));assert(numel(i)==1&&numel(j)==1);n=n+1;
 out(n,:)={channels{c},methods{m},S1.mean_throughput_Mbps(i),S2.mean_throughput_Mbps(j),S2.mean_throughput_Mbps(j)-S1.mean_throughput_Mbps(i),S1.bf_gain_report(i),S2.bf_gain_report(j),S2.bf_gain_report(j)-S1.bf_gain_report(i),S1.feedback_bits(i),S2.feedback_bits(j),S1.position_RMSE_m(i),S2.position_RMSE_m(j),S1.AP1_phase_RMSE_deg(i),S2.AP1_phase_RMSE_deg(j),S1.all_pair_phase_RMSE_deg(i),S2.all_pair_phase_RMSE_deg(j)};
end,end
comparison=cell2table(out,'VariableNames',{'channel_model','method_id','CPU_throughput_Mbps','UE_throughput_Mbps','UE_minus_CPU_throughput_Mbps','CPU_BF_gain','UE_BF_gain','UE_minus_CPU_BF_gain','CPU_feedback_bits','UE_feedback_bits','CPU_position_RMSE_m','UE_position_RMSE_m','CPU_AP1_phase_RMSE_deg','UE_AP1_phase_RMSE_deg','CPU_all_pair_phase_RMSE_deg','UE_all_pair_phase_RMSE_deg'});writetable(comparison,fullfile(root,'results','ULA_CPU_UE_B8_OBS30_comparison.csv'));
files=[dir(fullfile(cpu,'E*.png'));dir(fullfile(ue,'E*.png'))];assert(numel(files)==8,'Expected four report plots per ULA branch');
report=struct('passed',true,'ULA_CPU_trials',height(T1),'ULA_UE_trials',height(T2),'cases_per_branch',height(S1),'samples_per_case',400,'shared_baseline_max_received_power_difference',max_shared,'direction_integer_code_mismatches',max_direction,'beam_index_mismatches',max_beam,'max_CPU_UE_position_difference_m',max_position,'plot_count',numel(files),'comparison_csv',fullfile(root,'results','ULA_CPU_UE_B8_OBS30_comparison.csv'),'note','same 400 shared locations/noise and frozen v2 protocol');fid=fopen(fullfile(root,'results','ULA_full_audit.json'),'w');fwrite(fid,jsonencode(report),'char');fclose(fid);disp(report);disp(comparison);
end
function key=make_key(T)
key=strcat(string(T.sample_id),'|',string(T.B),'|',string(T.observation_SNR_dB),'|',string(T.channel_model),'|',string(T.method_id));
end
