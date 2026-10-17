function summary=summarize(trials)
ids=unique(trials.case_id,'stable');rows=cell(numel(ids),25);
for j=1:numel(ids)
 t=trials(strcmp(trials.case_id,ids{j}),:);p=t.position_error_m;mask=isfinite(p);position=NaN;if any(mask),position=sqrt(mean(p(mask).^2));end
 n1=sum(t.AP1_phase_valid_count);np=sum(t.all_pair_phase_valid_count);phase1=NaN;phasepair=NaN;
 if n1>0,phase1=sqrt(sum(t.AP1_phase_squared_error_sum_rad2,'omitnan')/n1)*180/pi;end
 if np>0,phasepair=sqrt(sum(t.all_pair_phase_squared_error_sum_rad2,'omitnan')/np)*180/pi;end
 rows(j,:)={ids{j},t.array{1},t.geometry_id{1},t.channel_model{1},t.panel_branch{1},t.phase_source{1},t.method_id{1},t.B(1),t.R(1),t.observation_SNR_dB(1),t.reference_data_SNR_dB(1),height(t),height(t),mean(t.received_power),mean(t.power_gain),sqrt(mean(t.power_gain)),mean(t.SE_bps_Hz),mean(t.throughput_Mbps),position,phase1,phasepair,t.feedback_bits(1),sum(strcmp(t.localization_status,'nonfinite')),sum(strcmp(t.localization_status,'outside_ROI')),sum(t.outside_unit_disk_count,'omitnan')};
end
summary=cell2table(rows,'VariableNames',{'case_id','array','geometry_id','channel_model','panel_branch','phase_source','method_id','B','R','observation_SNR_dB','reference_data_SNR_dB','samples','valid_samples','mean_received_power','mean_power_gain','bf_gain_report','mean_SE_bps_Hz','mean_throughput_Mbps','position_RMSE_m','AP1_phase_RMSE_deg','all_pair_phase_RMSE_deg','feedback_bits','localization_failure_count','outside_ROI_count','outside_unit_disk_count'});
end
