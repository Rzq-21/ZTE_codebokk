function row=result_row(net,sim,design,metric,channel,ids,snr,sample_id,p)
B=numel(ids);low=startsWith(design.method_id,'HAQSE_R')||startsWith(design.method_id,'DFT_R');cpu=low&&strcmp(design.phase_source,'CPU_model');ue=low&&strcmp(design.phase_source,'UE_estCSI');R=NaN;q=NaN;Q=NaN;dircodes=[];refs=[];phasecodes=[];
if low
 if cpu,panel='CPU';else,panel='UE';end
 refs=design.reference_ap_ids;R=numel(refs);dircodes=design.direction_message.integer_codes;
 if startsWith(design.method_id,'HAQSE'),q=sim.frequency_bits_per_dimension;end
 if ue,Q=sim.UE_phase_bits;phasecodes=design.phase_indices(2:end);end
else
 panel='SHARED';if ~strcmp(design.method_id,'SVD'),R=B;Q=sim.UE_phase_bits;dircodes=design.beam_indices;phasecodes=design.phase_indices(2:end);end
end
caseid=sprintf('V2_ULA64_%s_B%d_OBS%d_%s_%s',channel,B,snr,design.method_id,design.phase_source);status='not_applicable';if low,if design.localization.outside_user_region,status='outside_ROI';else,status='ok';end,end
phat=design.position;row={caseid,sim.protocol_version,sim.sample_bank_id,net.array,net.geometry_id,channel,panel,design.phase_source,design.method_id,B,R,q,Q,snr,sim.reference_data_SNR_dB,sample_id,jsonencode(ids),jsonencode(refs),jsonencode(dircodes),jsonencode(design.beam_indices),jsonencode(phasecodes),jsonencode(design.internal_amp_indices),jsonencode(design.internal_phase_indices),design.direction_bits,design.internal_bits,design.AP_phase_bits,design.feedback_bits,metric.norm_squared,metric.received_power,metric.ideal_power,metric.power_gain,metric.SE,metric.throughput,p(1),p(2),p(3),phat(1),phat(2),phat(3),metric.position_error,metric.AP1_SSE,metric.AP1_count,metric.pair_SSE,metric.pair_count,design.localization.condition,status,design.localization.outside_unit_disk_count};
assert(numel(row)==numel(upa.schema));
end
