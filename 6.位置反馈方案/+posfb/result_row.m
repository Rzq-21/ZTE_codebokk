function row=result_row(net,sim,design,metric,channel,ids,snr,sample_id,p,extra)
% Reuse frozen row schema through a direction-method alias, then label extension.
alias=design;alias.method_id=sprintf('HAQSE_R%d',numel(design.reference_ap_ids));
if strcmp(net.array,'ULA64'),row=ula.result_row(net,sim,alias,metric,channel,ids,snr,sample_id,p);else,row=upa.result_row(net,sim,alias,metric,channel,ids,snr,sample_id,p);end
row{9}=design.method_id;row{12}=NaN;
if strcmp(extra.kind,'FREQ'),row{12}=extra.budget;end
row{1}=sprintf('POSV1_%s_%s_B%d_OBS%d_%s_%s_%s_C%s_%s',net.array,channel,numel(ids),snr,design.method_id,design.phase_source,extra.experiment,num2str(extra.budget),extra.localizer);
row{2}='v2.1_POS_PHASE';pu=extra.pu;qe=norm(design.position-pu);
if strcmp(extra.kind,'FREQ'),qe=NaN;end
coordbudget=extra.budget;if strcmp(extra.kind,'FREQ'),coordbudget=NaN;end
row=[row,{extra.experiment,coordbudget,jsonencode(extra.bits),extra.allocation,extra.domain,extra.localizer,pu(1),pu(2),pu(3),norm(pu-p),qe,extra.overflow,extra.condition}];
end
