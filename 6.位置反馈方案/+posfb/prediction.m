function design=prediction(net,p,refs,active)
% All beams, including reference APs, are generated from the decoded position.
B=numel(active);d=1+strcmp(net.array,'UPA32x32');m=zeros(d,B);
for slot=1:B
 b=active(slot);delta=net.pb(:,b)-p;assert(norm(delta)>0);q=net.R(:,:,b)'*(delta/norm(delta));
 m(:,slot)=mod(round(net.M*q(2:d+1)/2),net.M);
end
if d==1,F=ula.beams(net,m);outside=any(p(1:2)<20|p(1:2)>30);else,F=upa.beams(net,m);outside=any(p<[20;20;0]|p>[30;30;5]);end
design=struct('method_id',sprintf('POS_R%d',numel(refs)),'phase_source','','F',F,'c',[],'w',[],...
 'beam_indices',m,'position',p,'localization',struct('condition',NaN,'outside_unit_disk_count',NaN,'outside_user_region',outside),...
 'direction_message',struct('integer_codes',[]),'phase_indices',[],'internal_amp_indices',[],'internal_phase_indices',[],...
 'direction_bits',NaN,'internal_bits',0,'AP_phase_bits',0,'feedback_bits',NaN,'reference_ap_ids',refs);
end
