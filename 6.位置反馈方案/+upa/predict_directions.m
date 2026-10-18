function design=predict_directions(net,message,active_ids)
% Shared public decoder: recover position and final per-AP beams from directions.
allowed={'method','reference_ap_ids','integer_codes','q_frequency'};
if ~isempty(setdiff(fieldnames(message),allowed)),error('upa:UnexpectedMessageField','Direction message contains fields outside the agreed integer schema');end
assert(isequal(active_ids(1),1));refs=message.reference_ap_ids;codes=double(message.integer_codes);assert(all(ismember(refs,active_ids))&&size(codes,1)==2&&size(codes,2)==numel(refs)&&all(codes(:)==round(codes(:))));
if strcmp(message.method,'DFT')
 assert(all(codes(:)>=0&codes(:)<net.M));frequency=-(mod(codes/net.M+.5,1)-.5);direction_bits=numel(refs)*net.index_bits;
elseif strcmp(message.method,'HAQSE')
 L=2^message.q_frequency;assert(all(codes(:)>=0&codes(:)<L));frequency=-.5+codes/(L-1);direction_bits=numel(refs)*2*message.q_frequency;
else,error('upa:Method','Unknown direction message');end
[position,local]=upa.localize(net,refs,frequency);B=numel(active_ids);m=zeros(2,B);
for slot=1:B
 b=active_ids(slot);q=net.R(:,:,b)'*((net.pb(:,b)-position)/norm(net.pb(:,b)-position));m(:,slot)=mod(round(net.M*q(2:3)/2),net.M);
 reference=find(refs==b,1);
 if ~isempty(reference)
  if strcmp(message.method,'DFT'),m(:,slot)=codes(:,reference);else,m(:,slot)=mod(round(-frequency(:,reference)*net.M),net.M);end
 end
end
F=upa.beams(net,m);design=struct('method_id',sprintf('%s_R%d',message.method,numel(refs)),'phase_source','','F',F,'c',[],'w',[],'beam_indices',m,'position',position,'localization',local,'direction_message',message,'phase_indices',[],'internal_amp_indices',[],'internal_phase_indices',[],'direction_bits',direction_bits,'internal_bits',0,'AP_phase_bits',0,'feedback_bits',direction_bits,'reference_ap_ids',refs);
end
