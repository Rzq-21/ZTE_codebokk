function design=predict_directions(net,message,active_ids)
allowed={'method','reference_ap_ids','integer_codes','q_frequency'};
if ~isempty(setdiff(fieldnames(message),allowed)),error('ula:UnexpectedMessageField','Direction message contains fields outside the agreed integer schema');end
assert(isequal(active_ids(1),1));refs=message.reference_ap_ids;codes=double(message.integer_codes);assert(all(ismember(refs,active_ids))&&isvector(codes)&&numel(codes)==numel(refs)&&all(codes(:)==round(codes(:))));codes=reshape(codes,1,[]);
if strcmp(message.method,'DFT')
 assert(all(codes>=0&codes<net.M));frequency=mod(-codes/net.M+.5,1)-.5;direction_bits=numel(refs)*net.index_bits;
elseif strcmp(message.method,'HAQSE')
 L=2^message.q_frequency;assert(all(codes>=0&codes<L));frequency=mod(codes/L+.5,1)-.5;direction_bits=numel(refs)*message.q_frequency;
else,error('ula:Method','Unknown direction message');end
[position,local]=ula.localize(net,refs,frequency);B=numel(active_ids);m=zeros(1,B);
for slot=1:B
 b=active_ids(slot);q=net.R(:,:,b)'*((net.pb(:,b)-position)/norm(net.pb(:,b)-position));m(slot)=mod(round(net.M*q(2)/2),net.M);
 reference=find(refs==b,1);
 if ~isempty(reference)
  if strcmp(message.method,'DFT'),m(slot)=codes(reference);else,m(slot)=mod(round(-frequency(reference)*net.M),net.M);end
 end
end
F=ula.beams(net,m);design=struct('method_id',sprintf('%s_R%d',message.method,numel(refs)),'phase_source','','F',F,'c',[],'w',[],'beam_indices',m,'position',position,'localization',local,'direction_message',message,'phase_indices',[],'internal_amp_indices',[],'internal_phase_indices',[],'direction_bits',direction_bits,'internal_bits',0,'AP_phase_bits',0,'feedback_bits',direction_bits,'reference_ap_ids',refs);
end
