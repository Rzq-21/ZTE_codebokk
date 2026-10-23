function [messages,diagnostic]=encode_directions(net,Hobs,active_ids,refs,q_frequency)
codes=zeros(1,numel(refs));frequency=zeros(size(codes));
for j=1:numel(refs)
 slot=find(active_ids==refs(j));assert(numel(slot)==1);ix=(slot-1)*net.Nap+(1:net.Nap);Hb=Hobs(:,ix);[~,~,v]=upa.principal(Hb);codes(j)=ula.spectrum(net,v);
 frequency(j)=haqse_ula_reference(Hb.');
end
L=2^q_frequency;k=mod(round(L*frequency),L);
messages.DFT=struct('method','DFT','reference_ap_ids',refs(:).','integer_codes',codes,'q_frequency',NaN);messages.HAQSE=struct('method','HAQSE','reference_ap_ids',refs(:).','integer_codes',k,'q_frequency',q_frequency);
diagnostic=struct('frequency_unquantized',frequency,'HAQSE_frequency_quantized',mod(k/L+.5,1)-.5);
end
