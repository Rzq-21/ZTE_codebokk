function [messages,diagnostic]=encode_directions(net,Hobs,active_ids,refs,q_frequency)
% UE only. Output messages contain integer fields, not measured CSI.
codes=zeros(2,numel(refs));frequency=zeros(size(codes));
for j=1:numel(refs)
 slot=find(active_ids==refs(j));assert(numel(slot)==1);ix=(slot-1)*net.Nap+(1:net.Nap);Hb=Hobs(:,ix);[~,~,v]=upa.principal(Hb);codes(:,j)=upa.spectrum(net,v);
 f=MMV_CF_mQSE2_SumOfRatios(reshape(Hb.',[32,32,4]),1,1,1,3);frequency(:,j)=mod(f(1,:).'+.5,1)-.5;
end
L=2^q_frequency;k=max(0,min(L-1,round((frequency+.5)*(L-1))));
messages.DFT=struct('method','DFT','reference_ap_ids',refs(:).','integer_codes',codes,'q_frequency',NaN);messages.HAQSE=struct('method','HAQSE','reference_ap_ids',refs(:).','integer_codes',k,'q_frequency',q_frequency);
diagnostic=struct('frequency_unquantized',frequency,'HAQSE_frequency_quantized',-.5+k/(L-1));
end
