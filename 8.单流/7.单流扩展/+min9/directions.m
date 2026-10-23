function [messages,frequency]=directions(net,Y,active,refs)
frequency=dmimo.frequency(net,Y,active,refs);d=1+strcmp(net.array,'UPA32x32');codes=zeros(d,numel(refs));
for j=1:numel(refs)
 slot=find(active==refs(j));ix=(slot-1)*2*net.Nap+(1:net.Nap);H=Y(:,[ix ix+net.Nap]);[V,~]=dmimo.modes(H,2);
 g=min9.spectrum(net,[V(1:net.Nap,:) V(net.Nap+1:end,:)]);[~,k]=max(sum(abs(g).^2,2));codes(:,j)=min9.codes(net,k);
end
q=14;L=2^q;
if d==1,k=mod(round(L*frequency),L);else,k=max(0,min(L-1,round((frequency+.5)*(L-1))));end
messages.DFT=struct('method','DFT','reference_ap_ids',refs,'integer_codes',codes,'q_frequency',NaN);
messages.HAQSE=struct('method','HAQSE','reference_ap_ids',refs,'integer_codes',k,'q_frequency',q);
end
