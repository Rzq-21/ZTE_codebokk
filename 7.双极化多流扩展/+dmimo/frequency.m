function f=frequency(net,Y,active,refs)
f=zeros(1+strcmp(net.array,'UPA32x32'),numel(refs));
for j=1:numel(refs)
 slot=find(active==refs(j));ix=(slot-1)*2*net.Nap+(1:2*net.Nap);H=Y(:,ix);
 S=[H(:,1:net.Nap).' H(:,net.Nap+1:end).'];
 [~,order]=sort(sum(abs(S).^2,1),'descend');S=S(:,order(1:4));
 if strcmp(net.array,'ULA64'),f(j)=haqse_ula_reference(S);
 else,z=MMV_CF_mQSE2_SumOfRatios(reshape(S,[32,32,4]),1,1,1,3);f(:,j)=mod(z(1,:).'+.5,1)-.5;end
end
end
