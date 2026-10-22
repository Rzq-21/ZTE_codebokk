function messages=baseline_messages(net,Y,B)
% Minimal extension of repository Type-I and composite two-beam Type-II.
% Both ranks use the same leading-two-mode spatial subspace for isolation.
[V,~]=dmimo.modes(Y,2);d=1+strcmp(net.array,'UPA32x32');one=zeros(d,B);two=zeros(d,2,B);amps=zeros(1,B);phases=zeros(1,B);
for b=1:B
 ix=(b-1)*2*net.Nap+(1:net.Nap);X=[V(ix,:) V(ix+net.Nap,:)];g=min9.spectrum(net,X);score=sum(abs(g).^2,2);
 [~,k]=max(score);one(:,b)=min9.codes(net,k);
 % Original groups contain orthogonal DFT beams. Keep two beams per AP.
 scores=score(net.groups);[values,order]=sort(scores,1,'descend');[~,group]=max(sum(values(1:2,:),1));chosen=net.groups(order(1:2,group),group);
 two(:,:,b)=min9.codes(net,chosen);P=g(chosen,:);Q=P*P';[A,D]=eig((Q+Q')/2,'vector');[~,j]=max(real(D));alpha=A(:,j);
 assert(abs(alpha(1))>1e-12);ratio=alpha(2)/alpha(1);assert(abs(ratio)<=1+1e-8);
 amps(b)=round(7*min(1,abs(ratio)));phases(b)=mod(round(angle(ratio)*4/(2*pi)),4);
 pair=min9.beams(net,two(:,:,b));assert(norm(pair'*pair-eye(2),'fro')<1e-10);
end
messages.TypeI=struct('version','BASE_DUAL_MIN1','kind','TypeI','codes',one,'amplitude_indices',[],'phase_indices',[]);
messages.TypeII=struct('version','BASE_DUAL_MIN1','kind','TypeII','codes',two,'amplitude_indices',amps,'phase_indices',phases);
end
