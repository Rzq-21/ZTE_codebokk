function p=weighted_position(net,p,frequency,refs,Hobs,active)
% Separate UE-only ablation: observed-energy weighted frequency least squares.
weights=zeros(1,numel(refs));
for j=1:numel(refs),slot=find(active==refs(j));Hb=Hobs(:,(slot-1)*net.Nap+(1:net.Nap));weights(j)=norm(Hb,'fro')^2;end
d=size(frequency,1);weights=repelem(weights/mean(weights),d).';target=frequency(:);dim=2+strcmp(net.array,'UPA32x32');
for it=1:20
 [f,J]=forward(net,p,refs,d,dim);e=mod(f-target+.5,1)-.5;loss=sum(weights.*e.^2);
 step=-(J'*(weights.*J)+1e-12*eye(dim))\(J'*(weights.*e));
 if norm(step)<1e-8,break;end
 accepted=false;
 for alpha=[1 .5 .25 .125 .0625]
  q=p;q(1:dim)=q(1:dim)+alpha*step;fc=forward(net,q,refs,d,dim);ec=mod(fc-target+.5,1)-.5;
  if sum(weights.*ec.^2)<loss,p=q;accepted=true;break;end
 end
 if ~accepted,break;end
end
end
function [f,J]=forward(net,p,refs,d,dim)
f=zeros(d*numel(refs),1);J=zeros(numel(f),3);
for j=1:numel(refs)
 b=refs(j);delta=net.pb(:,b)-p;r=norm(delta);q=delta/r;ix=(j-1)*d+(1:d);axes=net.R(:,2:d+1,b);
 f(ix)=-axes'*q/2;J(ix,:)=(axes'*(eye(3)-q*q')/(2*r));
end
J=J(:,1:dim);
end
