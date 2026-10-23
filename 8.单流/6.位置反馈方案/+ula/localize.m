function [position,diagnostic]=localize(net,refs,frequency)
A=zeros(2);rhs=zeros(2,1);outside=0;
for j=1:numel(refs)
 qy=-2*frequency(j);outside=outside+(abs(qy)>1);qx=-sqrt(max(0,1-qy*qy));q=[qx;qy;0];q=q/norm(q);g=-net.R(:,:,refs(j))*q;g=g(1:2);g=g/norm(g);P=eye(2)-g*g';A=A+P;rhs=rhs+P*net.pb(1:2,refs(j));
end
if rcond(A)<1e-14,error('ula:LocalizationSingular','Rank-deficient 2-D ray localization; do not drop sample');end
xy=A\rhs;position=[xy;0];assert(all(isfinite(position)));diagnostic=struct('condition',cond(A),'outside_unit_disk_count',outside,'outside_user_region',any(xy<[20;20])||any(xy>[30;30]));
end
