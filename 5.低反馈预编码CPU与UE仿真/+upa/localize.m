function [position,diagnostic]=localize(net,refs,frequency)
A=zeros(3);rhs=zeros(3,1);outside=0;
for j=1:numel(refs)
 y=-2*frequency(1,j);z=-2*frequency(2,j);outside=outside+(y*y+z*z>1);q=[-sqrt(max(0,1-y*y-z*z));y;z];q=q/norm(q);g=-net.R(:,:,refs(j))*q;P=eye(3)-g*g';A=A+P;rhs=rhs+P*net.pb(:,refs(j));
end
if rcond(A)<1e-14,error('upa:LocalizationSingular','Rank-deficient ray localization; do not drop sample');end
position=A\rhs;assert(all(isfinite(position)));diagnostic=struct('condition',cond(A),'outside_unit_disk_count',outside,'outside_user_region',any(position<[20;20;0])||any(position>[30;30;5]));
end
