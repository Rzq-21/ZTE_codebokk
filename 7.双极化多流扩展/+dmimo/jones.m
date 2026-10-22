function J=jones(net,c,p,b,ell,tilt)
% Public axes and path geometry; ell=0 denotes LOS. No CSI input.
t=tilt*pi/180;Ry=[cos(t) 0 sin(t);0 1 0;-sin(t) 0 cos(t)];Er=Ry*c.E;Et=net.R(:,:,b)*c.E;
if ell==0
 q=(p-net.pb(:,b));q=q/norm(q);J=Er'*(eye(3)-q*q')*Et;
else
 scatter=[10 40 10 40 25 25;10 10 40 40 5 45;zeros(1,6)];
 qt=scatter(:,ell)-net.pb(:,b);qr=p-scatter(:,ell);
 ph=c.path_phase(:,:,ell);M=exp(1j*ph).*[1 10^(-c.XPR_dB/20);10^(-c.XPR_dB/20) 1];
 J=Er'*transverse(qr)*M*transverse(qt)'*Et;
end
assert(norm(J,'fro')>1e-12);J=J*sqrt(2)/norm(J,'fro');
end
function T=transverse(q)
q=q/norm(q);a=[0;0;1];if abs(q'*a)>.9,a=[0;1;0];end
v=cross(q,a);v=v/norm(v);T=[v cross(q,v)];
end
