function E=cpu_model(net,p,active,F,c)
% CPU input: decoded coordinate, public axes, geometry and known beam indices.
% Optional known calibration is an explicitly labelled oracle diagnostic.
allowed={'E','ue_tilt_deg','port_phase_rad'};
assert(isempty(setdiff(fieldnames(c),allowed)),'dmimo:CPUBoundary','Unexpected CPU model input');
tilt=0;phi=zeros(1,10);if isfield(c,'ue_tilt_deg'),tilt=c.ue_tilt_deg;end
if isfield(c,'port_phase_rad'),phi=c.port_phase_rad;end
B=numel(active);E=zeros(8,2*B);
for slot=1:B
 b=active(slot);delta=net.pb(:,b)-p;r=norm(delta);q=delta/r;ql=net.R(:,:,b)'*q;
 a=exp(-1j*net.k*r)/r*exp(1j*net.k*net.ue_offsets.'*q);
 gain=exp(-1j*net.k*ql.'*net.local)*F(:,slot);
 J=dmimo.jones(net,c,p,b,0,tilt);Q=kron(J,a*gain);Q(:,2)=Q(:,2)*exp(1j*phi(b));E(:,2*slot-1:2*slot)=Q;
end
end
