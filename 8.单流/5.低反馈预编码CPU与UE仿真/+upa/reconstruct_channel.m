function H=reconstruct_channel(net,position,active_ids)
% CPU model: per-AP plane wave plus each AP's center-distance phase.
H=zeros(4,numel(active_ids)*net.Nap);
for slot=1:numel(active_ids)
 b=active_ids(slot);delta=net.pb(:,b)-position;r=norm(delta);assert(r>0);q=delta/r;ql=net.R(:,:,b)'*q;ix=(slot-1)*net.Nap+(1:net.Nap);
 H(:,ix)=exp(-1j*net.k*r)/r*exp(1j*net.k*net.ue_offsets.'*q)*exp(-1j*net.k*ql.'*net.local);
end
end
