function design=ue_phase_from_message(net,message,active_ids,Hobs,phase_bits)
% UE derives quantized relative co-phases from noisy measured CSI.
design=ula.predict_directions(net,message,active_ids);[~,~,vobs]=upa.principal(Hobs);B=numel(active_ids);a=zeros(B,1);
for slot=1:B
 ix=(slot-1)*net.Nap+(1:net.Nap);a(slot)=design.F(:,slot)'*vobs(ix);
end
[design.c,design.phase_indices]=upa.phase_weights(a,phase_bits);design.w=upa.assemble(design.F,design.c);
design.phase_source='UE_estCSI';design.AP_phase_bits=(B-1)*phase_bits;design.feedback_bits=design.direction_bits+design.AP_phase_bits;
assert(design.phase_indices(1)==0&&numel(design.phase_indices)==B);
end
