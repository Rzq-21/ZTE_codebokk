function design=ue_phase(net,prediction,vobs,Q)
% UE-only calculation, for the exact beams produced by the coordinate decoder.
design=prediction;B=size(design.F,2);a=sum(conj(design.F).*reshape(vobs,net.Nap,B),1).';
[design.c,design.phase_indices]=upa.phase_weights(a,Q);design.w=upa.assemble(design.F,design.c);design.phase_source='UE_estCSI';
design.AP_phase_bits=(B-1)*Q;design.feedback_bits=design.direction_bits+design.AP_phase_bits;
end
