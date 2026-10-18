function design=cpu_phase(net,prediction,active)
% CPU entry has no measured CSI or environment truth argument.
Hmodel=upa.reconstruct_channel(net,prediction.position,active);[~,~,vmodel]=upa.principal(Hmodel);
design=prediction;B=size(design.F,2);a=sum(conj(design.F).*reshape(vmodel,net.Nap,B),1).';
[design.c,design.phase_indices]=upa.phase_weights(a,Inf);design.w=upa.assemble(design.F,design.c);design.phase_source='CPU_model';
design.AP_phase_bits=0;design.feedback_bits=design.direction_bits;
end
