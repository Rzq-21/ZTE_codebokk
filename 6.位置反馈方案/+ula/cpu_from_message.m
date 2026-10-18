function design=cpu_from_message(net,message,active_ids)
% CPU reconstructs a geometry-only channel and derives continuous co-phases.
prediction=ula.predict_directions(net,message,active_ids);
Hmodel=ula.reconstruct_channel(net,prediction.position,active_ids);
[~,~,v]=upa.principal(Hmodel);B=numel(active_ids);a=zeros(B,1);
for slot=1:B
 ix=(slot-1)*net.Nap+(1:net.Nap);a(slot)=prediction.F(:,slot)'*v(ix);
end
[prediction.c,~]=upa.phase_weights(a,Inf);prediction.w=upa.assemble(prediction.F,prediction.c);
prediction.phase_source='CPU_model';prediction.AP_phase_bits=0;prediction.feedback_bits=prediction.direction_bits;design=prediction;
end
