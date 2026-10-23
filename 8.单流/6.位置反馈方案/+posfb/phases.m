function [cpu,ue]=phases(net,prediction,active,vobs,phase_bits)
% CPU uses only the public decoded prediction. vobs is used only by UE branch.
cpu=posfb.cpu_phase(net,prediction,active);ue=posfb.ue_phase(net,prediction,vobs,phase_bits);
assert(isequal(cpu.position,ue.position)&&isequal(cpu.beam_indices,ue.beam_indices));
end
