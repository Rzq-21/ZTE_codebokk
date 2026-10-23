function design=demo_cpu_message(method,R,channel)
% One complete UE -> integer message -> CPU -> evaluator example.
if nargin<1,method='DFT';end;if nargin<2,R=3;end;if nargin<3,channel='LOS';end
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'lib'));[net,sim,bank]=upa.load_inputs;
ids=net.active_order(1:8);sample=upa.sample_channels(net,sim,bank,1);[Htrue,Hobs]=upa.observe(sample,channel,ids,30);
if R==3,refs=net.reference3;elseif R==5,refs=net.reference5;else,error('R must be 3 or 5');end
messages=upa.encode_directions(net,Hobs,ids,refs,14);message=messages.(method);
% Only these two public objects cross the CPU boundary; the active AP IDs are public.
design=upa.cpu_from_message(net,message,ids);
[truth.u,truth.s2]=upa.principal(Htrue);metric=upa.evaluate(net,sim,Htrue,truth,design,sample.position);
disp(message);fprintf('Estimated position: [%.9f %.9f %.9f] m\n',design.position);fprintf('Feedback %d bits; sample G=%.6f; throughput=%.4f Mbps\n',design.feedback_bits,metric.power_gain,metric.throughput);
disp(table(ids(:),design.beam_indices(1,:).',design.beam_indices(2,:).',angle(design.c)*180/pi,'VariableNames',{'physical_AP_id','my','mz','CPU_phase_deg'}));
end
