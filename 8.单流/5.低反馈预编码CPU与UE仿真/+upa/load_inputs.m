function [net,sim,bank]=load_inputs
% Split public CPU geometry from simulation-only truth/sample inputs.
root=fileparts(fileparts(mfilename('fullpath')));base=fullfile(root,'shared');addpath(fullfile(root,'lib'));
allnet=jsondecode(fileread(fullfile(base,'network_geometry.json')));sim=jsondecode(fileread(fullfile(base,'simulation_settings.json')));m=allnet.models.UPA32x32;
assert(strcmp(m.geometry_id,'UPA32_G50_C4_DFT_FINAL'),'Wrong deployment');
net=struct('array','UPA32x32','geometry_id',m.geometry_id,'Nside',32,'Nap',1024,'Nu',4,'M',128,'index_bits',14,'lambda',3e8/7e9,'active_order',m.active_order(:).','reference3',m.reference_AP_sets.x3(:).','reference5',m.reference_AP_sets.x5(:).');
% jsondecode changes numeric field names "3" and "5" to x3 and x5.
net.k=2*pi/net.lambda;net.d=net.lambda/2;net.ue_offsets=[(-1.5:1.5)*net.d;zeros(2,4)];[ny,nz]=ndgrid(0:31,0:31);net.ny=ny(:);net.nz=nz(:);net.local=[zeros(1,1024);(ny(:).'-15.5)*net.d;(nz(:).'-15.5)*net.d];
net.pb=zeros(3,10);net.R=zeros(3,3,10);net.ap_elements=zeros(3,1024,10);
for b=1:10
 assert(m.APs(b).physical_AP_id==b);net.pb(:,b)=m.APs(b).position_m;net.R(:,:,b)=m.APs(b).R_local_to_global;assert(norm(net.R(:,:,b)'*net.R(:,:,b)-eye(3),'fro')<1e-12);net.ap_elements(:,:,b)=net.pb(:,b)+net.R(:,:,b)*net.local;
end
d=vecnorm(reshape(net.pb,3,10,1)-reshape(net.pb,3,1,10));d=reshape(d,10,10);d(1:11:end)=Inf;assert(min(d(:))>=4-1e-8);
net.groups=zeros(1024,16);column=0;
for oz=0:3,for oy=0:3,column=column+1;[my,mz]=ndgrid(oy+4*(0:31),oz+4*(0:31));net.groups(:,column)=my(:)+128*mz(:)+1;end,end
positions=readtable(fullfile(base,sim.UE_positions_file),'VariableNamingRule','preserve');coeff=readtable(fullfile(base,sim.scatter_coefficients_file),'VariableNamingRule','preserve');seeds=readtable(fullfile(base,sim.noise_seeds_file),'VariableNamingRule','preserve');
assert(isequal(positions.sample_id,(1:400).'));bank.positions=[positions.x_m positions.y_m positions.UPA_z_m].';bank.zeta=zeros(400,6);bank.noise_seed=zeros(400,10);
for s=1:400
 rows=coeff(coeff.sample_id==s,:);assert(isequal(rows.scatterer_id,(1:6).'));bank.zeta(s,:)=complex(rows.real,rows.imag).';
 for b=1:10,row=seeds(strcmp(seeds.array,'UPA32x32')&seeds.sample_id==s&seeds.physical_AP_id==b,:);assert(height(row)==1);bank.noise_seed(s,b)=row.matlab_twister_seed;end
end
bank.shared_dir=base;bank.manifest=readtable(fullfile(base,'主实验运行清单.csv'),'VariableNamingRule','preserve');bank.budget=readtable(fullfile(base,'反馈开销_确定值.csv'),'VariableNamingRule','preserve');
sim.protocol_version='v2.0';sim.rho_data=10^(sim.reference_data_SNR_dB/10)/sim.reference_power;
end
