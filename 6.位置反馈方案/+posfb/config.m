function cfg=config(net)
% Frozen public coordinate quantizer. No Monte Carlo observations enter here.
cfg.id='POS_V1_PHASE_AP1_B8';cfg.array=net.array;cfg.allocation='phase';
if strcmp(net.array,'ULA64')
 cfg.dim=2;cfg.lo=[19;19];cfg.hi=[31;31];cfg.main_bits=24;cfg.budgets=[16 20 24 27 30 36 42 45];
 [x,y]=ndgrid(linspace(20,30,7));grid=[x(:).';y(:).';zeros(1,numel(x))];
else
 cfg.dim=3;cfg.lo=[19;19;-1];cfg.hi=[31;31;6];cfg.main_bits=36;cfg.budgets=[24 30 36 42 48 54 60 70 84];
 [x,y,z]=ndgrid(linspace(20,30,7),linspace(20,30,7),linspace(0,5,5));grid=[x(:).';y(:).';z(:).'];
end
ids=net.active_order(1:8);cfg.sensitivity_weights=zeros(cfg.dim,1);
for i=1:size(grid,2)
 g=(grid(:,i)-net.pb(:,ids))./vecnorm(grid(:,i)-net.pb(:,ids));D=g(:,2:end)-g(:,1);
 cfg.sensitivity_weights=cfg.sensitivity_weights+mean(D(1:cfg.dim,:).^2,2);
end
cfg.sensitivity_weights=cfg.sensitivity_weights/size(grid,2);cfg.grid_points=size(grid,2);
cfg.design_AP_ids=ids;cfg.domain_id='ROI_PLUS_1M_FIXED';cfg.minimum_bits_per_axis=4;
cfg.main_allocation=posfb.allocate(cfg,cfg.main_bits);assert(all(cfg.sensitivity_weights>0));
end
