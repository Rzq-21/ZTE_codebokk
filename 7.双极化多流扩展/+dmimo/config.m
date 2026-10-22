function c=config
c.version='DUAL_V1_20261022';c.qA=3;c.qP=4;c.XPR_dB=10;
c.noise_seed_offset=8000000;c.pilot_seed=202610220;
c.E=[0 0;1 -1;1 1]/sqrt(2);
rng(20261022,'twister');c.path_phase=2*pi*rand(2,2,6);
v=randn(1,10);v=v-mean(v);c.port_phase_rad=v/sqrt(mean(v.^2))*15*pi/180;
c.scenarios={'nominal','tilt20','phase15','combined'};
c.axes=unique([repmat(8,5,1),(20:5:40)';(5:10)',repmat(30,6,1)],'rows');
c.hardware='AP: 2*Nap, UE: 2*4 ports, same spatial coordinates';
c.model='controlled transverse Jones extension; not full TR38.901';
c.power='total 1; all realizable precoders per-AP 1/B';
c.capacity='true-channel waterfilling, rank <= N, total power only';
c.frequency_snapshots='strongest 4 observed energy columns of 16; no truth selection';
end
