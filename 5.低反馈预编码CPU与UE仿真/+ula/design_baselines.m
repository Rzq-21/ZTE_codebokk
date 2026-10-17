function designs=design_baselines(net,Hobs,Htrue,active_ids,phase_bits)
% Common baselines: SVD, Type-II (two DFT beams/AP), and Type-I (one beam/AP).
B=numel(active_ids);[~,~,vtrue]=upa.principal(Htrue);[~,~,vobs]=upa.principal(Hobs);
F1=zeros(net.Nap,B);F2=zeros(net.Nap,B);m1=zeros(1,B);m2=zeros(2,B);amp_indices=zeros(1,B);internal_phase_indices=zeros(1,B);
for slot=1:B
 ix=(slot-1)*net.Nap+(1:net.Nap);[m1(slot),g]=ula.spectrum(net,vobs(ix));F1(:,slot)=ula.beams(net,m1(slot));
 gp=g(net.groups);[energy,order]=sort(abs(gp).^2,1,'descend');[~,group]=max(sum(energy(1:2,:),1));chosen=net.groups(order(1:2,group),group);values=g(chosen);m2(:,slot)=chosen-1;
 ratio=min(1,abs(values(2))/max(abs(values(1)),realmin));amp_indices(slot)=round(7*ratio);delta=angle(values(2))-angle(values(1));internal_phase_indices(slot)=mod(round(delta*4/(2*pi)),4);
 coefficient=amp_indices(slot)/7*exp(1j*2*pi*internal_phase_indices(slot)/4);f=ula.beams(net,m2(:,slot))*[1;coefficient];F2(:,slot)=f/norm(f);
end
template=struct('method_id','','phase_source','','F',[],'c',[],'w',[],'beam_indices',[],'position',nan(3,1),'localization',struct('condition',NaN,'outside_unit_disk_count',NaN,'outside_user_region',NaN),'direction_message',struct,'phase_indices',[],'internal_amp_indices',[],'internal_phase_indices',[],'direction_bits',NaN,'internal_bits',NaN,'AP_phase_bits',NaN,'feedback_bits',NaN,'reference_ap_ids',[]);
ideal=template;ideal.method_id='SVD';ideal.phase_source='ideal_truth';ideal.w=vtrue;
one=finish(template,'TypeI',F1,m1,vobs,net,phase_bits,B);one.direction_bits=9*B;one.internal_bits=0;one.feedback_bits=one.direction_bits+one.AP_phase_bits;
two=finish(template,'TypeII',F2,m2,vobs,net,phase_bits,B);two.direction_bits=18*B;two.internal_bits=5*B;two.feedback_bits=two.direction_bits+two.internal_bits+two.AP_phase_bits;two.internal_amp_indices=amp_indices;two.internal_phase_indices=internal_phase_indices;
designs={ideal,two,one};
end
function out=finish(template,name,F,indices,v,net,Q,B)
out=template;out.method_id=name;out.phase_source='UE_estCSI';out.F=F;out.beam_indices=indices;out.reference_ap_ids=[];a=zeros(B,1);
for slot=1:B,ix=(slot-1)*net.Nap+(1:net.Nap);a(slot)=F(:,slot)'*v(ix);end
[out.c,out.phase_indices]=upa.phase_weights(a,Q);out.w=upa.assemble(F,out.c);out.AP_phase_bits=(B-1)*Q;
end
