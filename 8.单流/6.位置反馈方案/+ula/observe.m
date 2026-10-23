function [Htrue,Hobs,Z,noise_variance]=observe(sample,channel,ids,snr)
if strcmp(channel,'LOS'),blocks=sample.LOS;elseif strcmp(channel,'LOS_NLOS_K10'),blocks=sample.MIX;else,error('ula:Channel','Unknown channel');end
Htrue=reshape(blocks(:,:,ids),4,[]);Z=reshape(sample.Z(:,:,ids),4,[]);noise_variance=mean(abs(Htrue(:)).^2)/10^(snr/10);Hobs=Htrue+sqrt(noise_variance)*Z;
end
