function [Td,ea]=downscaleDewpoint(pres_coarse,q_coarse,T_coarse,T_fine)
% Dewpoint in degrees K
% in: 
% pres - coarse air pressure, kPa
% q - coarse specific humidity, kg/kg
% T_coarse, coarse air temp, K
% T_fine, downscaled air temp, K
% out: 
% Td, dewpoint in deg K
% ea, vapor pressure of air, kPa

%fix zero and negative values
thresh=7e-5;
q_coarse( q_coarse < thresh) = thresh;
%vapor pressure, Pa, Peixoto and Oort (1996) Eq 3, p 3445
ea_coarse=q_coarse.*(pres_coarse*1000)/0.622;
%compute rh
%assume rh constant between coarse and fine pixels
rh=ea_coarse./(SaturationVaporPressure(T_coarse,'water')*1000);
%recalc if any pix below freezing
if any(T_coarse(:)<=0)
    idx=T_coarse <= 0;
    rh(idx)=ea_coarse(idx)./...
        SaturationVaporPressure(T_coarse(idx),'ice');
end

rh(rh < 0.01) = 0.01;
rh(rh > 1) = 1;

T_fineC=T_fine-273.15;

b=zeros(size(T_fineC));
c=zeros(size(T_fineC));
b(:,:)=22.587;
c(:,:)=273.86;
%recalc for below freezing
if any(T_fineC(:)<=0)
    b(T_fineC<=0)=17.625;
    c(T_fineC<=0)=243.04;
end
% Rayleigh et al (2013) doi:10.1002/2013WR013958
TdC=c.*(log(rh)+b.*T_fineC./(c+T_fineC))./...
(b-log(rh)-b.*T_fineC./(c+T_fineC));
%convert to K
Td=TdC+273.15;
%use Td in sat. formula to obtain actual vapor pressure
ea=SaturationVaporPressure(Td,'water');
%recalc if below freezing
if any(T_fineC(:)<=0);
    idx=T_fineC<=0;
    ea(idx)=SaturationVaporPressure(Td(idx),'ice');
end

% a=0.6112;% kPa
% b=17.67;
% c=243.5;% deg C
% 
% Td=(c.*log(ea./a))./(b-log(ea./a));
% Td=Td+273.15; % convert to K