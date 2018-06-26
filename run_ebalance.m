function [M,Tsfc,Lin,LinZ,Lout,sensible,latent,G,windspd,Td,ea,opt_out] = ...
    run_ebalance(todays_dateval,Lin_coarse,q_coarse,Zdiff,T_coarse,T_fine,...
    p_coarse,p_fine,albedo,Sin,topo,windS,fast_flag,mode,opt_in)
% calls solve ebalance in a parfor loop
% todays_dateval - MATLAB datenum, scalar integer
% Lin_coarse - coarse incoming longwave, W/m^2
% q_coarse - coarse specific humidity, kg/kg
% Zdiff - fine elevation difference from coarse elevation data
% corresponding to Lin_coarse (e.g. CERES Zdiff)
% T_coarse - coarse air temps, K
% T_fine - downscaled air temperatures, K
% p_coarse - coarse pressure, kPa
% p_fine - fine pressure, kPa
% albedo - snow or debris albedo, 0-1
% Sin - dowscaled global shortwave, W/m^2
% topo - structure with fine scale topographic information from GetTopography
% windS - wind structure containing wind flag and either:
% U,V (windflag true) or wind speed scalar (windflag false), m/s
% fast flag - true - only solve for M; false - solve for all outputs; only set for 'normal'
% mode, either normal - snow/ice melt;
% debris - ice melt under debris;
% or debris depth - solve for debris depth
% then follow with req'd inputs:
% for 'normal';
%   FOREST - forest structure
% for 'debris'
%   % d - debris cover depth, m
% for 'debris_depth'
%   Tsfc - surface temp, K
% output:
% M(Tsfc) - sum of melt energy, + indicates energy to melt, - indicates energy to
% freeze, W/m^2
% Tsfc - snow or debris surface temperature, K
% Lin - incoming longwave, W/m^2
% LinZ - incoming longwave, only corrected for elev., W/m^2
% Lout - outgoing longwave, W/m^2
% sensible - sensible heat flux, W/m^2
% latent - latent heat flux, W/m^2
% G - conduction, W/m^2
% Td - dewpoint temp, K
% ea - vapor pressure air, kPa
% windspd - wind speed, m/s
% optional output:
% opt_out - debris depth, m, only for 'debris_depth' mode; NaN otherwise

cc=zeros(size(Lin_coarse));
switch mode
    case 'normal'
        FOREST=opt_in;
        opt_input=FOREST.tau;
        cc=FOREST.cc;
        Tsfc=NaN(size(cc));
    case 'debris'
        d=opt_in;
        %debris depth as the optional input
        opt_input=d;
        Tsfc=NaN(size(opt_input));
    case 'debris depth'
        Tsfc=opt_in;
        %Tsfc as the optional input
        opt_input=Tsfc;
end

Vf=topo.view;
%if NLDAS U and V winds are present (windflag false)
if ~windS.windflag
    [u_fine,v_fine]=topo_winds(topo,FOREST,windS.U,windS.V,todays_dateval);
    windS.windspd=sqrt(u_fine.^2+v_fine.^2);
end
windspd=windS.windspd;

if fast_flag %call with matrices and don't solve for Tsfc
    [M,Tsfc,Lin,LinZ,Lout,sensible,latent,G,Td,ea]=...
        solve_ebalance(Lin_coarse,q_coarse,Zdiff,...
        T_coarse,T_fine,Vf,albedo,Sin,windspd,...
        p_coarse,p_fine,fast_flag,mode,opt_input,cc);
    opt_out=NaN(size(M));
else %call in a loop for each pixel and solve for Tsfc
    %reshape inputs
    opt_input=opt_input(:);
    cc=cc(:);
    sz=size(Lin_coarse);
    Lin_coarse=Lin_coarse(:);
    q_coarse=q_coarse(:);
    Zdiff=Zdiff(:);
    T_coarse=T_coarse(:);
    T_fine=T_fine(:);
    p_coarse=p_coarse(:);
    p_fine=p_fine(:);
    albedo=albedo(:);
    
    Vf=Vf(:);
    Sin=Sin(:);
    windspd=windspd(:);
    
    M=NaN(size(Lin_coarse));
    Lin=NaN(size(Lin_coarse));
    LinZ=NaN(size(Lin_coarse));
    Lout=NaN(size(Lin_coarse));
    sensible=NaN(size(Lin_coarse));
    latent=NaN(size(Lin_coarse));
    G=NaN(size(Lin_coarse));
    Td=NaN(size(Lin_coarse));
    ea=NaN(size(Lin_coarse));
    opt_output=NaN(size(Lin_coarse));
    
    for i=1:length(M)
        if ~isnan(albedo(i)) %skip areas w/o snow or debris
            [M(i),Tsfc(i),Lin(i),LinZ(i),Lout(i),sensible(i),latent(i),...
                G(i),Td(i),ea(i),opt_output(i)]=...
                solve_ebalance(Lin_coarse(i),q_coarse(i),Zdiff(i),...
                T_coarse(i),T_fine(i),Vf(i),albedo(i),Sin(i),windspd(i),...
                p_coarse(i),p_fine(i),fast_flag,mode,opt_input(i),cc(i));
        end
    end
    
    M=reshape(M,sz);
    Tsfc=reshape(Tsfc,sz);
    Lin=reshape(Lin,sz);
    LinZ=reshape(LinZ,sz);
    Lout=reshape(Lout,sz);
    sensible=reshape(sensible,sz);
    latent=reshape(latent,sz);
    windspd=reshape(windspd,sz);
    G=reshape(G,sz);
    Td=reshape(Td,sz);
    ea=reshape(ea,sz);
    if strcmp(mode,'debris depth')
        opt_out=reshape(opt_output,sz);
    else
        opt_out=NaN(sz);
    end
end
