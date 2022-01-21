function [M,Lin,LinZ,Lout,sensible,latent,G]=ebalance(Tsfc,Lin_coarse,...
    Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
    rho_air,De_h,B1,B2,Kd,mode,tau,cc,opt_arg)

%incoming longwave
LinZ=lapseLongWave(Lin_coarse,Zdiff);
% adjust for surrounding terrain, assume snow covered emissivity, but doesn't really
% matter
LinT = Vf.* LinZ + (1-Vf).*emissivity.ems.*sig.*Tsfc.^4;

switch mode
    case 'normal'
        G=zeros(size(Lin_coarse),'single');
        es0=ice_vapor_pressure(Tsfc);
        xLs=ones(size(De_h))*xLs;
        %use latent heat of vaporization instead of sublimation 
        ind=es0<ea;
        xLs(ind)=2.260e6; %J/kg
        % neutral stability for starters, then fill in for unstable/stable
        stability = ones(size(Lin_coarse),'single');
        % unstable
        ind=Tsfc > T_fine;
        B3= 1 + B2(ind).*sqrt(Tsfc(ind)-T_fine(ind));
        stability(ind) = 1 + B1(ind).*(Tsfc(ind)-T_fine(ind))./B3;
        % stable
        ind=Tsfc < T_fine;
        B8=B1(ind)./2;
        stability(ind) = 1./((1+B8.*(T_fine(ind)-Tsfc(ind))).^2);

        latent=rho_air.*xLs.*De_h.*stability.*(0.622*((ea-es0)./pres_fine));
        %  adjust for vegetation emitted longwave
        Lin = ((tau.*LinT + (1-tau).*emissivity.emc.*sig.*T_fine.^4).*cc)...
            + (LinT.*(1-cc));
        Lout = -(emissivity.ems.*sig.*Tsfc.^4+(1-emissivity.ems).*Lin);
    otherwise
        d=opt_arg;
        % switched to daily mean, linear temp gradient assumed
        G=-(Kd.*(Tsfc-273.15))./d; 
        Lin=LinT; % i.e. no veg correction
        Lout = -(emissivity.emd.*sig.*Tsfc.^4+(1-emissivity.ems).*Lin);
        stability = ones(size(Lin_coarse),'single'); % Rounce and McKinney 2014
        %assume dry debris
        latent=zeros(size(Lin_coarse),'single');
end

sensible=rho_air.*Cp.*De_h.*stability.*(T_fine-Tsfc);
L=Lin+Lout;
T=sensible+latent;
% total ebalance
M=Sin.*(1-albedo)+L+T+G;
end