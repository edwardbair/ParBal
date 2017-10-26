function [M,Tsfc,Lin,LinZ,Lout,sensible,latent,G,opt_out]=...
    solve_ebalance(Lin_coarse,q_coarse,Zdiff,T_coarse,T_fine,Vf,...
    albedo,Sin,windspd,pres_coarse,pres_fine,fast_flag,mode,varargin)

%solves energy balance for surface temperature for scalar (point) values
%input
% Lin_coarse - coarse incoming longwave, W/m^2
% q_coarse - coarse specific humidity, kg/kg
% T_coarse - coarse temp, K
% Zdiff - elevation diff between coarse and fine scale dem from Lin_coarse
% T_fine - downscaled air temperatures, K
% Vf - view factor, 0-1
% albedo - debris albedo, 0-1
% Sin - dowscaled global shortwave, W/m^2
% windspd -  wind speed, m/s
% pres_coarse - coarse pressure, kPa
% pres_fine - fine scale pressure, kPa
% fast flag - true - only solve for M; false - solve for all outputs; only set for 'normal'
% mode, either normal - snow/ice melt;
% debris - ice melt under debris;
% or debris depth - solve for debris depth
% then follow with req'd inputs:
% for 'normal';
%   tau - canopy transmissivity, 0-1
%   cc - canopy cover fraction, 0-1
% for 'debris'
%   % d - debris cover depth, m
% for 'debris depth'
%   Tsfc - surface temp, K
% NOTE in the mex version two inputs must always be given after mode since
% varargin is not implemented. Thus for debris of debris depth, provide a
% zero
% output:
% M- sum of melt energy, should be around 0
% Tsfc - debris surface temperature, K (solved for w/ energy balance or
% same as supplied)
% Lin - incoming longwave, W/m^2
% LinZ - incoming longwave, not corrected for veg or surrounding terrain,
% W/m^2
% Lout - outgoing longwave, W/m^2
% sensible, sensible heat flux, W/m^2
% latent, latent heat flux, W/m^2
% G, conductive heat flux, W/m^2
% optional, if mode is 'debris depth'
% opt_out - d, debris depth, m, NaN otherwise

% Constants
sig=single(5.6704E-8);% Stefan Boltzmann
emissivity.emd=single(0.94);% Emissivity of debris
emissivity.ems=single(0.99);% Emissivity of snow
emissivity.emc=single(0.98);% Emissivity of canopy

Cp=single(1010); % specific heat air, J/(kg deg)
Cs=single(2.09e3); %specific heat of ice, J deg C^-1
xLs=single(2.834e6); %latent heat of sublimation, J/(kg deg)
pres_0=single(101.325); % sea level pressure kPa.
rho_air=single(1.29).*(pres_fine./pres_0);%alt adj. density air,kg/m^3
%convert kPa to Pa
pres_fine=pres_fine*1000;
ht_wind_obs=single(2); %m
xkappa=single(0.41);% von Karmans constant
gravity=single(9.8); %m/s
z_0d=single(0.05); %debris roughness length, m, Lejeune and other, 2013
z_0s=single(0.0005); %snow roughness length, m
Kd=single(1.0); %W/(m deg K), debris avg from Schauwecker and other, 2015
%air vapor pressure
[~,ea]=downscaleDewpoint(pres_coarse,q_coarse,T_coarse,T_fine);
%convert kPa to Pa
ea=ea*1000;

Tf=273.15*ones(size(Lin_coarse),'single');%K, freezing temp

M=NaN(size(Lin_coarse),'single');
Lin=NaN(size(Lin_coarse),'single');
LinZ=NaN(size(Lin_coarse),'single');
Lout=NaN(size(Lin_coarse),'single');
sensible=NaN(size(Lin_coarse),'single');
latent=NaN(size(Lin_coarse),'single');
G=NaN(size(Lin_coarse),'single');

normalflag=false;
ddflag=false;

%debris case
z_0=z_0d;
Tsfc=NaN(size(Lin_coarse),'single');
tau=zeros(size(Lin_coarse),'single');
cc=zeros(size(Lin_coarse),'single');
ebalance_opt_arg=NaN(size(Lin_coarse),'single');
switch mode
    case 'normal'
        tau=varargin{1};
        cc=varargin{2};
        z_0=z_0s;
        normalflag=true;
    case 'debris'
        d=varargin{1};
        ebalance_opt_arg=d;
    case 'debris depth'
        Tsfc=varargin{1};
        ebalance_opt_arg=Tsfc;
        ddflag=true;
    otherwise
        error('mode must be normal, debris, or debris depth');
end

%compute turbulent coef
De_h=(xkappa.^2.*windspd)./((log(ht_wind_obs./z_0)).^2);
%additional stability function constants
C1=5.3.*9.4.*(xkappa./(log(ht_wind_obs./z_0))).^2.*...
    sqrt(ht_wind_obs./z_0);
C2=gravity.*ht_wind_obs./(T_fine.*windspd.^2);
B1=9.4.*C2;
B2=C1.*sqrt(C2);
opt_out=NaN(size(Lin_coarse),'single');
% if ~isnan(albedo)
    if ddflag % && ~isnan(Tsfc(:)) % Tsfc given
        [d,LinZ,Lout,sensible,G]=solve_debris_depth(Lin_coarse,Zdiff,T_fine,Vf,...
            albedo,Sin,Tsfc,sig,emissivity.ems,rho_air,Cp,De_h,Kd);
        opt_out=d;
    elseif ~ddflag %not solving for debris depth
        if ~fast_flag % solve w/ correct outputs
            x0=single(T_fine);
%             y0=ebalance(x0,Lin_coarse,...
%                 Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
%                 rho_air,De_h,B1,B2,Kd,mode,tau,cc,ebalance_opt_arg);
%             if ~isnan(x0) && ~isnan(y0) && isfinite(x0) && isfinite(y0)
                %solve for Tsfc
                y = newton_raphson_ebalance(x0,single(1e-5),single(20),Lin_coarse,...
                    Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
                    rho_air,De_h,B1,B2,Kd,mode,tau,cc,ebalance_opt_arg);
                %use that Tsfc to calculate ebalance
                Tsfc = y;
                if normalflag && y > Tf %don't adjust Tsfc for debris
                    Tsfc=Tf; %snow temp can't be greater than freezing
                end
                %note Tf as input since M is + for y > Tsfc and - for y < Tsfc
                [M,Lin,LinZ,Lout,sensible,latent,G]=ebalance(Tf,Lin_coarse,...
                    Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
                    rho_air,De_h,B1,B2,Kd,mode,tau,cc,ebalance_opt_arg);   
%             end
        elseif fast_flag
            %just solve using Tf - only the residual M is correct
            Tsfc=Tf;
            [M,Lin,LinZ,Lout,sensible,latent,G]=ebalance(Tsfc,Lin_coarse,...
                Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
                rho_air,De_h,B1,B2,Kd,mode,tau,cc,ebalance_opt_arg);  
        end
    end
end
% end