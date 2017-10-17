function Fin0 = partitionGlobalErbs1982(Gin0,T)
%Partition global into beam and diffuse using Erbs et al. (1982),Solar
%Energy.

%INPUT
% Gin0 - coarse resolution global incoming
% T - total transmittance
%OUTPUT
% Fin0 - coarse resolution diffuse

%Allocate diffuse 
Fin0=zeros(size(T),'single');

%Cloudy
cloudy = T <= 0.22;
Fin0(cloudy) =  Gin0(cloudy).*(1-0.09.*T(cloudy));

%Clear
clr = T > 0.80;
%Olyphant (1984) high altitude modification - 0.120 instead of 0.165
Fin0(clr) = Gin0(clr).*0.120.*T(clr);

%Partly cloudy
pc = ~cloudy & ~clr;
Fin0(pc) = Gin0(pc).*(0.9511 - 0.1604.*T(pc)+4.388.*T(pc).^2 ...
- 16.638.*T(pc).^3 + 12.366.*T(pc).^4);

% Check for negative Fin0
Fin0(Fin0<0 | ~isfinite(Fin0))=0;

% Check for infiniti values of Fin0
% Fin0(isinf(Fin0))=max(max(Fin0(~isinf(Fin0))));

% Check for Fin0 values greater than Gin0
% Fin0(Fin0>Gin0)=Gin0(Fin0>Gin0);

%Beam is remainder of global
% Bin0 = Gin0 - Fin0;