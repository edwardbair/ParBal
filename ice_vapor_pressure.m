function es0=ice_vapor_pressure(Tsfc)
%vapor pressure over ice, Pa (Buck 1981)
Tf=273.15;
A=6.1115*100;
B=22.452;
C=272.55;
%ice surface vap. pressure, Pa
es0=A.*exp((B.*(Tsfc-Tf))./(C+(Tsfc-Tf)));