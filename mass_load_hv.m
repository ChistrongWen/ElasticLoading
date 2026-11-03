function [UEN] = mass_load_hv(Gfun,mass_LLZ,obs,rho)
% 由格林函数算质量负荷源在obs点产生的负荷位移 输出为mm
% Gfun Farrell, 1972 给出的Green函数表格，规则化后的。

R_e    = 6371e3;
obslon = obs(1);
obslat = obs(2);
lon    = mass_LLZ.lon(1,:); lon = lon(:); Nlon = length(lon);
lat    = mass_LLZ.lat(:,1); lat = lat(:);
h      = mass_LLZ.rg; % 1 means mass

Dlon = lon(2)-lon(1);
Dlat = lat(2)-lat(1);
latarea = calArea( lat, Dlat, Dlon )*R_e^2; %m^2
area = repmat(latarea,1,Nlon); %面积矩阵
[arclen,az] = distance(mass_LLZ.lat,mass_LLZ.lon,obslat,obslon); % arclen 角距离, az第二点相对于第一点的方位角

GfunU = interp1( Gfun.theta, Gfun.ur2, arclen ) ./ (1.0E12*R_e*deg2rad(arclen)); % 2-D m/kg
GfunH = interp1( Gfun.theta, Gfun.v2,  arclen ) ./ (1.0E12*R_e*deg2rad(arclen)); % 2-D m/kg
% 
% *-sin(az)东向位移
% *-cos(az)北向位移

% GfunU [m/kg]
% GfunH [m/kg]
% h [m]
% area [m^2]
% rho [kg/m^3]

Udis =            GfunU.*(h.*area*rho); %    Up displacement (m)
Edis = -sind(az).*GfunH.*(h.*area*rho); %  East displacement (m)
Ndis = -cosd(az).*GfunH.*(h.*area*rho); % North displacement (m)

Udis = Udis(:);
Edis = Edis(:);
Ndis = Ndis(:);

Udis(isnan(Udis))=0;
Edis(isnan(Edis))=0;
Ndis(isnan(Ndis))=0;

Udis = sum(Udis);
Edis = sum(Edis);
Ndis = sum(Ndis);

UEN=[Udis*1000,Edis*1000,Ndis*1000]; %mm

end