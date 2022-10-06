function [sunel,TOA,LOD] = solgeo(lat,lon,time,lsm)
% solgeo drives sun elevation (sunel) [rad], top of atmosphere radiation
%(TOA) [W/m²] and length of day (LOD) [h] from given input. TOA is in UTC,
%not in local time
% Rainer Prinz, 19 Jan 2022
%   Input needed:
    % lat: latitude of location [dd.dddd]
    % lon: longitutde of location [dd.dddd]
    % time: time vector as datetime (must include minutes)
    % lsm: local standard meridian (°), e.g. Vienna (GMT+1) = 15°, Toronto
    % (GMT-5) = -75°

Io = 1362;  % solar constant, W/m²
doy = day(time,'dayofyear'); % day of year
mn = minute(time);
gamma = 2*pi*(doy-1)/365;      % Iqbal 1983 (1.2.2) day angle
Eo = 1.000110+0.034221.*cos(gamma)+0.001280.*sin(gamma)+0.000719...
    .*cos(2.*gamma)+0.000077.*sin(2.*gamma); % excentricity
delta = (0.006918-0.399912.*cos(gamma)+0.070257.*sin(gamma)-0.006758...
    .*cos(2.*gamma)+0.000907.*sin(2.*gamma)-0.002697.*cos(3.*gamma)...
    +0.00148.*sin(3.*gamma)); % declination
Et = ( 0.000075 + 0.001868.*cos(gamma) - 0.032077.*sin(gamma)...
     - 0.014615.*cos(2.*gamma) - 0.04089.*sin(2.*gamma) )...
         *229.18; % Equation of time
LATd = Et + 4*(lon-lsm); % Iqbal 1983 (1.4.2) difference to local apparent time
time_LT = time + minutes(LATd); % local apparent time
hro = hour(time_LT)+minute(time_LT)/60; % decimal hours of LAT [h]
omega = (hro-12).*deg2rad(15); % hour angle 0 = 12h LAT [h]
latrad = deg2rad(lat);
sunel = asin(sin(latrad).*sin(delta) + cos(latrad).*cos(delta).*cos(omega)); % sun elevation; complement angle to zenith angle
TOA = Io.*Eo.*sin(sunel); % Top of atmosphere radiation [rad]
omegas = acos(tan(latrad)*-1.*tan(delta));
omegas(omegas~=real(omegas)) = NaN;
LOD = 2/15.*rad2deg(omegas); % length of day [h]

end

