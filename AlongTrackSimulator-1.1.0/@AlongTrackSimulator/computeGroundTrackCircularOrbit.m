function [lat, lon] = computeGroundTrackCircularOrbit(semi_major_axis, e, incl, RAAN, omega, M0, t)
% computes the ground track assuming a circular orbit
%
% - Topic: Utilities — Groundtrack Algorithms
% - Declaration: [lat, lon] = computeGroundTrackCircularOrbit(semi_major_axis, e, incl, RAAN_0, omega, M0, t)
% - Parameter semi_major_axis: semi-major axis [km]
% - Parameter e: eccentricity (0 for circular, >0 for elliptical)
% - Parameter incl: inclination [degrees]
% - Parameter RAAN_0: right ascension of ascending node [rad]
% - Parameter omega: argument of perigee [rad]
% - Parameter M0: initial mean anomaly at t = 0 [rad]
% - Parameter t: time vector [s]
% - Returns lat: latitude (deg) vector corresponding to times in t
% - Returns lon: longitude (deg) vector corresponding to times in t
arguments
    semi_major_axis 
    e 
    incl 
    RAAN 
    omega 
    M0 
    t 
end

% Define constants
mu = 398600.4418;         % Earth's gravitational parameter [km^3/s^2]
Re = 6378;                % Earth's radius [km]
omega_e = 7.2921159e-5;   % Earth's rotation rate [rad/s]

N = length(t);
lat = zeros(1, N);
lon = zeros(1, N);


r = semi_major_axis;
T_orbit = 2*pi*sqrt(r^3/mu); % Orbital period [s]
for i = 1:length(t)
    % Compute the orbital angle (true anomaly) at time t(i)
    theta = M0 + (2*pi/T_orbit) * t(i);
    
    % Satellite position in the orbital plane (ECI frame before inclination)
    r_orbital = [r*cos(theta); r*sin(theta); 0];  % [km]
    
    % Rotate to account for the orbital inclination (rotation about the x-axis)
    R_inc = [1      0           0;
             0  cos(incl) -sin(incl);
             0  sin(incl)  cos(incl)];
    r_ECI = R_inc * r_orbital;
    
    % Account for Earth rotation: Transform from ECI to ECEF frame.
    % The transformation rotates about the z-axis by the Earth rotation angle.
    theta_e = omega_e * t(i);      % Earth's rotation angle [rad]
    R_ecef = [cos(theta_e)  sin(theta_e) 0;
             -sin(theta_e)  cos(theta_e) 0;
                  0             0       1];
    r_ECEF = R_ecef * r_ECI;
    
    % Convert ECEF coordinates to geodetic latitude and longitude
    % (Assuming a spherical Earth)
    x = r_ECEF(1);
    y = r_ECEF(2);
    z = r_ECEF(3);
    lat(i) = rad2deg(asin(z / norm(r_ECEF)));  % Latitude [deg]
    lon(i) = rad2deg(atan2(y, x));              % Longitude [deg]
end

end
