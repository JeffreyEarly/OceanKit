function [lat, lon] = computeGroundTrackWithNodalPrecessionSimple(semi_major_axis, e, incl, RAAN_0, omega, M0, t)
% computes the ground track of a satellite given including the nodal precession (unoptimized version)
%
% - Topic: Utilities — Groundtrack Algorithms
% - Declaration: [lat, lon] = computeGroundTrackWithNodalPrecessionSimple(semi_major_axis, e, incl, RAAN_0, omega, M0, t)
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
    RAAN_0 
    omega 
    M0 
    t 
end

% Constants
mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)
J2 = 1.08263e-3;      % Earth's J2 coefficient
RE = 6378.137;        % Earth's equatorial radius (km)
omega_E = 7.2921159e-5; % Earth's rotation rate (rad/s)

% Convert inclination to radians
RAAN_0 = deg2rad(RAAN_0);
incl = deg2rad(incl);
omega = deg2rad(omega);
M0 = deg2rad(M0);

% Compute Mean Motion
n = sqrt(mu / semi_major_axis^3); % Mean motion (rad/s)

% Compute Nodal Precession Rate
p = semi_major_axis * (1 - e^2); % Semi-latus rectum
RAAN_dot = - (3/2) * J2 * (RE / p)^2 * n * cos(incl); % RAAN precession rate (rad/s)

% Initialize outputs
lat = zeros(size(t));
lon = zeros(size(t));

R1_i = [1, 0, 0; 0, cos(incl), -sin(incl); 0, sin(incl), cos(incl)]; % Rotate by inclination
R3_w = [cos(omega), -sin(omega), 0; sin(omega), cos(omega), 0; 0, 0, 1]; % Rotate by argument of perigee
R13_iw = R1_i * R3_w;

% Loop through each time step
for k = 1:length(t)
    % Compute updated RAAN due to nodal precession
    RAAN = RAAN_0 + RAAN_dot * t(k);

    % Compute Mean Anomaly at time t
    M = M0 + n * t(k);

    % kepler2 appears to be the fastest
    E = AlongTrackSimulator.kepler2(M,e);

    % Compute True Anomaly
    theta = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Compute Position in Orbital Plane
    r = p / (1 + e * cos(theta)); % Radius (km)
    x_orb = r * cos(theta);
    y_orb = r * sin(theta);

    % Rotation Matrices to transform to ECI
    R3_W = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1]; % Rotate by RAAN
    % R1_i = [1, 0, 0; 0, cos(incl), -sin(incl); 0, sin(incl), cos(incl)]; % Rotate by inclination
    % R3_w = [cos(omega), -sin(omega), 0; sin(omega), cos(omega), 0; 0, 0, 1]; % Rotate by argument of perigee

    % Compute ECI coordinates
    % r_eci = R3_W * R1_i * R3_w * [x_orb; y_orb; 0];
    r_eci = R3_W * R13_iw * [x_orb; y_orb; 0];

    % Convert to ECEF (Earth rotates at ω_E)
    theta_GMST = omega_E * t(k); % Earth's rotation angle

    R3_E = [cos(theta_GMST), sin(theta_GMST), 0; -sin(theta_GMST), cos(theta_GMST), 0; 0, 0, 1]; % Earth rotation matrix
    r_ecef = R3_E * r_eci;

    % Convert ECEF to Latitude and Longitude
    lon(k) = atan2(r_ecef(2), r_ecef(1));
    lat(k) = atan2(r_ecef(3), sqrt(r_ecef(1)^2 + r_ecef(2)^2));
end

% Convert output to degrees
lat = rad2deg(lat);
lon = rad2deg(lon);

% Ensure longitude is within -180 to 180 range
lon = mod(lon + 180, 360) - 180;
end
