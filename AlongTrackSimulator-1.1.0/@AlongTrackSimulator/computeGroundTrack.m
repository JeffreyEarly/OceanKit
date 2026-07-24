function [lat, lon] = computeGroundTrack(semi_major_axis, e, incl, RAAN_0, argPerigee, M0, t)
% computes the ground track of a satellite given including the nodal precession (unoptimized version)
%
% - Topic: Utilities — Groundtrack Algorithms
% - Declaration: [lat, lon] = computeGroundTrack(semi_major_axis, e, incl, RAAN, omega, M0, t)
% - Parameter semi_major_axis: Semi-major axis [km]
% - Parameter e: Eccentricity (0 for circular, >0 for elliptical)
% - Parameter incl: Inclination [degrees]
% - Parameter RAAN_0: Right Ascension of Ascending Node [rad]
% - Parameter argPerigee: Argument of perigee [rad]
% - Parameter M0: Initial mean anomaly at t = 0 [rad]
% - Parameter t: Time vector [s]
% - Returns lat: Latitude (deg) vector corresponding to times in t
% - Returns lon: Longitude (deg) vector corresponding to times in t
arguments
    semi_major_axis 
    e 
    incl 
    RAAN_0 
    argPerigee 
    M0 
    t 
end

% Define constants
mu = 398600.4418;         % Earth's gravitational parameter [km^3/s^2]
Re = 6378.137;                % Earth's radius [km]
omega_e = 7.2921159e-5;   % Earth's rotation rate [rad/s]

% Convert inclination to radians
incl = deg2rad(incl);

J2 = 1.08263e-3;      % Earth's J2 coefficient
p = semi_major_axis * (1 - e^2); % Semi-latus rectum
n = sqrt(mu / semi_major_axis^3); % Mean motion (rad/s)
RAAN_dot = - (3/2) * J2 * (Re / p)^2 * n * cos(incl); % RAAN precession rate (rad/s)

N = length(t);
lat = zeros(1, N);
lon = zeros(1, N);

a = semi_major_axis;
for i = 1:N
    % Compute mean anomaly at time t(i)
    M = M0 + n * t(i);
    
    RAAN = RAAN_0 + 0*RAAN_dot*t(i);

    % Solve Kepler's Equation: M = E - e*sin(E)
    % Using Newton-Raphson method for iterative solution
    E = M;            % Initial guess
    tol = 1e-8;       % Convergence tolerance
    maxIter = 100;    % Maximum iterations
    for iter = 1:maxIter
        f = E - e*sin(E) - M;
        fprime = 1 - e*cos(E);
        deltaE = -f / fprime;
        E = E + deltaE;
        if abs(deltaE) < tol
            break;
        end
    end
    
    % Compute true anomaly from eccentric anomaly
    theta = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    
    % Compute distance from Earth's center at this point in the orbit
    r = a * (1 - e*cos(E));
    
    % Position in the perifocal coordinate system (orbital plane)
    r_pf = [r*cos(theta); r*sin(theta); 0];
    
    % Transformation matrix from perifocal to ECI coordinates
    R_peri2eci = [ cos(RAAN)*cos(argPerigee)-sin(RAAN)*sin(argPerigee)*cos(incl), -cos(RAAN)*sin(argPerigee)-sin(RAAN)*cos(argPerigee)*cos(incl),  sin(RAAN)*sin(incl);
                   sin(RAAN)*cos(argPerigee)+cos(RAAN)*sin(argPerigee)*cos(incl), -sin(RAAN)*sin(argPerigee)+cos(RAAN)*cos(argPerigee)*cos(incl), -cos(RAAN)*sin(incl);
                   sin(argPerigee)*sin(incl),                              cos(argPerigee)*sin(incl),                               cos(incl) ];
               
    % Satellite position in ECI coordinates
    r_ECI = R_peri2eci * r_pf;
    
    % Transform from ECI to ECEF by accounting for Earth's rotation
    theta_e = omega_e * t(i);
    R_ecef = [ cos(theta_e)  sin(theta_e) 0;
              -sin(theta_e)  cos(theta_e) 0;
                   0             0         1];
    r_ECEF = R_ecef * r_ECI;
    
    % Convert ECEF coordinates to geodetic latitude and longitude
    x = r_ECEF(1);
    y = r_ECEF(2);
    z = r_ECEF(3);
    % Using a spherical Earth approximation
    lat(i) = rad2deg(asin(z / norm(r_ECEF)));
    lon(i) = rad2deg(atan2(y, x));
end
end
