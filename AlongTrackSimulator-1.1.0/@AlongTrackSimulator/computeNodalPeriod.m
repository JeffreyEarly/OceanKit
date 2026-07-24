function T_nodal = computeNodalPeriod(a, e, i)
% Computes the nodal precession period from orbital parameters
%
% - Topic: Utilities — Computing orbital parameters
% - Declaration: T_nodal = computeNodalPeriod(a, e, i)
% - Parameter a: semi-major axis [km]
% - Parameter e: eccentricity (0 for circular, >0 for elliptical)
% - Parameter i: inclination [degrees]
% - Returns T_nodal: nodal period [s]
arguments
    a
    e
    i
end
% Constants
mu = 398600.4418; % Earth's gravitational parameter (km^3/s^2)
J2 = 1.08262668e-3;      % Earth's J2 coefficient
RE = 6378.1363;        % Earth's equatorial radius (km)
omega = 0;

% Convert inclination to radians
i = deg2rad(i);

% Compute Mean Motion
n = sqrt(mu / a^3); % Mean motion (rad/s)

% Compute Semi-Latus Rectum
p = a * (1 - e^2);

% Compute Nodal Precession Rate
RAAN_dot = - (3/2) * J2 * (RE / p)^2 * n * cos(i); % RAAN precession rate (rad/s)

A = 3 * J2 * (4 - 5 * sin(i) * sin(i)) / ( 4 * (a / RE)^2 * sqrt(1-e^2) * (1 + e*cos(omega))^2 );
B = 3 * J2 * (1 + e*cos(omega))^3 / ( 2 * (a / RE)^2 * (1-e^2)^3 );
factor = 1 - A - B;

% Compute Nodal Period
T_nodal = 2 * pi * factor / n;
end
