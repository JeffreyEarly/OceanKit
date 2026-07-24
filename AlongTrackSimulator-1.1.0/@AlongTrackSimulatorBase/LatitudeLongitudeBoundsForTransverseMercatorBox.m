function [x0, y0, minLat, minLon, maxLat, maxLon] = LatitudeLongitudeBoundsForTransverseMercatorBox(options)
arguments
    options.lat0 {mustBeNonempty}
    options.lon0 {mustBeNonempty}
    options.Lx {mustBeNonempty}
    options.Ly {mustBeNonempty}
end
lat0 = options.lat0;
lon0 = options.lon0;
Lx = options.Lx;
Ly = options.Ly;

[x0, y0] = AlongTrackSimulator.LatitudeLongitudeToTransverseMercator(lat0, lon0, lon0=lon0);
x = zeros(6,1);
y = zeros(6,1);

x(1) = x0 - Lx / 2;
y(1) = y0 - Ly / 2;

x(2) = x0 - Lx / 2;
y(2) = y0 + Ly / 2;

x(3) = x0;
y(3) = y0 + Ly / 2;

x(4) = x0;
y(4) = y0 - Ly / 2;

x(5) = x0 + Lx / 2;
y(5) = y0 - Ly / 2;

x(6) = x0 + Lx / 2;
y(6) = y0 + Ly / 2;

[lats, lons] = AlongTrackSimulator.TransverseMercatorToLatitudeLongitude(x, y, lon0=lon0);
minLat = min(lats);
maxLat = max(lats);
minLon = min(lons);
maxLon = max(lons);

end