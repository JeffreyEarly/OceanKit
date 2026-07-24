classdef AlongTrackSimulatorBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    methods (Abstract)
        missions = missions;
        [lat,lon,time] = groundTrackForMissionWithName;
    end

    methods

        function alongtrack = projectedPointsForReferenceOrbit(self,options)
            arguments
                self AlongTrackSimulatorBase
                options.Lx
                options.Ly
                options.lat0
                options.lon0
            end
            optionsArgs = namedargs2cell(options);
            [x0, y0, minLat, minLon, maxLat, maxLon] = AlongTrackSimulatorBase.LatitudeLongitudeBoundsForTransverseMercatorBox(optionsArgs{:});
            [lat,lon,time] = self.groundTrackForMissionWithName("s6a");
            withinBox = lat >= minLat & lat <= maxLat & lon >= minLon & lon <= maxLon;
            lat(~withinBox) = [];
            lon(~withinBox) = [];
            time(~withinBox) = [];
            use options;
            [x,y] = AlongTrackSimulatorBase.LatitudeLongitudeToTransverseMercator(lat,lon,lon0=lon0);
            out_of_bounds = (x < x0 - Lx/2) | (x > x0 + Lx/2) | (y < y0 - Ly/2) | (y > y0 + Ly/2);
            alongtrack.x = x(~out_of_bounds);
            alongtrack.y = y(~out_of_bounds)-y0;
            alongtrack.time = time(~out_of_bounds);
            [alongtrack.time,I] = sort(alongtrack.time);
            alongtrack.x = alongtrack.x(I);
            alongtrack.y = alongtrack.y(I);
        end
    end

    methods (Static)
        [x,y] = LatitudeLongitudeToTransverseMercator(lat, lon, options)
        [lat,lon] = TransverseMercatorToLatitudeLongitude(x, y, options)
        [x0, y0, minLat, minLon, maxLat, maxLon] = LatitudeLongitudeBoundsForTransverseMercatorBox(options)

        y = MeridionalArcPROJ4(phi)
        phi = InverseMeridionalArcPROJ4(y)
    end

end