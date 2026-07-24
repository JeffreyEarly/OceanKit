classdef AlongTrackSimulatorEmpirical < AlongTrackSimulatorBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        path
    end

    methods
        function self = AlongTrackSimulatorEmpirical(options)
            arguments
                options.path = "./orbital_paths.nc"
            end
            self.path = options.path;
        end

        function [lat,lon,time] = groundTrackForMissionWithName(self,mission)
            arguments
                self AlongTrackSimulatorEmpirical
                mission {mustBeText}
            end
            lat = ncread("@AlongTrackSimulatorEmpirical/" + self.path, mission + "/latitude");
            lon = ncread("@AlongTrackSimulatorEmpirical/" + self.path, mission + "/longitude");
            time = days(ncread("@AlongTrackSimulatorEmpirical/" + self.path, mission + "/time")) + datetime(1950,01,01);
        end

        function missions = missions(self)
            info = ncinfo("@AlongTrackSimulatorEmpirical/" + self.path);
            missions = reshape(string({info.Groups.Name}),[],1);
        end
        
        % function alongtrack = projectedPointsForReferenceOrbit(self,options)
        %     arguments
        %         self AlongTrackSimulatorEmpirical
        %         options.Lx
        %         options.Ly
        %         options.lat0
        %         options.lon0
        %     end
        %     optionsArgs = namedargs2cell(options);
        %     [x0, y0, minLat, minLon, maxLat, maxLon] = AlongTrackSimulatorEmpirical.LatitudeLongitudeBoundsForTransverseMercatorBox(optionsArgs{:});
        %     [lat,lon,time] = self.pathForMission("s6a");
        %     withinBox = lat >= minLat & lat <= maxLat & lon >= minLon & lon <= maxLon;
        %     lat(~withinBox) = [];
        %     lon(~withinBox) = [];
        %     time(~withinBox) = [];
        %     use options;
        %     [x,y] = AlongTrackSimulatorEmpirical.LatitudeLongitudeToTransverseMercator(lat,lon,lon0=lon0);
        %     out_of_bounds = (x < x0 - Lx/2) | (x > x0 + Lx/2) | (y < y0 - Ly/2) | (y > y0 + Ly/2);
        %     alongtrack.x = x(~out_of_bounds);
        %     alongtrack.y = y(~out_of_bounds)-y0;
        %     alongtrack.time = time(~out_of_bounds);
        %     [alongtrack.time,I] = sort(alongtrack.time);
        %     alongtrack.x = alongtrack.x(I);
        %     alongtrack.y = alongtrack.y(I);
        % end
    end

end