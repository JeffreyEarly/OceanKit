classdef AlongTrackSimulator < AlongTrackSimulatorBase
    % Simulate and manipulate satellite along-track ground tracks for altimetry missions.
    %
    % AlongTrackSimulator provides utilities to compute orbit ground tracks (optionally with
    % nodal precession), enforce exact repeat cycles for supported missions, and project
    % tracks into a local Cartesian box for along-track sampling and modeling workflows.
    %
    % - Topic: Initialization
    % - Topic: Constant orbital parameters
    % - Topic: Mission catalog
    % - Topic: Mission catalog — Fetching parameters
    % - Topic: Ground tracks
    % - Topic: Utilities
    % - Topic: Utilities — Working with alongtrack data
    % - Topic: Utilities — WaveVortexModel integration
    % - Topic: Utilities — Computing orbital parameters
    % - Topic: Utilities — Forcing exact repeat orbits
    % - Topic: Utilities — Groundtrack algorithms
    % - Topic: Utilities — Solving Kepler's equation
    %
    % - Declaration: classdef AlongTrackSimulator < AlongTrackSimulatorBase
    properties
        % Mission parameter dictionary keyed by mission abbreviation.
        %
        % A dictionary (string -> struct) describing orbit geometry and metadata
        % for supported altimetry missions. The constructor loads a default catalog and
        % may adjust some parameters (e.g., semi-major axis) to enforce exact repeat cycles.
        %
        % - Topic: Mission catalog
        missionParameters
    end

    properties (Constant)
        % Earth's standard gravitational parameter μ ($$km^3/s^2$$).
        %
        % Used for orbital period computations via Kepler's third law.
        %
        % - Topic: Constant orbital parameters
        mu = 398600.4418

        % Earth's second zonal harmonic coefficient J2.
        %
        % Used to approximate nodal precession and nodal period corrections.
        %
        % - Topic: Constant orbital parameters
        J2 = 1.08262668e-3

        % Earth's equatorial radius (km).
        %
        % Used for simple Earth geometry in ground-track calculations.
        %
        % - Topic: Constant orbital parameters
        RE = 6378.1363

        % Length of a sidereal day (seconds).
        %
        % Used when mapping inertial motion to Earth-fixed longitude.
        %
        % - Topic: Constant orbital parameters
        T_sidereal = 86164
    end

    methods

        function self = AlongTrackSimulator()
            % Create an AlongTrackSimulator instance.
            %
            % Loads the mission parameter catalog and (when applicable) adjusts mission
            % parameters so that repeat-cycle missions have an exact repeat by solving for
            % a semi-major axis consistent with the requested passes-per-cycle.
            %
            % - Topic: Initialization
            % - Declaration: self = AlongTrackSimulator()
            % - Returns self: AlongTrackSimulator instance
            self.missionParameters = AlongTrackSimulator.missionParametersCatalog();
            missions = self.missionParameters.keys;
            for i=1:numel(missions)
                % We force perfect repeat cycles by adjusting the semi-major axis.
                if isfinite(self.missionParameters(missions(i)).passes_per_cycle)
                    semi_major_axis = AlongTrackSimulator.semimajorAxisForExactRepeatForMission(missions(i));
                    % pctChange = abs(1-self.missionParameters(missions(i)).semi_major_axis/semi_major_axis)*100;
                    % disp("Changed semi-major axis from " + self.missionParameters(missions(i)).semi_major_axis + " to " + semi_major_axis + " (" + pctChange + "%) " + "for mission " + missions(i));
                    self.missionParameters(missions(i)).semi_major_axis = semi_major_axis;
                    T_orbit = self.orbitalPeriodForMissionWithName(missions(i));
                    self.missionParameters(missions(i)).repeat_cycle = T_orbit*self.missionParameters(missions(i)).passes_per_cycle/2/86400;
                end
            end
        end

        function [lat,lon,time] = groundTrackForMissionWithName(self,mission,options)
            % Compute the ground track for a mission over one or more orbits.
            %
            % Returns latitude/longitude as a function of time for the specified mission's
            % orbital elements. By default, the track starts at the mission start date and
            % spans options.N_orbits orbits; you may override the time vector.
            %
            % - Topic: Ground tracks
            % - Declaration: [lat,lon,time] = groundTrackForMissionWithName(mission,options)
            % - Parameter mission: text (mission key in missionParameters)
            % - Parameter options.time: numeric or datetime vector (optional) — sample times
            % - Parameter options.N_orbits: double (default 1) — number of orbits to simulate when options.time is not provided
            % - Returns lat: double column vector — latitude [deg]
            % - Returns lon: double column vector — longitude [deg]
            % - Returns time: same type as options.time or numeric seconds — time stamps for each point
            arguments
                self AlongTrackSimulator
                mission {mustBeText}
                options.time
                options.N_orbits = 1
            end
            p = self.missionParameters(mission);

            T_orbit = 2*pi*sqrt((p.semi_major_axis)^3/self.mu);  % Orbital period [s]
            if isfield(options,"time")
                time = options.time;
            else
                dt = 1;                    % Time step [s]
                time = 0:dt:options.N_orbits*T_orbit;           % Time vector
            end
            time_shifted = time-T_orbit/4;

            % RAAN = pi/2;              % Right Ascension of Ascending Node [rad], or longitude of the ascending node
            % argPerigee = -pi/2;       % Argument of perigee [rad]. -pi/2 starts the orbit on the ascending track
            % M0 = pi/2;                   % Initial mean anomaly [rad]
            RAAN = p.longitude_at_equator;    % Right Ascension of Ascending Node [rad], or longitude of the ascending node
            argPeriapsis = 0;         % Argument of periapsis [rad].
            M0 = 0;                   % Initial mean anomaly [rad]
            % [lat, lon] = AlongTrackSimulator.computeGroundTrack(p.semi_major_axis, p.eccentricity, p.inclination, RAAN, argPeriapsis, M0, time_shifted);
            [lat, lon] = AlongTrackSimulator.computeGroundTrackWithNodalPrecession(p.semi_major_axis, p.eccentricity, p.inclination, RAAN, argPeriapsis, M0, time_shifted);
        end

        function [lat,lon,time] = repeatGroundTrackForMissionWithName(self,mission)
            % Compute the full repeat-cycle ground track for a mission.
            %
            % For missions with a defined repeat cycle, returns the ground track over one
            % complete repeat period (repeat_cycle days). For missions
            % without a defined repeat cycle, this function will throw an
            % error.
            %
            % - Topic: Ground tracks
            % - Declaration: [lat,lon,time] = repeatGroundTrackForMissionWithName(mission)
            % - Parameter mission: text (mission key in missionParameters)
            % - Returns lat: double column vector — latitude [deg]
            % - Returns lon: double column vector — longitude [deg]
            % - Returns time: numeric column vector — elapsed seconds from start of repeat cycle
            missionParametersDict = self.missionParameters;
            p = missionParametersDict(mission);
            N_orbits = p.passes_per_cycle/2;
            if isinf(N_orbits)
                error("The mission " + mission + " does not have a repeat cycle.")
            end
            [lat,lon,time] = self.groundTrackForMissionWithName(mission,N_orbits=N_orbits);
        end

        function T = orbitalPeriodForMissionWithName(self,mission)
            % Return the orbital period for a mission.
            %
            % - Topic: Mission catalog — Fetching parameters
            % - Declaration: T = orbitalPeriodForMissionWithName(mission)
            % - Parameter mission: text (mission key in missionParameters)
            % - Returns T: double — orbital period [s]
            arguments
                self AlongTrackSimulator
                mission {mustBeText}
            end
            a = self.missionParameters(mission).semi_major_axis;
            T = 2*pi*sqrt(a^3/self.mu);  % Orbital period [s]
        end

        function T = nodalPeriodForMissionWithName(self,mission)
            % Return the nodal period for a mission (including J2 nodal precession).
            %
            % Uses a J2-based correction to compute the period in the Earth-fixed frame.
            %
            % - Topic: Mission catalog — Fetching parameters
            % - Declaration: T = nodalPeriodForMissionWithName(mission)
            % - Parameter mission: text (mission key in missionParameters)
            % - Returns T: double — nodal period [s]
            arguments
                self AlongTrackSimulator
                mission {mustBeText}
            end
            p = self.missionParameters(mission);
            T = self.computeNodalPeriod(p.semi_major_axis,p.eccentricity, p.inclination);  % Orbital period [s]
        end

        function T = repeatCycleForMissionWithName(self,mission)
            % Return the mission repeat cycle length in seconds.
            %
            % - Topic: Mission catalog — Fetching parameters
            % - Declaration: T = repeatCycleForMissionWithName(mission)
            % - Parameter mission: text (mission key in missionParameters)
            % - Returns T: double — repeat cycle duration [s]
            arguments
                self AlongTrackSimulator
                mission {mustBeText}
            end
            p = self.missionParameters(mission);
            T = p.repeat_cycle*86400;
        end

        function missions = missions(self)
            % List all mission keys in the catalog.
            %
            % - Topic: Mission catalog
            % - Declaration: missions = missions()
            % - Returns missions: string column vector — mission keys
            arguments
                self AlongTrackSimulator
            end
            missions = self.missionParameters.keys;
        end

        function missions = currentMissions(self)
            % List mission keys for currently operating satellites.
            %
            % - Topic: Mission catalog
            % - Declaration: missions = currentMissions()
            % - Returns missions: string column vector — mission keys
            arguments
                self AlongTrackSimulator
            end
            allMissions = self.missionParameters.keys;
            missions = string.empty(0,0);
            for i=1:numel(allMissions)
                if isinf(self.missionParameters(allMissions(i)).end_date)
                    missions(end+1) = allMissions(i);
                end
            end
            missions = reshape(missions,[],1);
        end

        function summarizeMissionWithName(self,missionNames)
            % Display a table summarizing one or more missions.
            %
            % Prints mission name, abbreviation, start/end dates, and repeat cycle.
            %
            % - Topic: Mission catalog
            % - Declaration: summarizeMissionWithName(missionNames)
            % - Parameter missionNames: string or string array — mission keys to summarize
            arguments
                self AlongTrackSimulator
                missionNames string
            end
            Name = strings(numel(missionNames),1);
            Abbreviation = strings(numel(missionNames),1);
            Start = strings(numel(missionNames),1);
            End = strings(numel(missionNames),1);
            RepeatCycle = strings(numel(missionNames),1);
            for i=1:numel(missionNames)
                Abbreviation(i) = missionNames(i);
                Name(i) = self.missionParameters(missionNames(i)).name;
                Start(i) = string(self.missionParameters(missionNames(i)).start_date);
                End(i) = string(self.missionParameters(missionNames(i)).end_date);
                RepeatCycle(i) = string(self.missionParameters(missionNames(i)).repeat_cycle);
            end
            T = table(Name,Abbreviation,Start,End,RepeatCycle);
            disp(T);
        end

        function outputGroups = addMissionsToOutputFile(self,outputFile,options)
            % Add along-track output groups for one or more missions.
            %
            % By default this method adds every mission returned by `currentMissions` and
            % samples sea-surface height. Supply any registered WaveVortexModel variables
            % with dimensions `{'x','y'}` through `fieldNames`.
            %
            % ```matlab
            % fields = ["ssh","ssu","ssv"];
            % outputGroups = ats.addMissionsToOutputFile(outputFile,fieldNames=fields);
            % ```
            %
            % - Topic: Utilities — WaveVortexModel integration
            % - Declaration: outputGroups = addMissionsToOutputFile(outputFile,options)
            % - Parameter outputFile: WVModelOutputFile scalar — output file receiving the mission groups
            % - Parameter options.missionNames: string array — mission keys; default is every current mission
            % - Parameter options.fieldNames: string array — registered horizontal fields to sample; default `"ssh"`
            % - Returns outputGroups: WVModelOutputGroupAlongTrack row vector — groups added to the output file
            arguments
                self AlongTrackSimulator
                outputFile WVModelOutputFile
                options.missionNames string
                options.fieldNames string = "ssh"
            end
            if isfield(options,"missionNames")
                missionNames = reshape(options.missionNames,1,[]);
            else
                missionNames = reshape(self.currentMissions,1,[]);
            end
            outputGroupCells = cell(1,length(missionNames));
            for iMission = 1:length(missionNames)
                outputGroup = WVModelOutputGroupAlongTrack(outputFile.model,missionNames(iMission),self,fieldNames=options.fieldNames);
                outputFile.addOutputGroup(outputGroup);
                outputGroupCells{iMission} = outputGroup;
            end
            if isempty(outputGroupCells)
                outputGroups = WVModelOutputGroupAlongTrack.empty(1,0);
            else
                outputGroups = [outputGroupCells{:}];
            end
        end

        function alongtrack = projectedPointsForMissionWithName(self,missionName,requiredOptions,options)
            % Project mission ground-track points into a local Transverse Mercator box.
            %
            % Computes the ground track for the requested mission within a bounding box
            % centered at (lon0,lat0) with size (Lx,Ly). Points are
            % returned in projected coordinates (using a transverse
            % Mercator projector). By default the origin is set to the
            % lower-left of the bounding box, but can be optionally set to
            % the center.
            % 
            % The function returns a struct with fields (x,y,t)
            %   - x: double column vector — projected x-coordinates
            %   - y: double column vector — projected y-coordinates
            %   - t: column vector — times corresponding to each point (sorted)
            %
            % - Topic: Ground tracks
            % - Declaration: alongtrack = projectedPointsForMissionWithName(missionName,requiredOptions,options)
            % - Parameter missionName: string — mission key
            % - Parameter requiredOptions.Lx: double — box width [m]
            % - Parameter requiredOptions.Ly: double — box height  [m]
            % - Parameter requiredOptions.lat0: double — box reference latitude [deg]
            % - Parameter requiredOptions.lon0: double — central meridian for projection [deg]
            % - Parameter options.time: numeric/datetime vector — times at which to sample the orbit
            % - Parameter options.origin: {'lower-left','center'} (default 'lower-left') — coordinate origin convention for returned (x,y)
            % - Returns: alongtrack struct with field (x,y,t)
            arguments
                self AlongTrackSimulator
                missionName string
                requiredOptions.Lx
                requiredOptions.Ly
                requiredOptions.lat0
                requiredOptions.lon0
                options.time
                options.origin {mustBeMember(options.origin,{'lower-left','center'})} = 'lower-left'
            end
            optionsArgs = namedargs2cell(requiredOptions);
            [x0, y0, minLat, minLon, maxLat, maxLon] = self.LatitudeLongitudeBoundsForTransverseMercatorBox(optionsArgs{:});
            [lat,lon,time] = self.groundTrackForMissionWithName(missionName,time=options.time);

            % 1) apply crude filter
            withinBox = lat >= minLat & lat <= maxLat & lon >= minLon & lon <= maxLon;
            lat(~withinBox) = [];
            lon(~withinBox) = [];
            time(~withinBox) = [];

            % 2) project points
            lon0 = requiredOptions.lon0;
            Lx = requiredOptions.Lx;
            Ly = requiredOptions.Ly;
            [x,y] = AlongTrackSimulatorBase.LatitudeLongitudeToTransverseMercator(lat,lon,lon0=lon0);

            % 3) apply more precise filter
            out_of_bounds = (x < x0 - Lx/2) | (x > x0 + Lx/2) | (y < y0 - Ly/2) | (y > y0 + Ly/2);
            alongtrack.x = x(~out_of_bounds);
            alongtrack.y = y(~out_of_bounds)-y0;
            alongtrack.t = reshape(time(~out_of_bounds),[],1);

            [alongtrack.t,I] = sort(alongtrack.t);

            if options.origin == "lower-left"
                alongtrack.x = alongtrack.x(I) + Lx/2;
                alongtrack.y = alongtrack.y(I) + Ly/2;
            else
                alongtrack.x = alongtrack.x(I);
                alongtrack.y = alongtrack.y(I);
            end
        end

        
        function alongtrack = projectedPointsForRepeatMissionWithName(self,missionName,requiredOptions,options)
            % Project repeat-cycle ground-track points into a local Transverse Mercator box.
            %
            % Similar to -projectedPointsForMissionWithName, but works only
            % for missions with repeat cycles and returns exactly one
            % repeat cycle through the box.
            %
            % Computes the ground track for the requested mission within a bounding box
            % centered at (lon0,lat0) with size (Lx,Ly). Points are
            % returned in projected coordinates (using a transverse
            % Mercator projector). By default the origin is set to the
            % lower-left of the bounding box, but can be optionally set to
            % the center.
            % 
            % The function returns a struct with fields (x,y,t)
            %   - x: double column vector — projected x-coordinates
            %   - y: double column vector — projected y-coordinates
            %   - t: column vector — times corresponding to each point (sorted)
            %   - repeatCycle: time of the repeat cycle [s]
            %
            % - Topic: Ground tracks
            % - Declaration: alongtrack = projectedPointsForMissionWithName(missionName,requiredOptions,options)
            % - Parameter missionName: string — mission key
            % - Parameter requiredOptions.Lx: double — box width [m]
            % - Parameter requiredOptions.Ly: double — box height  [m]
            % - Parameter requiredOptions.lat0: double — box reference latitude [deg]
            % - Parameter requiredOptions.lon0: double — central meridian for projection [deg]
            % - Parameter options.time: numeric/datetime vector — times at which to sample the orbit
            % - Parameter options.origin: {'lower-left','center'} (default 'lower-left') — coordinate origin convention for returned (x,y)
            % - Returns: alongtrack struct with field (x,y,t,repeatCycle)
            arguments
                self AlongTrackSimulator
                missionName string
                requiredOptions.Lx
                requiredOptions.Ly
                requiredOptions.lat0
                requiredOptions.lon0
                options.origin {mustBeMember(options.origin,{'lower-left','center'})} = 'lower-left'
            end
            optionsArgs = namedargs2cell(requiredOptions);
            [x0, y0, minLat, minLon, maxLat, maxLon] = self.LatitudeLongitudeBoundsForTransverseMercatorBox(optionsArgs{:});
            [lat,lon,time] = self.repeatGroundTrackForMissionWithName(missionName);

            % 1) apply crude filter
            withinBox = lat >= minLat & lat <= maxLat & lon >= minLon & lon <= maxLon;
            lat(~withinBox) = [];
            lon(~withinBox) = [];
            time(~withinBox) = [];

            % 2) project points
            lon0 = requiredOptions.lon0;
            Lx = requiredOptions.Lx;
            Ly = requiredOptions.Ly;
            [x,y] = AlongTrackSimulatorBase.LatitudeLongitudeToTransverseMercator(lat,lon,lon0=lon0);

            % 3) apply more precise filter
            out_of_bounds = (x < x0 - Lx/2) | (x > x0 + Lx/2) | (y < y0 - Ly/2) | (y > y0 + Ly/2);
            alongtrack.x = x(~out_of_bounds);
            alongtrack.y = y(~out_of_bounds)-y0;
            alongtrack.t = reshape(time(~out_of_bounds),[],1);

            [alongtrack.t,I] = sort(alongtrack.t);

            if options.origin == "lower-left"
                alongtrack.x = alongtrack.x(I) + Lx/2;
                alongtrack.y = alongtrack.y(I) + Ly/2;
            else
                alongtrack.x = alongtrack.x(I);
                alongtrack.y = alongtrack.y(I);
            end
            alongtrack.repeatCycle = self.repeatCycleForMissionWithName(missionName);
        end

    end

    methods (Static)

        missionData = missionParametersCatalog();

        [lat, lon] = computeGroundTrack(semi_major_axis, e, incl, RAAN_0, argPerigee, M0, t)

        [lat, lon] = computeGroundTrackWithNodalPrecession(semi_major_axis, e, incl, RAAN, omega, M0, t)

        [lat, lon] = computeGroundTrackWithNodalPrecessionSimple(semi_major_axis, e, incl, RAAN, omega, M0, t)

        [lat, lon] = computeGroundTrackCircularOrbit(semi_major_axis, e, incl, RAAN, omega, M0, t)

        T_nodal = computeNodalPeriod(a, e, i)

        function outputGroup = wvmOutputGroupForRepeatMissionWithName(model,missionName,options)
            % Build a WaveVortexModel output group for an along-track sampling pattern.
            %
            % This compatibility helper delegates to `WVModelOutputGroupAlongTrack`.
            %
            % - Topic: Utilities — WaveVortexModel integration
            % - Declaration: outputGroup = wvmOutputGroupForRepeatMissionWithName(model,missionName,options)
            % - Parameter model: WVModel
            % - Parameter missionName: string — mission key
            % - Parameter options.fieldNames: string array — registered horizontal fields to sample; default `"ssh"`
            % - Returns outputGroup: WVModelOutputGroupAlongTrack
            arguments
                model WVModel
                missionName string
                options.fieldNames string = "ssh"
            end
            ats = AlongTrackSimulator();
            outputGroup = WVModelOutputGroupAlongTrack(model,missionName,ats,fieldNames=options.fieldNames);
        end

        % function addObservingSystemToModelForRepeatMissionWithName(model, missionName)
        %     outputGroup = self.wvmOutputGroupForRepeatMissionWithName(model, missionName);
        % 
        % end

        function tracks = convertAlongTrackStructureToPass(alongtrack)
            % Convert a flat alongtrack struct into a cell array of passes.
            %
            % Splits an alongtrack struct (with fields x, y, t) whenever there is a
            % discontinuity in time (diff(t) > 1).
            %
            % - Topic: Utilities — Working with alongtrack data
            % - Declaration: tracks = convertAlongTrackStructureToPass(alongtrack)
            % - Parameter alongtrack: struct — must contain fields x, y, t
            % - Returns tracks: cell array — each cell contains a struct for one pass
            arguments
                alongtrack struct
            end
            trackIndices = find(diff(alongtrack.t)>1);
            trackIndices(end+1) = length(alongtrack.t);
            startIndex = 1;
            tracks = cell(length(trackIndices),1);
            for i=1:length(trackIndices)
                endIndex = trackIndices(i);
                tracks{i}.x = alongtrack.x(startIndex:endIndex);
                tracks{i}.y = alongtrack.y(startIndex:endIndex);
                tracks{i}.t = alongtrack.t(startIndex:endIndex);
                startIndex = endIndex+1;
            end
        end

        function inclination = inclinationForExactRepeatForMission(mission)
            % Compute inclination required for an exact repeat orbit for the mission.
            %
            % - Topic: Utilities — Forcing exact repeat orbits
            % - Declaration: inclination = inclinationForExactRepeatForMission(mission)
            % - Parameter mission: text — mission key
            % - Returns inclination: double — inclination (units per implementation)
            arguments
                mission {mustBeText}
            end
            mu = AlongTrackSimulator.mu;
            J2 = AlongTrackSimulator.J2;
            RE = AlongTrackSimulator.RE;
            T_sidereal = AlongTrackSimulator.T_sidereal;

            missionParametersDict = AlongTrackSimulator.missionParametersCatalog;
            p = missionParametersDict(mission);

            a = p.semi_major_axis; e = p.eccentricity; i = p.inclination;
            N_orbits = p.passes_per_cycle/2;
            T_orbit = 2*pi*sqrt(a^3 / mu);
            p = a * (1 - e^2);
            gamma = (3/2) * J2 * (RE / p)^2;

            M = round(N_orbits * (T_orbit/T_sidereal + gamma*cosd(i)));
            inclination = acosd((M/N_orbits - T_orbit/T_sidereal)/gamma);
        end

        function a = semimajorAxisForExactRepeatForMission(mission)
            % Solve for semi-major axis that yields an exact repeat cycle for the mission.
            %
            % - Topic: Utilities — Forcing exact repeat orbits
            % - Declaration: a = semimajorAxisForExactRepeatForMission(mission)
            % - Parameter mission: text — mission key
            % - Returns a: double — semi-major axis
            arguments
                mission {mustBeText}
            end
            mu = AlongTrackSimulator.mu;
            J2 = AlongTrackSimulator.J2;
            RE = AlongTrackSimulator.RE;
            T_sidereal = AlongTrackSimulator.T_sidereal;

            missionParametersDict = AlongTrackSimulator.missionParametersCatalog;
            p = missionParametersDict(mission);

            a = p.semi_major_axis; e = p.eccentricity; i = p.inclination;
            N_orbits = p.passes_per_cycle/2;
            T_orbit = 2*pi*sqrt(a^3 / mu);
            p = a * (1 - e^2);
            gamma = (3/2) * J2 * (RE / p)^2;
            M = round(N_orbits * (T_orbit/T_sidereal + gamma*cosd(i)));

            f = @(a) abs(((2*pi/T_sidereal) * sqrt(a^3/mu) + (3/2) * cosd(i) * J2 * (RE / (a * (1 - e^2)))^2)*N_orbits - M);

            a = fminsearch(f,a);
        end

        function J2 = j2ForExactRepeatForMission(mission)
            % Return the J2 value used in exact-repeat calculations.
            %
            % - Topic: Utilities — Forcing exact repeat orbits
            % - Declaration: J2 = j2ForExactRepeatForMission(mission)
            % - Parameter mission: text — mission key
            % - Returns J2: double — J2 coefficient
            arguments
                mission {mustBeText}
            end
            mu = AlongTrackSimulator.mu;
            J2 = AlongTrackSimulator.J2;
            RE = AlongTrackSimulator.RE;
            T_sidereal = AlongTrackSimulator.T_sidereal;

            missionParametersDict = AlongTrackSimulator.missionParametersCatalog;
            p = missionParametersDict(mission);

            a = p.semi_major_axis; e = p.eccentricity; i = p.inclination;
            N_orbits = p.passes_per_cycle/2;
            T_orbit = 2*pi*sqrt(a^3 / mu);
            p = a * (1 - e^2);
            gamma = (3/2) * J2 * (RE / p)^2;
            M = round(N_orbits * (T_orbit/T_sidereal + gamma*cosd(i)));

            f = @(J2) abs(((2*pi/T_sidereal) * sqrt(a^3/mu) + (3/2) * cosd(i) * (J2*1e-4) * (RE / (a * (1 - e^2)))^2)*N_orbits - M);

            J2 = fminsearch(f,J2*1e4)/1e4;
        end

        function e = eccentricityForExactRepeatForMission(mission)
            % Return the eccentricity used in exact-repeat calculations.
            %
            % - Topic: Utilities — Forcing exact repeat orbits
            % - Declaration: e = eccentricityForExactRepeatForMission(mission)
            % - Parameter mission: text — mission key
            % - Returns e: double — eccentricity
            arguments
                mission {mustBeText}
            end
            mu = AlongTrackSimulator.mu;
            J2 = AlongTrackSimulator.J2;
            RE = AlongTrackSimulator.RE;
            T_sidereal = AlongTrackSimulator.T_sidereal;

            missionParametersDict = AlongTrackSimulator.missionParametersCatalog;
            p = missionParametersDict(mission);

            a = p.semi_major_axis; e = p.eccentricity; i = p.inclination;
            N_orbits = p.passes_per_cycle/2;
            T_orbit = 2*pi*sqrt(a^3 / mu);
            p = a * (1 - e^2);
            gamma = (3/2) * J2 * (RE / p)^2;
            M = round(N_orbits * (T_orbit/T_sidereal + gamma*cosd(i)));

            f = @(e) abs(((2*pi/T_sidereal) * sqrt(a^3/mu) + (3/2) * cosd(i) * J2 * (RE / (a * (1 - (1e-5*e)^2)))^2)*N_orbits - M);

            e = fminsearch(f,e*1e5)/1e5;
        end

        function E = kepler1(M,e)
            % Solve Kepler's equation for eccentric anomaly E.
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler1(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]

            % Solve Kepler's Equation for Eccentric Anomaly using Newton-Raphson
            E = M; % Initial guess
            for j = 1:10 % Iterate to improve accuracy
                E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
            end
        end

        function E = kepler2(M,e)
            % Solve Kepler's equation for eccentric anomaly E.
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler2(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]

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
        end

        function E = kepler2vec(M, e)
            % Solve Kepler's equation for eccentric anomaly E (fastest version)
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler2vec(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]

            %KEPLER2  Vectorized solution of Kepler's equation: M = E - e*sin(E)
            %   E = kepler2(M, e) solves for the eccentric anomaly E given mean anomaly M
            %   and eccentricity e.  M and e may be arrays of the same size.

            % Convergence parameters
            tol = 1e-8;       % tolerance on deltaE
            maxIter = 100;    % maximum number of Newton-Raphson iterations

            % Initial guess (good for e < 1)
            E = M;
            % Track which elements have converged
            converged = false(size(M));

            for iter = 1:maxIter
                % Compute residuals and derivative
                f = E - e .* sin(E) - M;
                fprime = 1 - e .* cos(E);
                % Compute update
                deltaE = -f ./ fprime;
                % Update only unconverged elements
                E(~converged) = E(~converged) + deltaE(~converged);
                % Check convergence for each element
                newlyConverged = abs(deltaE) < tol;
                converged = converged | newlyConverged;
                % If all have converged, exit early
                if all(converged, 'all')
                    break;
                end
            end
        end


        function E = kepler3(M,e)
            % Solve Kepler's equation for eccentric anomaly E.
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler3(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]

            % Define Kepler's Equation as a function to minimize
            kepler_eq = @(E) abs(E - e * sin(E) - M); % Absolute error in Kepler's equation

            % Use fminsearch to find E that minimizes the function
            % optimset('TolFun', 1e-10, 'TolX', 1e-10,'MaxIter',10)
            E = fminsearch(kepler_eq, M,optimset('TolFun', 1e-10, 'TolX', 1e-10));
        end

        function E = kepler4(M,e)
            % Solve Kepler's equation for eccentric anomaly E.
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler4(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]
            E = M; % Initial guess
            tol = 1e-8;
            max_iter = 10;

            for j = 1:max_iter
                f = E - e * sin(E) - M;
                fp = 1 - e * cos(E);
                fpp = e * sin(E);

                dE = -f / (fp - 0.5 * f * fpp / fp); % Householder's update
                E = E + dE;

                if abs(dE) < tol
                    break;
                end
            end

        end


        function E = kepler5(M,e)
            % Solve Kepler's equation for eccentric anomaly E.
            %
            % Given mean anomaly M and eccentricity e, solves M = E - e*sin(E) using
            % a particular iterative scheme. Different kepler* variants implement different
            % initial guesses and/or update formulas; see code for details.
            %
            % - Topic: Utilities — Solving Kepler's equation
            % - Declaration: E = kepler5(M,e)
            % - Parameter M: double (scalar or vector) — mean anomaly [rad]
            % - Parameter e: double (scalar) — eccentricity
            % - Returns E: double (same size as M) — eccentric anomaly [rad]
            
            % Hybrid Newton-Raphson & Halley's Method
            if e < 0.8
                E = M + e * sin(M); % Good initial guess for low e
            else
                E = pi; % Better guess for high e
            end

            for j = 1:4  % Converges even faster than Newton
                f_E = E - e * sin(E) - M;
                f_E_deriv = 1 - e * cos(E);
                f_E_2nd_deriv = e * sin(E);

                % Halley's method update
                E = E - (2 * f_E * f_E_deriv) / (2 * f_E_deriv^2 - f_E * f_E_2nd_deriv);
            end
        end

    end

end
