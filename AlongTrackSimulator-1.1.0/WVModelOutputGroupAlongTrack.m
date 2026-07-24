classdef WVModelOutputGroupAlongTrack < WVModelOutputGroup
    % Represent WaveVortexModel output for satellite along-track sampling.
    %
    % WVModelOutputGroupAlongTrack manages the time sampling and NetCDF output for a
    % single satellite altimetry mission using an AlongTrackSimulator. The group
    % precomputes mission pass-overs through the model domain and, when the model
    % reaches a pass-over time, writes the full along-track sample sequence for that
    % pass-over into the corresponding NetCDF group.
    %
    % Typical usage:
    % - Create a WVModel NetCDF output file and initialize an AlongTrackSimulator.
    % - Add the requested missions and registered horizontal fields to the output file.
    %
    % The following code adds output groups for all current satellites
    % ```matlab
    % outputFile = model.createNetCDFFileForModelOutput('ModelOutput.nc',outputInterval=86400);
    % ats = AlongTrackSimulator();
    % ats.addMissionsToOutputFile(outputFile,fieldNames=["ssh","ssu","ssv"]);
    % ```
    %
    % Major responsibilities:
    % - Store mission metadata and repeat-cycle information.
    % - Determine pass-over output times for a model integration window.
    % - Write the complete along-track time series for each pass-over into NetCDF.
    % - Provide class annotation metadata for property introspection.
    %
    % - Topic: Initialization
    % - Topic: Mission metadata
    % - Topic: Output scheduling
    % - Topic: NetCDF output
    % - Topic: Class annotations
    %
    % - Declaration: classdef WVModelOutputGroupAlongTrack < WVModelOutputGroup
    properties
        % Mission identifier used to configure the along-track sampling.
        %
        % Key used with the AlongTrackSimulator mission catalog (e.g., "s6a").
        %
        % - Topic: Mission metadata
        missionName

        % AlongTrackSimulator used to compute and project mission ground tracks.
        %
        % The simulator defines the orbit geometry and provides methods that project
        % latitude/longitude into the model's local Cartesian domain.
        %
        % - Topic: Mission metadata
        ats

        % Mission repeat cycle (s).
        %
        % For repeating missions this is the repeat-cycle duration in seconds. For
        % non-repeating missions this is Inf.
        %
        % - Topic: Mission metadata
        repeatCycle

        % Pass-over tracks through the model domain.
        %
        % Cell array of structs with fields x, y, and t, containing the projected
        % along-track coordinates (m) and model times (s) for each pass-over.
        %
        % - Topic: Output scheduling
        tracks

        % Model time of first sample for each pass-over (s).
        %
        % Column vector of the first time stamp in each element of tracks. This is
        % used to determine when the model should stop to write an along-track pass-over.
        %
        % - Topic: Output scheduling
        firstPassoverTime
    end

    properties (SetAccess=private)
        % Registered horizontal WaveVortexModel fields sampled along track.
        %
        % Each field must be registered with the transform and have dimensions
        % `{'x','y'}`. The default is `"ssh"`.
        %
        % - Topic: Mission metadata
        fieldNames (1,:) string
    end

    methods
        function self = WVModelOutputGroupAlongTrack(model,missionName,ats,options)
            % Create an along-track output group for a satellite mission.
            %
            % Initializes the output group and precomputes projected along-track pass-overs through
            % the model domain for repeat-cycle missions. The resulting tracks and first-passover
            % times are used to schedule output and to write pass-over samples into NetCDF. Each
            % requested field must be a registered two-dimensional WaveVortexModel variable with
            % dimensions `{'x','y'}`.
            %
            % - Topic: Initialization
            % - Declaration: self = WVModelOutputGroupAlongTrack(model,missionName,ats,options)
            % - Parameter model: WVModel scalar — parent model instance providing domain geometry and output file context
            % - Parameter missionName: text scalar — mission key used by AlongTrackSimulator
            % - Parameter ats: AlongTrackSimulator scalar — simulator used to compute and project tracks into the model domain
            % - Parameter options.fieldNames: string array — registered horizontal fields to sample; default `"ssh"`
            % - Returns self: WVModelOutputGroupAlongTrack instance
            arguments
                model WVModel
                missionName {mustBeText}
                ats AlongTrackSimulator
                options.fieldNames string = "ssh"
            end
            fieldNames = reshape(options.fieldNames,1,[]);
            if isempty(fieldNames)
                error("WVModelOutputGroupAlongTrack:EmptyFieldNames","fieldNames must contain at least one registered horizontal WaveVortexModel variable.");
            end
            if numel(unique(fieldNames)) ~= numel(fieldNames)
                error("WVModelOutputGroupAlongTrack:DuplicateFieldNames","fieldNames must not contain duplicate variable names.");
            end
            for iField = 1:length(fieldNames)
                if ~model.wvt.hasVariableWithName(char(fieldNames(iField)))
                    error("WVModelOutputGroupAlongTrack:UnknownField","The WaveVortexModel transform does not have a registered variable named '%s'.",fieldNames(iField));
                end
                annotation = model.wvt.propertyAnnotationWithName(fieldNames(iField));
                if ~isequal(annotation.dimensions,{'x','y'})
                    error("WVModelOutputGroupAlongTrack:InvalidFieldDimensions","The variable '%s' has dimensions {%s}; along-track fields must have dimensions {'x','y'}.",fieldNames(iField),strjoin(string(annotation.dimensions),", "));
                end
            end
            self@WVModelOutputGroup(model,name=missionName);
            self.missionName = missionName;
            self.ats = ats;
            self.fieldNames = fieldNames;
            self.repeatCycle = ats.repeatCycleForMissionWithName(missionName);

            if ~isinf(self.repeatCycle)
                wvt = self.model.wvt;
                alongtrack = ats.projectedPointsForRepeatMissionWithName(missionName,Lx=wvt.Lx,Ly=wvt.Ly,lat0=wvt.latitude,lon0=0);
                self.tracks = self.convertTrackVectorToPassoverCellArray(alongtrack);

                % only stop the model at the first t in the passover.
                self.firstPassoverTime = zeros(length(self.tracks),1);
                for iPassover = 1:length(self.tracks)
                    self.firstPassoverTime(iPassover) = self.tracks{iPassover}.t(1);
                end
            end
            self.addObservingSystem(WVAlongTrackObservingSystem(model,self));
        end

        function aString = description(self)
            % Describe the sampling pattern represented by this output group.
            %
            % Returns a concise, human-readable description string including repeat-cycle information
            % when the mission is configured as a repeat mission.
            %
            % - Topic: Mission metadata
            % - Declaration: aString = description()
            % - Returns aString: string scalar — description of mission sampling pattern
            arguments
                self WVModelOutputGroupAlongTrack
            end
            
            if ~isinf(self.repeatCycle)
                aString = "Sampling pattern for the " + self.missionName + " mission with a repeat cycle of " + string(self.repeatCycle) + "s.";
            else
                aString = "Sampling pattern for the " + self.missionName + " mission.";
            end
        end

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % Return output times within an integration window.
            %
            % For repeat-cycle missions, returns one time per repeat cycle corresponding to the first
            % sample of the pass-over through the model domain. For non-repeating missions, computes
            % pass-overs over [initialTime, finalTime] and returns the first sample time for each.
            %
            % - Topic: Output scheduling
            % - Declaration: t = outputTimesForIntegrationPeriod(initialTime,finalTime)
            % - Parameter initialTime: double scalar — integration window start time (s)
            % - Parameter finalTime: double scalar — integration window end time (s)
            % - Returns t: double column vector — pass-over trigger times within the integration window (s)
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end

            if ~isinf(self.repeatCycle)
                iCycle = floor(initialTime/self.repeatCycle);
                fCycle = floor(finalTime/self.repeatCycle);
                t = reshape(self.firstPassoverTime + (iCycle:fCycle)*self.repeatCycle,[],1);

                t(t<initialTime) = [];
                t(t>finalTime) = [];
                t = setdiff(t,self.timeOfLastIncrementWrittenToGroup);
            else
                wvt = self.model.wvt;
                % increment by 1s, typical along track sampling frequency
                t_full = ceil(initialTime):1:floor(finalTime);
                alongtrack = self.ats.projectedPointsForMissionWithName(self.missionName,Lx=wvt.Lx,Ly=wvt.Ly,lat0=wvt.latitude,lon0=0,time=t_full);
                self.tracks = self.convertTrackVectorToPassoverCellArray(alongtrack);

                self.firstPassoverTime = zeros(length(self.tracks),1);
                for iPassover = 1:length(self.tracks)
                    self.firstPassoverTime(iPassover) = self.tracks{iPassover}.t(1);
                end
                t = self.firstPassoverTime;
            end
        end

        function writeTimeStepToNetCDFFile(self,ncfile,t)
            % Write a model time step to the NetCDF output group.
            %
            % Overrides the WVModelOutputGroup behavior. When the model reaches a scheduled pass-over
            % time, this method writes the complete along-track time series for that pass-over into
            % the NetCDF group. The group time dimension and bookkeeping counters are advanced by the
            % number of samples written, while the model stop time remains the triggering time.
            %
            % - Topic: NetCDF output
            % - Declaration: writeTimeStepToNetCDFFile(t)
            % - Parameter ncfile: NetCDFFile — netcdf file being written to
            % - Parameter t: double scalar — model time at which output is triggered (s)
            arguments
                self WVModelOutputGroup
                ncfile NetCDFFile
                t double
            end
            if ~self.didInitializeStorage
                self.initializeOutputGroup(ncfile);
            end
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToGroup )
                if ~isinf(self.repeatCycle)
                    iPassover = find( abs(self.firstPassoverTime - mod(t,self.repeatCycle)) < 1, 1, 'first');
                else
                    iPassover = find( abs(self.firstPassoverTime - t) < 1, 1, 'first');
                end
                firstOutputIndex = self.incrementsWrittenToGroup + 1;
                lastOutputIndex = firstOutputIndex + length(self.tracks{iPassover}.t) - 1;
                outputIndices = firstOutputIndex:lastOutputIndex;
                outputTimes = t + self.tracks{iPassover}.t - self.tracks{iPassover}.t(1);

                self.group.variableWithName('t').setValueAlongDimensionAtIndex(outputTimes,'t',outputIndices);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).writeTimeStepToFile(self.group,outputIndices);
                end

                self.incrementsWrittenToGroup = outputIndices(end);
                self.timeOfLastIncrementWrittenToGroup = t(1);
            end
        end

    end

    methods (Static)

        function vars = classRequiredPropertyNames()
            % Return names of required properties for this class.
            %
            % Provides the list of properties that must be present for a valid output group instance
            % in workflows that perform validation or serialization.
            %
            % - Topic: Class annotations
            % - Declaration: vars = classRequiredPropertyNames()
            % - Returns vars: cell row vector of char — required property names
            vars = {'missionName','name','repeatCycle'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Define class property annotations for introspection and UI layers.
            %
            % Returns an array of CAPropertyAnnotation objects describing key properties, their units,
            % and semantic meaning for downstream tooling.
            %
            % - Topic: Class annotations
            % - Declaration: propertyAnnotations = classDefinedPropertyAnnotations()
            % - Returns propertyAnnotations: CAPropertyAnnotation vector — property annotation objects
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('missionName','name the mission');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of output group');
            propertyAnnotations(end+1) = CAPropertyAnnotation('fieldNames','registered horizontal fields sampled along track');
            propertyAnnotations(end+1) = CANumericProperty('repeatCycle', {}, 's','orbital repeat cycle (Inf indicates that the orbit is non-repeat)');
            % propertyAnnotations(end+1) = CANumericProperty('firstPassoverTime', {}, 's','model time of first passover into the model domain');
        end

        function tracks = convertTrackVectorToPassoverCellArray(alongtrack)
            % Convert a track time series into a cell array of pass-overs.
            %
            % Splits a continuous along-track struct (with fields t, x, y) into individual pass-overs
            % using gaps in time. Each output cell contains a struct with fields x, y, and t.
            %
            % - Topic: Output scheduling
            % - Declaration: tracks = convertTrackVectorToPassoverCellArray(alongtrack)
            % - Parameter alongtrack: struct — input track with fields:
            %   - t: double vector — model time stamps (s)
            %   - x: double vector — projected x coordinate (m)
            %   - y: double vector — projected y coordinate (m)
            % - Returns tracks: cell column vector — pass-over structs with fields x, y, t
            % alongtrack is a struct with fields (t,x,y)
            arguments
                alongtrack (1,1) struct
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
    end
end
