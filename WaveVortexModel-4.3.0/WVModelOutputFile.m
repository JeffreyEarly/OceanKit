classdef WVModelOutputFile < handle & matlab.mixin.Heterogeneous
    % Organize one NetCDF output file for a WVModel.
    %
    % A `WVModelOutputFile` owns the writable NetCDF handle and one or more
    % output groups. The file is opened when its first group initializes, so
    % `ncfile` may be empty before output begins.
    %
    % Most users create a configured file through `WVModel`:
    %
    % ```matlab
    % outputFile = model.createNetCDFFileForModelOutput("myfile.nc",outputInterval=86400);
    % ```
    %
    % This adds an evenly spaced output group containing the wave-vortex
    % coefficients. For explicit control, first create an empty file
    % description and then add one or more groups:
    %
    % ```matlab
    % outputFile = model.addNewOutputFile("myfile.nc");
    % outputGroup = outputFile.addNewEvenlySpacedOutputGroup( ...
    % "daily",initialTime=model.t,outputInterval=86400);
    % ```
    %
    % An empty output file has no schedule and writes nothing until a group is
    % added. The physical NetCDF file is created lazily at the first scheduled
    % output time, so `ncfile` remains empty beforehand. Close the file through
    % the model or output-file facade when writing is complete.
    %
    %
    % - Topic: Initializing
    % - Topic: Properties
    % - Topic: Internal
    %
    % - Declaration: WVModelOutputFile < handle
    properties (WeakHandle)
        % reference to the WVModel being used
        %
        % - Topic: Properties
        model WVModel
    end

    properties
        % current (or future) path of the NetCDF file
        %
        % - Topic: Properties
        path

        % reference to the NetCDFFile being used for model output
        %
        % This property may be empty if the file is not yet created
        %
        % - Topic: Properties
        ncfile NetCDFFile

        % boolean indicating whether or not the internal structure of the NetCDF file has been created
        %
        % - Topic: Properties
        didInitializeStorage = false

        % time at which the NetCDF file will be created
        %
        % - Topic: Properties
        tInitialize = Inf
    end

    properties (Dependent)
        % array of `WVModelOutputGroup`s that will be written to file
        %
        % - Topic: Properties
        outputGroups

        % name of the current (or future) NetCDF file
        %
        % - Topic: Properties
        filename

        % pass-through of the wvt instance
        %
        % - Topic: Properties
        wvt
    end

    properties (Access=private)
        outputGroupNameMap = configureDictionary("string","WVModelOutputGroup")
        outputGroupNameOutputTimeMap = configureDictionary("string","cell")
    end

    methods
        function self = WVModelOutputFile(model,path,options)
            % initialize a WVModelOutputFile
            %
            % - Topic: Initialization
            % - Declaration: self = WVModelOutputFile(model,path,options)
            % - Parameter model: a WVModel instance
            % - Parameter path: path where the file is to be written
            % - Parameter ncfile: (optional) handle to existing NetCDFFile 
            % - Parameter shouldOverwriteExisting: (optional) whether to overwrite an existing file (default: false).
            % - Returns self: a WVModelOutputFile instance
            arguments
                model WVModel
                path {mustBeText}
                options.ncfile
                options.shouldOverwriteExisting logical = false
            end
            if isfield(options,"ncfile")
                self.ncfile = options.ncfile;
                self.didInitializeStorage = true;
            else
                if options.shouldOverwriteExisting == 1
                    if isfile(path)
                        delete(path);
                    end
                else
                    if isfile(path)
                        error('A file already exists with that name.')
                    end
                end
            end
            self.model = model;
            self.path = path;
        end

        function filename = get.filename(self)
            [~,name,ext] = fileparts(self.path);
            filename = strcat(name,ext);
        end

        function wvt = get.wvt(self)
            wvt = self.model.wvt;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add/remove output groups
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function outputGroups = get.outputGroups(self)
            % return array of output groups associated with this file
            %
            % - Topic: Output groups
            outputGroups = [self.outputGroupNameMap(self.outputGroupNameMap.keys)];
        end

        function names = outputGroupNames(self)
            % retrieve the names of all output group names
            %
            % - Topic: Output groups
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
            end
            arguments (Output)
                names string
            end
            names = self.outputGroupNameMap.keys;
        end

        function val = outputGroupWithName(self,name)
            % retrieve a WVModelOutputGroup by name
            %
            % - Topic: Output groups
            arguments (Input)
                self WVModelOutputFile {mustBeNonempty}
                name {mustBeText,mustBeNonempty}
            end
            arguments (Output)
                val WVModelOutputGroup
            end
            name = string(name);
            if ~isKey(self.outputGroupNameMap,name)
                error('No output group named %s is registered with this file.',name);
            end
            val = self.outputGroupNameMap(name);
        end

        function addOutputGroup(self,outputGroup)
            % add an output group to this file
            %
            % - Topic: Output groups
            arguments
                self WVModelOutputFile {mustBeNonempty}
                outputGroup (1,1) WVModelOutputGroup
            end
            if outputGroup.model ~= self.model
                error('The output group %s was not initialized for this output file''s model.',outputGroup.name);
            end
            if isKey(self.outputGroupNameMap,outputGroup.name)
                registeredGroup = self.outputGroupNameMap(outputGroup.name);
                if registeredGroup == outputGroup
                    return
                end
                error('A different output group named %s is already registered with this file.',outputGroup.name);
            end
            self.outputGroupNameMap(outputGroup.name) = outputGroup;
        end

        function outputGroup = addNewEvenlySpacedOutputGroup(self,name,options)
            % add an evenly-spaced output group to this file
            %
            % - Topic: Output groups
            arguments (Input)
                self WVModelOutputFile
                name {mustBeText}
                options.outputInterval (1,1) double {mustBePositive}
                options.initialTime (1,1) double = -Inf
                options.finalTime (1,1) double = Inf
            end
            arguments (Output)
                outputGroup WVModelOutputGroup
            end
            options.name = name;
            optionCell = namedargs2cell(options);
            outputGroup = WVModelOutputGroupEvenlySpaced(self.model,optionCell{:});
            self.addOutputGroup(outputGroup);
        end

        function addObservingSystem(self,observingSystem)
            % add an observing system to the ouput group (if there is only one group)
            %
            % If there are multiple output groups this will throw an error,
            % as you must decide which group you want to add the observing
            % system to.
            %
            % - Topic: Output groups
            arguments
                self WVModelOutputFile {mustBeNonempty}
                observingSystem WVObservingSystem
            end
            if length(self.outputGroups) ~= 1
                error("This output file has %d output groups, thus you must add the observing system the group you want.",length(self.outputGroups))
            end
            self.outputGroups(1).addObservingSystem(observingSystem);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % returns a unique, ordered array of the aggregate output times during the requested integration period.
            %
            % This will be called exactly once by the model before an
            % integration begins.
            %
            % - Topic: Internal
            arguments (Input)
                self WVModelOutputFile
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            t = [];
            outputGroups_ = self.outputGroups;
            for iGroup = 1:length(outputGroups_)
                t_group = outputGroups_(iGroup).outputTimesForIntegrationPeriod(initialTime,finalTime);
                self.outputGroupNameOutputTimeMap{outputGroups_(iGroup).name} = t_group;
                t = cat(1,t,t_group);
            end
            t = unique(t,"sorted");
            if self.didInitializeStorage == false && ~isempty(t)
                self.tInitialize = t(1);
            end
        end

        function writeTimeStepToOutputFile(self,t)
            % tells the output groups to write data at time t
            %
            % This will only be called if the
            % `outputTimesForIntegrationPeriod` previously returned this
            % time `t`. Thus, the actual NetCDF file will be created if
            % it does not yet exist, and then the output groups can write
            % to file.
            %
            % - Topic: Internal
            arguments
                self WVModelOutputFile
                t (1,1) double
            end

            try
                % 1) initialize the netcdf file if necessary
                if self.didInitializeStorage == false && t == self.tInitialize
                    self.initializeOutputFile();
                end

                % 2) inform the appropriate groups that they need to write a time step.
                outputGroupNames = self.outputGroupNameOutputTimeMap.keys;
                didWriteToFile = false;
                for i = 1:length(outputGroupNames)
                    t_group = self.outputGroupNameOutputTimeMap{outputGroupNames(i)};
                    if ~isempty(t_group) && t == t_group(1)
                        self.outputGroupWithName(outputGroupNames(i)).writeTimeStepToNetCDFFile(self.ncfile,t);
                        t_group(1) = [];
                        self.outputGroupNameOutputTimeMap{outputGroupNames(i)} = t_group;
                        didWriteToFile = true;
                    end
                end
                if didWriteToFile
                    self.ncfile.sync();
                end
            catch exception
                self.closeNetCDFFile();
                rethrow(exception)
            end
        end

        function initializeOutputFile(self)
            % tells the output groups to initialize themselves in the NetCDF file
            %
            % One noteworthy piece of logic handle by this function: we
            % want all files to able to re-initialize the model, so we make
            % sure that the NetCDFFile has all the required properties of
            % the transform. However, we are writing time series data, so
            % we need to exclude `t`. Additionally, if we happen to also be
            % writing the wave-vortex coefficients, we can neglect writing
            % those to file.
            %
            % Initialization is transactional. If any group or observing
            % system fails, the new file is closed and deleted and every
            % storage object returns to its uninitialized state.
            %
            % - Topic: Internal
            if self.observingSystemWillWriteWaveVortexCoefficients == true
                properties = setdiff(self.wvt.requiredProperties,{'Ap','Am','A0','t'});
            else
                properties = setdiff(self.wvt.requiredProperties,{'t'});
            end
            % in theory we already removed the file if the user requested
            if isfile(self.path)
                error('A file already exists at this path.');
            end
            newFile = NetCDFFile.empty(0,0);
            try
                newFile = self.wvt.writeToFile(self.path,properties{:},shouldOverwriteExisting=false,shouldAddRequiredProperties=false);
                newFile.addAttribute('WVModelIsDynamicsLinear',uint8(self.model.isDynamicsLinear));
                outputGroups_ = self.outputGroups;
                for iGroup = 1:length(outputGroups_)
                    outputGroups_(iGroup).initializeOutputGroup(newFile);
                end
            catch exception
                outputGroups_ = self.outputGroups;
                for iGroup = 1:length(outputGroups_)
                    outputGroups_(iGroup).group = NetCDFGroup.empty(0,0);
                    outputGroups_(iGroup).incrementsWrittenToGroup = 0;
                    outputGroups_(iGroup).timeOfLastIncrementWrittenToGroup = -Inf;
                    outputGroups_(iGroup).didInitializeStorage = false;
                end
                if ~isempty(newFile) && isvalid(newFile) && ~isempty(newFile.id)
                    newFile.close();
                end
                self.ncfile = NetCDFFile.empty(0,0);
                self.didInitializeStorage = false;
                if isfile(self.path)
                    delete(self.path);
                end
                rethrow(exception)
            end
            self.ncfile = newFile;
            self.didInitializeStorage = true;
        end

        function bool = observingSystemWillWriteWaveVortexCoefficients(self)
            % A simple check to see if one of the observing systems will be writing wave-vortex coefficients
            %
            % - Topic: Internal
            outputGroups_ = self.outputGroups;
            bool = false;
            for iGroup = 1:length(outputGroups_)
                observingSystems = outputGroups_(iGroup).observingSystems;
                for iObs = 1:length(observingSystems)
                    if isa(observingSystems(iObs),'WVEulerianFields')
                        bool = bool | all(ismember(intersect({'Ap','Am','A0'},self.wvt.variableNames),observingSystems(iObs).netCDFOutputVariables));
                    end
                end
            end
        end

        function recordNetCDFFileHistory(self,options)
            % tells the output groups to log this time step in the NetCDF history
            %
            % - Topic: Internal
            arguments
                self WVModelOutputFile {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end

            if ~isempty(self.ncfile)
                arrayfun( @(outputGroup) outputGroup.recordNetCDFFileHistory(didBlowUp=options.didBlowUp), self.outputGroups);
            end
        end


        function closeNetCDFFile(self)
            % closes the netcdf file after informing the output groups
            %
            % The file handle is released even if an output-group cleanup
            % notification fails.
            %
            % - Topic: Internal
            if ~isempty(self.ncfile)
                ownedFile = self.ncfile;
                self.ncfile = NetCDFFile.empty(0,0);
                cleanup = onCleanup(@()WVModelOutputFile.closeIfOpen(ownedFile));
                outputGroups_ = self.outputGroups;
                for iGroup = 1:length(outputGroups_)
                    outputGroups_(iGroup).closeNetCDFFile();
                end
                WVModelOutputFile.closeIfOpen(ownedFile);
                clear cleanup
            end
        end

    end

    methods (Static)
        function outputFile = modelOutputFileFromFile(file,model)
            % create a WVModelOutputFile from an existing NetCDFFile
            %
            % - Topic: Initialization
            arguments
                file NetCDFFile {mustBeNonempty}
                model WVModel {mustBeNonempty}
            end
            outputFile = WVModelOutputFile(model,file.path,ncfile=file);

            for iGroup=1:length(file.groups)
                group = file.groups(iGroup);
                if group.hasAttributeWithName("AnnotatedClass")
                    className = group.attributes('AnnotatedClass');
                    if exist(className,'class') && ismember("WVModelOutputGroup",superclasses(className))
                        outputFile.addOutputGroup(WVModelOutputGroup.modelOutputGroupFromGroup(group,model));                        
                    end
                end
            end
            WVModelOutputFile.canonicalizeIntegratedObservingSystems(outputFile);
        end
    end

    methods (Static, Access=private)
        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end

        function canonicalizeIntegratedObservingSystems(outputFile)
            outputGroups_ = outputFile.outputGroups;
            processedNames = configureDictionary("string","logical");
            for iGroup = 1:length(outputGroups_)
                for iObserver = 1:length(outputGroups_(iGroup).observingSystems)
                    observer = outputGroups_(iGroup).observingSystems(iObserver);
                    observerName = string(observer.name);
                    if observer.nFluxComponents == 0 || isKey(processedNames,observerName)
                        continue
                    end
                    processedNames(observerName) = true;

                    [candidates,candidateGroups] = WVModelOutputFile.integratedObserverCandidates(outputGroups_,observerName);
                    WVModelOutputFile.validateEquivalentObserverConfigurations(candidates);
                    if isa(observer,'WVCoefficients')
                        canonicalObserver = outputFile.model.wvCoefficientFluxedObservingSystem();
                        if isempty(canonicalObserver)
                            error('The restored output contains coefficients, but the model is configured for linear dynamics.');
                        end
                        canonicalObserver.absTolerance = candidates{1}.absTolerance;
                    else
                        candidateTimes = cellfun(@(group) group.timeOfLastIncrementWrittenToGroup,candidateGroups);
                        timeTolerance = 8*eps(max(1,abs(outputFile.model.t)));
                        currentStateIndices = find(abs(candidateTimes-outputFile.model.t) <= timeTolerance);
                        if isempty(currentStateIndices)
                            error('The integrated observing system %s has no saved state at the model restart time t=%g.',observerName,outputFile.model.t);
                        end
                        canonicalObserver = candidates{currentStateIndices(1)};
                        canonicalState = canonicalObserver.initialConditions();
                        for iCandidate = reshape(currentStateIndices(2:end),1,[])
                            if ~isequaln(candidates{iCandidate}.initialConditions(),canonicalState)
                                error('The output groups contain inconsistent saved states for observing system %s at t=%g.',observerName,outputFile.model.t);
                            end
                        end
                        outputFile.model.addFluxedObservingSystem(canonicalObserver);
                    end

                    for iCandidate = 1:length(candidates)
                        candidateGroup = candidateGroups{iCandidate};
                        candidateIndex = find(arrayfun(@(existing) existing == candidates{iCandidate},candidateGroup.observingSystems),1);
                        candidateGroup.observingSystems(candidateIndex) = canonicalObserver;
                    end
                end
            end
        end

        function [candidates,candidateGroups] = integratedObserverCandidates(outputGroups,observerName)
            maximumCandidateCount = sum(arrayfun(@(group)length(group.observingSystems),outputGroups));
            candidates = cell(1,maximumCandidateCount);
            candidateGroups = cell(1,maximumCandidateCount);
            candidateCount = 0;
            for iGroup = 1:length(outputGroups)
                for iObserver = 1:length(outputGroups(iGroup).observingSystems)
                    candidate = outputGroups(iGroup).observingSystems(iObserver);
                    if candidate.nFluxComponents > 0 && string(candidate.name) == observerName
                        candidateCount = candidateCount+1;
                        candidates{candidateCount} = candidate;
                        candidateGroups{candidateCount} = outputGroups(iGroup);
                    end
                end
            end
            candidates = candidates(1:candidateCount);
            candidateGroups = candidateGroups(1:candidateCount);
        end

        function validateEquivalentObserverConfigurations(candidates)
            reference = candidates{1};
            requiredProperties = feval(strcat(class(reference),'.classRequiredPropertyNames'));
            for iCandidate = 2:length(candidates)
                candidate = candidates{iCandidate};
                if ~strcmp(class(candidate),class(reference))
                    error('Integrated observing systems named %s must have the same class in every output group.',reference.name);
                end
                for iProperty = 1:length(requiredProperties)
                    propertyName = requiredProperties{iProperty};
                    if ~isequaln(candidate.(propertyName),reference.(propertyName))
                        error('Integrated observing systems named %s have inconsistent persisted configuration for %s.',reference.name,propertyName);
                    end
                end
            end
        end
    end
end
