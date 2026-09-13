classdef WVModelOutputGroup < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Schedule observing systems into one NetCDF output group.
    %
    % A `WVModelOutputGroup` contains one or more observing systems and
    % defines the model times at which their state is written.
    %
    % The simplest output group is the
    % [`WVModelOutputGroupEvenlySpaced`](/classes/model-output/wvmodeloutputgroupevenlyspaced/)
    % which, as the name suggests, writes outputs at an evenly spaced
    % interval.
    %
    % Separate groups allow different observing systems to use different
    % schedules or bounded output windows within the same file. This matters
    % when an observing system does not sample evenly, such as a satellite
    % along-track simulator, or when one diagnostic needs a shorter,
    % higher-frequency window than the rest of a long model run. For example,
    % a mooring may need to resolve the buoyancy frequency, while a tracer
    % experiment may run for only 24 hours in the middle of the simulation.
    %
    % ### Usage
    %
    % ```matlab
    % outputFile = model.addNewOutputFile("myfile.nc");
    % outputGroup = WVModelOutputGroupEvenlySpaced(model,name="high-temporal-resolution",initialTime=wvt.t,outputInterval=wvt.inertialPeriod/20);
    % outputFile.addOutputGroup(outputGroup);
    % ```
    %
    % 
    % - Topic: Initializing
    % - Topic: Properties
    % - Topic: Observing systems
    % - Topic: Required subclass overrides
    % - Topic: Internal
    %
    % - Declaration: WVModelOutputGroup < handle
    properties (WeakHandle)
        % Reference to the WVModel being used
        %
        % - Topic: Properties
        model WVModel

        % Reference to the NetCDFGroup being used for model output
        % - Topic: Writing to NetCDF files
        % Empty indicates no file output. The output group creates the
        % NetCDFGroup, but the NetCDFFile owns it, hence a WeakHandle.
        group NetCDFGroup
    end

    properties
        % name of the current (or future) group in the NetCDF file
        %
        % - Topic: Properties
        name string

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        incrementsWrittenToGroup (1,1) uint64 = 0

        % output index of the current/most recent step.
        % - Topic: Integration
        % If stepsTaken=0, outputIndex=1 means the initial conditions get written at index 1
        timeOfLastIncrementWrittenToGroup (1,1) double = -Inf

        % boolean indicating whether or not the internal structure of the NetCDF file has been created
        %
        % - Topic: Properties
        didInitializeStorage = false
    end

    properties
        % array of WVObservingSystem that will be written to the group
        %
        % - Topic: Observing systems
        observingSystems
    end

    methods
        function self = WVModelOutputGroup(model,options)
            % initialize a WVModelOutputGroup
            %
            % - Topic: Initialization
            % - Declaration: self = WVModelOutputGroup(model,path,options)
            % - Parameter model: a WVModel instance
            % - Parameter name: name of the group
            % - Returns self: a WVModelOutputFile instance
            arguments
                model WVModel
                options.name {mustBeText}
            end
            if ~isfield(options,"name")
                error("You must specify an output group name");
            end
            self.model = model;
            self.name = options.name;
            self.observingSystems = WVObservingSystem.empty(1,0);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add/remove observing systems
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addObservingSystem(self,observingSystem)
            % add an observing system to this file
            %
            % - Topic: Observing systems
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                observingSystem WVObservingSystem
            end
            if self.didInitializeStorage
                error('Storage already initialized! You cannot add a new observing system after the storage has been initialized.');
            end

            registeredSystems = self.observingSystems;
            fluxedSystems = WVObservingSystem.empty(1,0);
            coefficientSystems = {};
            coefficientTolerances = {};
            for iObs = 1:length(observingSystem)
                observer = observingSystem(iObs);
                if observer.model ~= self.model
                    error('The observing system %s was not initialized for this output group''s model.',observer.name);
                end
                if isa(observer,'WVCoefficients')
                    % Model construction creates the canonical coefficient observer
                    % before persisted output groups are restored.
                    canonicalCoefficients = self.model.wvCoefficientFluxedObservingSystem();
                    if ~isempty(canonicalCoefficients) && canonicalCoefficients ~= observer
                        coefficientSystems{end+1} = canonicalCoefficients;
                        coefficientTolerances{end+1} = observer.absTolerance;
                        observer = canonicalCoefficients;
                    end
                end
                if any(strcmp(string({registeredSystems.name}),string(observer.name)))
                    error('An observing system named %s already exists in this output group.',observer.name);
                end
                if isempty(registeredSystems)
                    registeredSystems = observer;
                else
                    registeredSystems(end+1) = observer;
                end
                if observer.nFluxComponents > 0
                    fluxedSystems(end+1) = observer;
                end
            end
            if ~isempty(fluxedSystems)
                self.model.addFluxedObservingSystem(fluxedSystems);
            end
            for iCoefficient = 1:length(coefficientSystems)
                coefficientSystems{iCoefficient}.absTolerance = coefficientTolerances{iCoefficient};
            end
            self.observingSystems = registeredSystems;
        end

        function removeObservingSystem(self, observingSystem)
            % remove an observing system to this file
            %
            % - Topic: Observing systems
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                observingSystem WVObservingSystem
            end

            registeredSystems = self.observingSystems;
            removeMask = false(size(registeredSystems));
            fluxedSystems = WVObservingSystem.empty(1,0);
            for iObs = 1:length(observingSystem)
                observer = observingSystem(iObs);
                if observer.model ~= self.model
                    error('The observing system %s was not initialized for this output group''s model.',observer.name);
                end
                sameHandle = cellfun(@(existing) existing == observer,num2cell(registeredSystems));
                if ~any(sameHandle)
                    error('The observing system %s is not registered with this output group.',observer.name);
                end
                removeMask = removeMask | sameHandle;
                if observer.nFluxComponents > 0
                    fluxedSystems(end+1) = observer;
                end
            end
            if ~isempty(fluxedSystems)
                self.model.removeFluxedObservingSystem(fluxedSystems);
            end
            registeredSystems(removeMask) = [];
            self.observingSystems = registeredSystems;
        end

        function observingSystem = observingSystemWithName(self,name)
            % retrieve an observing system by name
            %
            % - Topic: Observing systems
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                name {mustBeText}
            end
            idx = find(strcmp(string({self.observingSystems.name}),string(name)),1);
            if isempty(idx)
                error('No observing system named %s is registered with this output group.',name);
            end
            observingSystem = self.observingSystems(idx);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % returns a unique, ordered array of the aggregate output times during the requested integration period.
            %
            % Any new subclass of `WVModelOutputGroup` must override this
            % methods and return the appropriate output times.
            %
            % - Topic: Required subclass overrides
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            t = [];
        end

        function initializeOutputGroup(self,ncfile)
            % initializes a new output group in the NetCDF file
            %
            % This will only be called only once. This creates a new group,
            % adds the required properties, creates a time dimension, then
            % tells the observing systems to also initialize their storage.
            %
            % - Topic: Internal
            arguments (Input)
                self WVModelOutputGroup {mustBeNonempty}
                ncfile NetCDFGroup {mustBeNonempty}
            end
            if self.didInitializeStorage
                error('Storage already initialized!');
            end
            try
                self.group = ncfile.addGroup(self.name);
                self.writeToGroup(self.group,self.propertyAnnotationWithName(self.requiredProperties));

                varAnnotation = self.model.wvt.propertyAnnotationWithName('t');
                varAnnotation.attributes('units') = varAnnotation.units;
                varAnnotation.attributes('long_name') = varAnnotation.description;
                varAnnotation.attributes('standard_name') = 'time';
                varAnnotation.attributes('long_name') = 'time';
                varAnnotation.attributes('units') = 'seconds since 1970-01-01 00:00:00';
                varAnnotation.attributes('axis') = 'T';
                varAnnotation.attributes('calendar') = 'standard';
                self.group.addDimension(varAnnotation.name,length=Inf,type="double",attributes=varAnnotation.attributes);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).initializeStorage(self.group);
                end
            catch exception
                self.group = NetCDFGroup.empty(0,0);
                self.incrementsWrittenToGroup = 0;
                self.timeOfLastIncrementWrittenToGroup = -Inf;
                self.didInitializeStorage = false;
                rethrow(exception)
            end

            self.incrementsWrittenToGroup = 0;
            self.didInitializeStorage = true;
        end

        function writeTimeStepToNetCDFFile(self,ncfile,t)
            % writes data at time t
            %
            % This is called by the `WVModelOutputFile` when the model
            % reaches time `t`. The new time is written to file, adn the
            % observing systems are also told to write to file.
            %
            % - Topic: Internal
            arguments
                self WVModelOutputGroup
                ncfile NetCDFFile
                t double
            end
            if ~self.didInitializeStorage
                self.initializeOutputGroup(ncfile);
            end
            if ( ~isempty(self.group) && t > self.timeOfLastIncrementWrittenToGroup )
                outputIndex = self.incrementsWrittenToGroup + 1;

                self.group.variableWithName('t').setValueAlongDimensionAtIndex(t,'t',outputIndex);

                for iObs = 1:length(self.observingSystems)
                    self.observingSystems(iObs).writeTimeStepToFile(self.group,outputIndex);
                end

                self.incrementsWrittenToGroup = outputIndex;
                self.timeOfLastIncrementWrittenToGroup = t;
            end
        end

        function recordNetCDFFileHistory(self,options)
            % losg this time step in the NetCDF history
            %
            % - Topic: Internal
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                options.didBlowUp {mustBeNumeric} = 0
            end
            if isempty(self.group)
                return
            end

            if options.didBlowUp == 1
                a = sprintf('%s: wrote %d time points to file. Terminated due to model blow-up.',datetime('now'),self.incrementsWrittenToGroup);
            else
                a = sprintf('%s: wrote %d time points to file',datetime('now'),self.incrementsWrittenToGroup);
            end
            if isKey(self.group.attributes,'history')
                history = reshape(self.group.attributes('history'),1,[]);
                history =cat(2,squeeze(history),a);
            else
                history = a;
            end
            self.group.addAttribute('history',history);
        end

        function closeNetCDFFile(self)
            % notification that the NetCDF file will close
            %
            % This gives the output group an opportunity to display some
            % relevant data or do other necessary clean up.
            %
            % - Topic: Internal
            if ~isempty(self.group)
                fprintf('Ending simulation. Wrote %d time points to %s group\n',self.incrementsWrittenToGroup,self.name);
            end
        end

        function initObservingSystemsFromGroup(self,outputGroup)
            % asks the output group to load the observing systems in the NetCDF file
            %
            % Called by the static method `modelOutputGroupFromGroup` during the init from file process, this asks the output group to load the observing systems in the NetCDF file 
            %
            % - Topic: Internal
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                outputGroup NetCDFGroup {mustBeNonempty}
            end

            f = @(className,group) feval(strcat(className,'.observingSystemFromGroup'),group, self.model, self);
            vars = CAAnnotatedClass.propertyValuesFromGroup(outputGroup,{"observingSystems"},classConstructor=f);
            self.observingSystems = vars.observingSystems;
        end

    end

    methods (Static)
        function outputGroup = modelOutputGroupFromGroup(group,model)
            %initialize a WVModelOutputGroup instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: outputGroup = modelOutputGroupFromGroup(group,model)
            % - Parameter group: the NetCDFGroup to be used
            % - Parameter model: the WVModel to be used
            % - Returns outputGroup: a new instance of WVModelOutputGroup
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
            end
            className = group.attributes('AnnotatedClass');
            requiredProperties = feval(strcat(className,'.classRequiredPropertyNames'));
            requiredProperties(ismember(requiredProperties,'observingSystems')) = [];
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,requiredProperties);
            if isempty(vars)
                outputGroup = feval(className,model);
            else
                options = namedargs2cell(vars);
                outputGroup = feval(className,model,options{:});
            end
            outputGroup.group = group;
            nPoints = group.dimensionWithName("t").nPoints;
            outputGroup.incrementsWrittenToGroup = nPoints;
            if nPoints > 0
                outputGroup.timeOfLastIncrementWrittenToGroup = group.readVariablesAtIndexAlongDimension('t',nPoints,'t');
            end
            outputGroup.initObservingSystemsFromGroup(group);
            outputGroup.didInitializeStorage = true;
        end
    end
end
