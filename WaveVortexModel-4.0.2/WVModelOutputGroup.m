classdef WVModelOutputGroup < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % A WVModelOutputGroup represents a group of observing system to be written to file at particular output times
    %
    % A `WVModelOutputGroup` encapsulates a netcdf group with particular
    % output times `t`; it has one or more observing systems that get written to the group at those times.
    %
    % The simplest output group is the
    % [`WVModelOutputGroupEvenlySpaced`](/classes/model-output/wvmodeloutputgroupevenlyspaced/)
    % which, as the name suggests, writes outputs at an evenly spaced
    % interval.
    %
    % The reason for this abstraction is that some observing systems do not
    % have evenly spaced output intervals, e.g., the
    % [AlongTrackSimulator](https://satmapkit.github.io/AlongTrackSimulator/),
    % and it is often the case that one might want to sample the model for
    % a short interval at higher frequency, e.g., sampling a mooring at
    % high frequency to resolve the buoyancy frequency, or running a tracer
    % experiment for 24 hours in the middle of a long model run.
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
    % - Declaration: WVModelOutputFile < handle
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

            for iObs = 1:length(observingSystem)
                anObservingSystem = observingSystem(iObs);
                if anObservingSystem.wvt ~= self.model.wvt
                    error('This observing system was not initialized with the same wvt that it is being added to!')
                end

                if ~isempty(find(strcmp({self.observingSystems.name}, anObservingSystem.name), 1))
                    error('An observing system named %s already exists.\n',anObservingSystem.name)
                end

                self.observingSystems(end+1) = anObservingSystem;
                if anObservingSystem.nFluxComponents > 0
                    self.model.addFluxedObservingSystem(anObservingSystem);
                end
            end
        end

        function removeObservingSystem(self, observingSystem)
            % remove an observing system to this file
            %
            % - Topic: Observing systems
            arguments
                self WVModelOutputGroup {mustBeNonempty}
                observingSystem WVObservingSystem
            end

            for iObs = 1:length(observingSystem)
                anObservingSystem = observingSystem(iObs);

                % Verify that the observing system belongs to the same wvt
                if anObservingSystem.wvt ~= self
                    error('This observing system does not belong to the same wvt!');
                end

                % Find index of the observing system with the matching name
                idx = find(strcmp({self.observingSystems.name}, anObservingSystem.name));
                if isempty(idx)
                    error('No observing system named %s exists to remove.', anObservingSystem.name);
                end

                % Remove the observing system from the collection
                self.observingSystems(idx) = [];

                % If the system includes flux components, remove it from the fluxed systems in the model
                if anObservingSystem.nFluxComponents > 0
                    self.model.removeFluxedObservingSystem(anObservingSystem);
                end
            end
        end

        function observingSystem = observingSystemWithName(self,name)
            % retrieve an observing system by name
            %
            % - Topic: Observing systems
            idx = find(strcmp({self.observingSystems.name}, name),1);
            if isempty(idx)
                error('No observing system named %s exists to remove.', anObservingSystem.name);
            else
                observingSystem = self.observingSystems(idx);
            end
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
            self.addObservingSystem(vars.observingSystems);
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