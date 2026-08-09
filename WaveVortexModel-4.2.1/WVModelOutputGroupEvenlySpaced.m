classdef WVModelOutputGroupEvenlySpaced < WVModelOutputGroup
    %A WVModelOutputGroupEvenlySpaced is an evenly spaced output group with optional start and end times
    %
    % This is the primary `WVModelOutputGroup` subclass which implements
    % evenly spaced output, standard in numerical modeling. It also
    % includes optional `initialTime` and `finalTime` properties, which
    % allow you to customize when this group is active in model time.
    %
    % 
    % - Topic: Initializing
    % - Topic: Properties
    %
    % - Declaration: WVModelOutputGroupEvenlySpaced < [WVModelOutputGroup](/classes/model-output/wvmodeloutputgroup/)

    properties
        % model output interval (seconds)
        %
        % The model output interval written to the group
        %
        % - Topic: Properties
        outputInterval = [] % (1,1) double

        % initial model time that the output group is active (seconds)
        %
        % This optional properties determines when the output group will
        % become active. By default it is set to `-Inf`.
        %
        % - Topic: Properties
        initialTime

        % final model time that the output group is active (seconds)
        %
        % This optional properties determines when the output group will
        % de-actives. By default it is set to `Inf`, indicating that it
        % will always write to file.
        %
        % - Topic: Properties
        finalTime
    end

    methods
        function self = WVModelOutputGroupEvenlySpaced(model,options)
            % initialize a WVModelOutputGroupEvenlySpaced
            %
            % 
            % - Topic: Initialization
            % - Declaration: self = WVModelOutputGroupEvenlySpaced(model,options)
            % - Parameter model: a WVModel instance
            % - Parameter name: name of the group
            % - Parameter outputInterval: model output interval (seconds)
            % - Parameter initialTime: (optional) initial model time that the output group is active (seconds), default `-Inf`
            % - Parameter finalTime: (optional) final model time that the output group is active (seconds), default `Inf`
            % - Returns self: a WVModelOutputGroupEvenlySpaced instance
            arguments
                model WVModel
                options.name {mustBeText}
                options.outputInterval (1,1) double {mustBePositive}
                options.initialTime (1,1) double = -Inf
                options.finalTime (1,1) double = Inf
            end
            self@WVModelOutputGroup(model,name=options.name);
            if ~isfield(options,"name")
                error("You must specify an output group name");
            end
            if ~isfield(options,"outputInterval")
                error("You must specify an output interval");
            end
            self.outputInterval = options.outputInterval;
            if options.initialTime == -Inf
                options.initialTime = model.wvt.t;
            end
            self.initialTime = options.initialTime;
            self.finalTime = options.finalTime;
        end

        function aString = description(self)
            if self.finalTime == Inf
                aString = "evenly spaced output with an interval of " + string(self.outputInterval) + " starting at time t=" + string(self.initialTime) + " and continuing indefinitely";
            else
                aString = "evenly spaced output with an interval of " + string(self.outputInterval) + " start at time t=" + string(self.initialTime) + "and ending at time t=" + string(self.finalTime);
            end
        end

        function t = outputTimesForIntegrationPeriod(self,initialTime,finalTime)
            % Return unwritten output times on the group's fixed time lattice.
            %
            % Times are generated from `initialTime + k*outputInterval`,
            % bounded by the group and integration windows, and strictly
            % later than the last written increment. Anchoring every call
            % to the original lattice prevents segmented-run timing drift.
            %
            % - Topic: Integration
            arguments (Input)
                self WVModelOutputGroup
                initialTime (1,1) double
                finalTime (1,1) double
            end
            arguments (Output)
                t (:,1) double
            end
            lowerTime = max(initialTime,self.initialTime);
            upperTime = min(finalTime,self.finalTime);
            if upperTime < lowerTime
                t = zeros(0,1);
                return
            end

            lowerIndex = (lowerTime-self.initialTime)/self.outputInterval;
            lowerIndexTolerance = 8*eps(max(1,abs(lowerIndex)));
            if abs(lowerIndex-round(lowerIndex)) <= lowerIndexTolerance
                lowerIndex = round(lowerIndex);
            end
            firstIndex = max(0,ceil(lowerIndex));

            if self.timeOfLastIncrementWrittenToGroup ~= -Inf
                lastWrittenIndex = (self.timeOfLastIncrementWrittenToGroup-self.initialTime)/self.outputInterval;
                lastWrittenIndexTolerance = 8*eps(max(1,abs(lastWrittenIndex)));
                if abs(lastWrittenIndex-round(lastWrittenIndex)) <= lastWrittenIndexTolerance
                    lastWrittenIndex = round(lastWrittenIndex);
                end
                firstIndex = max(firstIndex,floor(lastWrittenIndex)+1);
            end

            lastIndex = (upperTime-self.initialTime)/self.outputInterval;
            lastIndexTolerance = 8*eps(max(1,abs(lastIndex)));
            if abs(lastIndex-round(lastIndex)) <= lastIndexTolerance
                lastIndex = round(lastIndex);
            end
            lastIndex = floor(lastIndex);
            if lastIndex < firstIndex
                t = zeros(0,1);
            else
                t = self.initialTime+(firstIndex:lastIndex).'*self.outputInterval;
            end
        end

    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'observingSystems','name','outputInterval','initialTime','finalTime'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAObjectProperty('observingSystems','array of WVObservingSystem objects');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of output group');
            propertyAnnotations(end+1) = CANumericProperty('outputInterval', {}, 's','output interval');
            propertyAnnotations(end+1) = CANumericProperty('initialTime', {}, 's','model time of first allowed write to file');
            propertyAnnotations(end+1) = CANumericProperty('finalTime', {}, 's','model time of last allowed write to file');
        end
    end
end
