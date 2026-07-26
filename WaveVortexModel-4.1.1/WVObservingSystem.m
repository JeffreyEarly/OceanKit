classdef WVObservingSystem < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    %A WVObservingSystem is an abstract class that defines different ways of observing the model
    %
    %
    %
    % 
    % - Topic: Initializing
    % - Topic: Properties
    % - Topic: Required subclass overrides
    % - Topic: Required subclass overrides for integrated (fluxed) components
    % - Topic: Internal
    %
    % - Declaration: WVObservingSystem < handle
    properties (WeakHandle)
        % reference to the WVModel being used
        %
        % - Topic: Properties
        model WVModel
    end

    properties (GetAccess=public, SetAccess=protected)
        % name of the observing system
        %
        % - Topic: Properties
        name

        % number of components that need to be integrated in time.
        %
        % Setting a value greater than zero will require that you
        % implement,
        %   -absErrorTolerance
        %   -initialConditions
        %   -fluxAtTime
        %   -updateIntegratorValues
        %
        % - Topic: Properties
        nFluxComponents uint8 = 0

        % reference to the WVModel being used
        %
        % - Topic: Properties
        wvt % this should be weak
    end

    methods
        function self = WVObservingSystem(model,name)
            %create a new observing system
            %
            % This class is intended to be subclassed, so it generally
            % assumed that this initialization will not be called directly.
            %
            % - Topic: Initializing
            % - Declaration: self = WVObservingSystem(wvt,name)
            % - Parameter wvt: the WVTransform instance
            % - Parameter name: name of the observing system
            % - Returns self: a new instance of WVObservingSystem
            arguments
                model WVModel
                name {mustBeText}
            end
            % Do we actually want to inherit the properties from the
            % WVTransform? I'm not sure. I think this should be optional.
            % If an OS does, then its output can go in the wave-vortex
            % group.
            % self@CAAnnotatedClass();
            self.model = model;
            self.name = name;
            self.wvt = model.wvt;
        end

        function aString = description(self)
            aString = "An instance of " + class(self) + " named " + self.name;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integrated variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function nArray = lengthOfFluxComponents(self)
            % return an array containing the numel of each flux component.
            %
            % - Topic: Required subclass overrides for integrated (fluxed) components
            nArray = [];
        end

        function Y0 = absErrorTolerance(self)
            % return a cell array of the absolute tolerances of the
            % variables being integrated. You can pass either scalar
            % values, or an array of the same size as the variable.
            %
            % this will only be called when the time-stepping is run with
            % an adaptive integrator.
            %
            % - Topic: Required subclass overrides for integrated (fluxed) components
            Y0 = {};
        end

        function Y0 = initialConditions(self)
            % return a cell array of variables that need to be integrated
            %
            % - Topic: Required subclass overrides for integrated (fluxed) components
            Y0 = {};
        end

        function F = fluxAtTime(self,t,y0)
            % return a cell array of the flux of the variables being
            % integrated. You may want to call -updateIntegratorValues.
            %
            % - Topic: Required subclass overrides for integrated (fluxed) components
            F = {};
        end

        function updateIntegratorValues(self,t,y0)
            % passes updated values of the variables being integrated.
            % Called 1) during fluxAtTime and 2) during interpolation
            %
            % - Topic: Required subclass overrides for integrated (fluxed) components
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read and write to file
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function initializeStorage(self,group)
            % called once to allow the observing system to initialize its storage space in the NetCDFGroup
            %
            % Any observing system must implement this method to initialize
            % the NetCDFGroup with the appropriate attributes, dimensions
            % and variables.
            %
            % - Topic: Required subclass overrides
        end

        function writeTimeStepToFile(self,group,outputIndex)
            % called at each time for the observing system to write to file
            %
            % Any observing system must implement this method to initialize
            % and write to the NetCDFGroup at the particular outputIndex.
            %
            % - Topic: Required subclass overrides
        end
    end

    methods (Static)
        function os = observingSystemFromGroup(group,model, outputGroup)
            %initialize a WVObservingSystem instance from NetCDF file
            %
            % Subclasses to should override this method to enable model
            % restarts. This method works in conjunction with -writeToFile
            % to provide restart capability.
            %
            % - Topic: Initialization
            % - Declaration: os = observingSystemFromGroup(group,wvt)
            % - Parameter model: the WVModel to be used
            % - Returns os: a new instance of WVObservingSystem
            arguments
                group NetCDFGroup {mustBeNonempty}
                model WVModel {mustBeNonempty}
                outputGroup WVModelOutputGroup
            end
            className = group.attributes('AnnotatedClass');
            vars = CAAnnotatedClass.requiredPropertiesFromGroup(group);
            if isempty(vars)
                os = feval(className,model);
            else
                options = namedargs2cell(vars);
                os = feval(className,model,options{:});
            end
        end
    end
end