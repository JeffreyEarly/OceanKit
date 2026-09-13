classdef WVVariableAnnotation < CANumericProperty
    % Describe a variable computed from a WVTransform.
    % 
    % In addition to adding a name, description and detailed description of
    % a given variable, you also specify its dimensions, units, and whether
    % or whether it has an imaginary part. These annotations are used for
    % online documentation and for writing to NetCDF files.
    %
    % Setting the two properties `isVariableWithLinearTimeStep` and
    % `isVariableWithNonlinearTimeStep` are important for determining
    % how the variable is cached, and when it is saved to a NetCDF file.
    %
    % A matching Markdown sidecar named for the variable may provide the
    % longer mathematical or scientific description. `CAPropertyAnnotation`
    % locates that sidecar when detailed documentation is requested, and the
    % website builder merges the canonical sidecar into each generated class
    % reference that exposes the variable. This keeps equations and tables out
    % of constructor calls without disconnecting them from the annotation.
    %
    % - Declaration: classdef WVVariableAnnotation < WVAnnotation
    properties
        % WVOperation responsible for computing this variable
        %
        % This property will be automatically populated when the variable
        % annotation is passed to the WVOperation.
        % - Topic: Properties
        modelOp

        % boolean indicating whether the variable changes value with a linear time step
        %
        % This information is used when caching variables and when writing
        % to NetCDF file.
        % - Topic: Properties
        isVariableWithLinearTimeStep logical = true

        % boolean indicating whether the variable changes value with a non-linear time step
        %
        % This information is used when caching variables and when writing
        % to NetCDF file.
        % - Topic: Properties
        isVariableWithNonlinearTimeStep logical = true

        % boolean indicating whether the variable depends on Ap, Am, or A0
        %
        % This information is used when caching variables and when writing
        % to NetCDF file.
        % - Topic: Properties
        isDependentOnApAmA0 logical = true
    end

    methods
        function self = WVVariableAnnotation(name,dimensions,units,description,options)
            % create a new instance of WVVariableAnnotation
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: variableAnnotation = WVVariableAnnotation(name,dimensions,units,description,options)
            % - Parameter name: name of the variable
            % - Parameter dimensions: ordered list of the dimensions, or empty cell array
            % - Parameter units: abbreviated SI units of the variable
            % - Parameter description: short description of the variable
            % - Parameter isComplex: (optional) indicates whether the variable has an imaginary part (default 0)
            % - Parameter detailedDescription: (optional) detailed description of the variable
            % - Returns variableAnnotation: a new instance of WVVariableAnnotation
            arguments
                name char {mustBeNonempty}
                dimensions
                units char {mustBeNonempty}
                description char {mustBeNonempty}
                options.isComplex double {mustBeMember(options.isComplex,[0 1])} = 0
                options.detailedDescription char = ''
            end
            self@CANumericProperty(name,dimensions,units,description,isComplex=options.isComplex,detailedDescription=options.detailedDescription);
        end

    end
end
