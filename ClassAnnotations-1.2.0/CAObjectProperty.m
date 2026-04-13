classdef CAObjectProperty < CAPropertyAnnotation
    % Describe an annotated persisted object property.
    %
    % `CAObjectProperty` documents a property whose value is another
    % `CAAnnotatedClass` object or object array. The optional `className`
    % and `sizeText` fields are documentation-only metadata used by
    % `class-docs` when reflected MATLAB property validation is absent or
    % incomplete.
    %
    % Like other property annotations, this class can populate its detailed
    % description from a same-named markdown sidecar when documentation
    % builders enable sidecar loading.
    %
    % - Declaration: classdef CAObjectProperty < CAPropertyAnnotation

    properties (GetAccess=public, SetAccess=public)
        % Documented class name for the object property.
        %
        % This field is documentation-only metadata and does not affect
        % annotated persistence.
        % - Topic: Properties
        className

        % Human-readable size text for the object property.
        %
        % Use this for freeform shapes such as `nonempty vector` when the
        % property does not have reflected MATLAB size validation.
        % - Topic: Properties
        sizeText
    end

    methods
        function self = CAObjectProperty(name,description,options)
            % Create a new object-property annotation.
            %
            % If a markdown file of the same name is in the same directory
            % or child directory, it will be loaded as the detailed
            % description upon initialization.
            %
            % - Topic: Initialization
            % - Declaration: propAnnotation = CAObjectProperty(name,description,options)
            % - Parameter name: name of the property
            % - Parameter description: short description of the property
            % - Parameter detailedDescription: (optional) detailed description of the property
            % - Parameter className: (optional) documentation-only object class name
            % - Parameter sizeText: (optional) documentation-only size text
            % - Returns propAnnotation: a new instance of CAObjectProperty
            arguments
                name char {mustBeNonempty}
                description char {mustBeNonempty}
                options.detailedDescription char = ''
                options.className = ''
                options.sizeText = ''
            end

            self@CAPropertyAnnotation(name,description,detailedDescription=options.detailedDescription);
            self.className = string(options.className);
            self.sizeText = string(options.sizeText);
        end
    end
end
