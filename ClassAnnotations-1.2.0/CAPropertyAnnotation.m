classdef CAPropertyAnnotation < handle & matlab.mixin.Heterogeneous
    % annotates properties
    %
    % The purpose of this class is twofold. First, it lets us add
    % descriptions and detailedDescriptions for any method or property both
    % inline in code, and externally in a markdown file. Second, its
    % subclasses support annotating units and dimensions, which is not only
    % useful for help documentation, but is also necessary for adding those
    % variables to NetCDF files.
    %
    % This class can load a detailedDescription from a .md file with the
    % same name when documentation builders explicitly enable sidecar
    % loading during initialization.
    %
    % - Declaration: classdef CAPropertyAnnotation < handle
    properties (GetAccess=public, SetAccess=private)
        % name of the method, property, or variable
        % 
        % The name should be an exact match to the method, property, or
        % variable that it is describing.
        % - Topic: Properties
        name

        % short description of the method, property, or variable
        % 
        % The short description is used in the table-of-contents of
        % documentation and as the long_name property in NetCDF files.
        % - Topic: Properties
        description
    end

    properties (GetAccess=public, SetAccess=public)
        % a detailed description of the method, property, or variable
        % 
        % The detailed description may be written in markdown and may also
        % use MathJax latex notation. The detailed description will be used
        % for the help documentation.
        %
        % The detailed description can be populated using a markdown
        % sidecar file. Specifically, when documentation builders enable
        % sidecar loading, initialization looks in the call-site directory
        % and its descendants for a markdown file with the name of the
        % annotation and the file extension `.md`.
        % - Topic: Properties
        detailedDescription
    end

    properties (GetAccess=public, SetAccess=public)
        attributes
    end

    methods
        function self = CAPropertyAnnotation(name,description,options)
            % create a new instance of CAPropertyAnnotation
            %
            % Creates a new instance of CAPropertyAnnotation with a name,
            % description and optional detailed description.
            %
            % If documentation builders enable sidecar loading and a
            % markdown file of the same name is in the same directory or a
            % child directory, it will be loaded during initialization.
            %
            % - Topic: Initialization
            % - Declaration: annotation = CAPropertyAnnotation(name,description,options)
            % - Parameter name: name of the method, property, or variable
            % - Parameter description: short description of the method, property, or variable
            % - Parameter detailedDescription: (optional) a detailed description of the method, property, or variable
            % - Returns annotation: a new instance of CAPropertyAnnotation
            arguments
                name char {mustBeNonempty}
                description char {mustBeNonempty}
                options.detailedDescription char = ''
                options.attributes = containers.Map();
            end
            self.name = name;
            self.description = description;
            if isappdata(0, 'CAEnableDetailedDescriptionSidecars_isEnabled') && ...
                    isequal(getappdata(0, 'CAEnableDetailedDescriptionSidecars_isEnabled'), true)
                self.loadDetailedDescriptionIfAvailable();
            end
            self.attributes = options.attributes;
            if ~isempty(options.detailedDescription)
                self.detailedDescription = options.detailedDescription;
            end
        end
    end

    methods (Access=protected)
        function loadDetailedDescriptionIfAvailable(self)
            st = dbstack;
            annotationRoot = fileparts(mfilename('fullpath'));
            searchRoot = '';
            for iStack = 2:length(st)
                callerFile = which(st(iStack).file);
                if isempty(callerFile)
                    callerFile = st(iStack).file;
                end
                if isempty(callerFile)
                    continue
                end
                if startsWith(callerFile, annotationRoot)
                    continue
                end
                searchRoot = fileparts(callerFile);
                break
            end
            if isempty(searchRoot)
                return
            end

            exactPath = fullfile(searchRoot, strcat(self.name, ".md"));
            if isfile(exactPath)
                self.detailedDescription = fileread(exactPath);
                return
            end

            files = dir(fullfile(searchRoot, "**", strcat(self.name, ".md")));
            if ~isempty(files)
                self.detailedDescription = fileread(fullfile(files(1).folder, files(1).name));
            end
        end
    end
end
