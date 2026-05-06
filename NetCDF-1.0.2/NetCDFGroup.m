classdef NetCDFGroup < handle
    % NetCDFGroup A group of NetCDF variables, dimensions, attributes, and child groups.
    %
    % `NetCDFGroup` represents a single NetCDF group (including the root group) and
    % provides a uniform interface for:
    %
    % - Creating and registering coordinate dimensions (and their coordinate variables)
    % - Creating, registering, and reading variables (real and complex)
    % - Creating and traversing a group hierarchy
    % - Managing group-level (global) attributes
    %
    % Typical usage is through `NetCDFFile`, which constructs the root group and then
    % exposes methods that delegate to this class.
    %
    % - Topic: Initializing
    % - Topic: Accessing group properties
    % - Topic: Working with dimensions
    % - Topic: Working with variables
    % - Topic: Working with groups
    % - Topic: Working with global attributes
    %
    % - Declaration: classdef NetCDFGroup < handle

    properties (WeakHandle)
        % parent group (may be nil)
        % - Topic: Accessing group properties
        parentGroup NetCDFGroup
    end

    properties
        % id of the group
        % - Topic: Accessing group properties
        id

        % name of the group
        % - Topic: Accessing group properties
        name


        % array of NetCDFGroup objects
        % - Topic: Working with groups
        groups

        % path of group
        % - Topic: Working with groups
        groupPath = "";

        % array of NetCDFDimension objects
        %
        % An array of NetCDFDimension objects for each coordinate dimension
        % defined in the NetCDF file.
        %
        % Usage
        % ```matlab
        % dim = ncfile.dimensions(dimID+1); % get the dimension with dimID
        % ```
        % - Topic: Working with dimensions
        dimensions

        % array of NetCDFVariable objects
        % - Topic: Working with variables
        realVariables

        % key-value Map of group attributes
        %
        % A `containers.Map` type that contains the key-value pairs of all
        % attributes for this group. This is intended to be
        % *read only*. If you need to add a new attribute to file, use
        % [`addAttribute`](#addattribute).
        %
        % Usage
        % ```matlab
        % model = ncfile.attributes('model');
        % ```
        % - Topic: Working with global attributes
        attributes

        % array of NetCDFComplexVariable objects
        % - Topic: Working with variables
        complexVariables            
    end

    properties (Access=private)
        dimensionIDMap
        dimensionNameMap
        realVariableIDMap
        realVariableNameMap
        complexVariableNameMap
        groupIDMap
        groupNameMap
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = NetCDFGroup(options)
            % NetCDFGroup Create and manage a NetCDF group.
            %
            % - Topic: Initializing
            % - Declaration: grp = NetCDFGroup(options)
            % - Parameter options.name: (char) Group name.
            % - Parameter options.id: (double) NetCDF group ID to load from an existing file.
            % - Parameter options.parentGroup: (NetCDFGroup) Parent group. Use [] for the root group.
            % - Returns grp: (NetCDFGroup) New group instance.
            arguments
                options.name char {mustBeNonempty}
                options.id double {mustBeNonempty}
                options.parentGroup = []
            end
            
            self.dimensions = NetCDFDimension.empty(0,0);
            self.realVariables = NetCDFRealVariable.empty(0,0);
            self.groups = NetCDFGroup.empty(0,0);
            self.attributes = containers.Map;
            self.complexVariables = NetCDFComplexVariable.empty(0,0);

            self.dimensionIDMap = configureDictionary('double','NetCDFDimension');
            self.dimensionNameMap = configureDictionary('string','NetCDFDimension');

            self.realVariableIDMap = configureDictionary('double','NetCDFRealVariable');
            self.realVariableNameMap = configureDictionary('string','NetCDFRealVariable');
            self.complexVariableNameMap = configureDictionary('string','NetCDFComplexVariable');

            self.groupIDMap = configureDictionary('double','NetCDFGroup');
            self.groupNameMap = configureDictionary('string','NetCDFGroup');

            if isfield(options,'id')
                self.initializeGroupFromFile(options.parentGroup,options.id);
            else
                requiredFields = {'name'};
                for iField=1:length(requiredFields)
                    if ~isfield(options,requiredFields{iField})
                        error('You must specify %s to initialize a new group.',requiredFields{iField});
                    end
                end
                self.initGroup(options.parentGroup,options.name);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dimensions
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function dim = dimensionWithName(self,name)  
            % dimensionWithName Retrieve a dimension by name.
            %
            % Usage
            % ```matlab
            % xDim = ncfile.dimensionWithName('x');
            % ```
            %

            % - Topic: Working with dimensions
            % - Declaration: dim = dimensionWithName(name)
            % - Parameter name: (char|string) Dimension name.
            % - Returns dim: (NetCDFDimension) Matching dimension.
            arguments (Input)
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Input, Repeating)
                name char
            end
            
            dim = NetCDFDimension.empty(0,0);
            
            % The logical here is a little complicated.
            % If the user provides an exact path, but the variable isn't
            % there, then we error.
            for iArg=1:length(name)
                namePath = split(name{iArg},"/");
                if length(namePath) > 1
                    grp = self;
                    for iGroup=1:(length(namePath)-1)
                        grp = grp.groupWithName(namePath(iGroup));
                    end
                    dim(end+1) = grp.dimensionWithName(namePath(end));
                else
                    if isKey(self.dimensionNameMap,name{iArg})
                        dim(end+1) = self.dimensionNameMap(name{iArg});
                    else
                        dimPaths = self.dimensionPathWithName(name{iArg});
                        numPaths = length(dimPaths);
                        if numPaths == 0
                            error('Unable to find a dimension with the name %s',name{iArg});
                        elseif numPaths == 1
                            dim(end+1) = self.variableWithName(dimPaths);
                        else
                            error('Found more than one dimension with the name %s',name{iArg});
                        end
                    end
                end
            end
        end

        function dims = dimensionWithID(self,dimids)
            % dimensionWithID Retrieve dimensions by NetCDF dimension ID.
            %
            % - Topic: Working with dimensions
            % - Declaration: dims = dimensionWithID(dimids)
            % - Parameter dimids: (double) NetCDF dimension IDs.
            % - Returns dims: (NetCDFDimension) Matching dimensions (in the same order as dimids).
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                dimids (1,:) double
            end
            dims = self.dimensionIDMap(dimids);
        end

        function [dim,var] = addDimension(self,name,value,options)
            % addDimension Add a new coordinate dimension and associated coordinate variable.
            %
            % ```matlab
            % ncfile.addDimension('x',0:9);
            % ```
            %
            % which will immediately add the dimension and variable data to
            % file, OR you may initialize the dimension and an empty
            % variable, e.g.,
            %
            % ```matlab
            % ncfile.addDimension('t',length=10,type='double');
            % ```
            %
            % which may be useful when initializing an unlimited
            % dimension. In this case you may also want to specify whether
            % or not the variable will be complex valued (it will default
            % to assuming the variable will be real valued).
            %
            % - Topic: Working with dimensions
            % - Declaration: [dim,var] = addDimension(name,value,options)
            % - Parameter name: (char|string) Dimension name.
            % - Parameter value: (:) Coordinate values. Omit or pass [] to create an empty coordinate variable.
            % - Parameter options.isUnlimited: (logical) Mark the dimension as unlimited. Default false.
            % - Parameter options.shouldWriteImmediately: (logical) Write metadata and data to file immediately when possible. Default true.
            % - Returns dim: (NetCDFDimension) Created dimension.
            % - Returns var: (NetCDFVariable) Created coordinate variable.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                value = []
                options.length = 0
                options.type char = []
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
            end
            if isKey(self.dimensionNameMap,name) || isKey(self.realVariableNameMap,name)
                error('A dimension with that name already exists.');
            end
            if ~( (options.length > 0 && ~isempty(options.type)) || ~isempty(value))
                error('You must specify either the value of the data, or the length and type.');
            end

            if ~isempty(value)
                n = length(value);
            else
                n = options.length;
            end
            dim = NetCDFDimension(self,name=name,nPoints=n);
            self.addDimensionPrimitive(dim)
            % TODO: This needs to pass this down to child groups

            if ~isempty(value)
                options.type = class(value);
                options.isComplex = ~isreal(value);
            else
                options.isComplex = false;
            end
            var = self.addVariable(name,{dim.name},value,type=options.type,isComplex=options.isComplex,attributes=options.attributes);           
        end

        function dimPath = dimensionPathWithName(self,dimensionName)
            % dimensionPathWithName Return the full group-qualified path of a dimension name.
            %
            % - Topic: Working with dimensions
            % - Declaration: dimPath = dimensionPathWithName(dimensionName)
            % - Parameter dimensionName: (char|string) Dimension name.
            % - Returns dimPath: (string) Full path to the dimension within the file/group hierarchy.
            arguments
                self NetCDFGroup {mustBeNonempty}
                dimensionName char
            end
            dimPath = string.empty(0,0);
            if isKey(self.dimensionNameMap,dimensionName)
                dimPath = cat(1,dimPath,self.dimensionNameMap(dimensionName).namePath);
            end
            for iGroup=1:length(self.groups)
                dimPath = cat(1,dimPath,self.groups(iGroup).dimensionPathWithName(dimensionName));
            end
        end

        function bool = hasDimensionWithName(self,dimName)
            % hasDimensionWithName Determine whether a dimension exists in this group.
            %
            % ```matlab
            % bool = ncfile.hasDimensionWithName('x');
            % ```
            %
            % - Topic: Working with dimensions
            % - Declaration: tf = hasDimensionWithName(dimName)
            % - Parameter dimName: (char|string) Dimension name.
            % - Returns tf: (logical) True if the dimension exists.
            arguments
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Repeating)
                dimName char
            end
            bool = zeros(length(dimName),1);

            for iArg=1:length(dimName)
                bool(iArg) = any(~isempty(self.dimensionPathWithName(dimName{iArg})));
            end
            bool = logical(bool);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variables
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function var = addVariable(self,name,dimNames,value,options)
            % addVariable Add a real or complex variable to the group.
            %
            % ```matlab
            % ncfile.addVariable('fluid-tracer', {'x','y','t'}, myVariableData);
            % ```
            %
            % or intializes an variable without setting the data,
            % ```matlab
            % ncfile.addVariable('fluid-tracer', {'x','y','t'});
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: var = addVariable(name,dimNames,value,options)
            % - Parameter name: (char|string) Variable name.
            % - Parameter dimNames: (cell|string) Dimension names (coordinate dimensions) defining the variable.
            % - Parameter value: (:) Variable data. Omit or pass [] to create metadata without writing data.
            % - Parameter options.isComplex: (logical) Create a complex variable wrapper. Default false.
            % - Parameter options.shouldWriteImmediately: (logical) Write metadata and data to file immediately when possible. Default true.
            % - Returns var: (NetCDFVariable) Created variable object.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                dimNames = {}
                value {mustBeNumericOrLogical} = [] 
                options.type char = []
                options.isComplex logical
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
            end
            if isKey(self.realVariableNameMap,name) || isKey(self.complexVariableNameMap,name)
                error('A variable with that name already exists.');
            end
            if ~( (~isempty(options.type) && ~isempty(options.isComplex)) || ~isempty(value))
                error('You must specify either the value of the data, or the ncType and isComplex.');
            end

            if ~isempty(value) && isempty(options.type)
                options.type = class(value);
            end
            if ~isempty(value) && ~isfield(options,'isComplex')
                options.isComplex = ~isreal(value);
            end
            
            if ~isempty(value) && isfield(options,'isComplex') && options.isComplex == false && ~isreal(value)
                warning('The variable %s is complex, but the property annotation indicates that it is real.\n',name);
            end

            if options.isComplex==1
                var = NetCDFComplexVariable(group=self,name=name,dimensions=self.dimensionWithName(dimNames{:}),attributes=options.attributes,type=options.type);
                self.addComplexVariablePrimitive(var);
            else
                var = NetCDFRealVariable(self,name=name,dimensions=self.dimensionWithName(dimNames{:}),attributes=options.attributes,type=options.type);
                self.addRealVariablePrimitive(var);
            end
            if ~isempty(value)
                var.value = value;
            end

        end

        function var = addFunctionHandle(self,name,value,options)
            % addFunctionHandle Add a lazily-evaluated variable backed by a function handle.
            %
            % ```matlab
            % m = 1; b = 2;
            % f = @(x) m*x + b;
            % ncfile.addFunctionHandle('f',f);
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: var = addFunctionHandle(name,value,options)
            % - Parameter name: (char|string) Variable name.
            % - Parameter value: (function_handle) Function handle used to generate variable data on demand.
            % - Parameter options.dimNames: (cell|string) Dimension names defining the variable.
            % - Parameter options.isComplex: (logical) Create a complex variable wrapper. Default false.
            % - Returns var: (NetCDFVariable) Created variable object.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                name {mustBeText}
                value function_handle
                options.attributes containers.Map = containers.Map(KeyType='char',ValueType='any');
            end
            if isKey(self.realVariableNameMap,name) || isKey(self.complexVariableNameMap,name)
                error('A variable with that name already exists.');
            end

            [binaryData, warningMessage, warningIdentifier] = NetCDFGroup.serializeFunctionHandle(name, value);
            if strlength(string(warningMessage)) > 0
                warning('NetCDFGroup:PotentiallyUnserializableFunctionHandle', '%s', NetCDFGroup.functionHandleWarningMessage(name, value, warningIdentifier, warningMessage));
            end

            dim = NetCDFDimension(self,name=name,nPoints=length(binaryData));
            self.addDimensionPrimitive(dim)
            options.attributes(NetCDFVariable.GLNetCDFSchemaIsFunctionHandleTypeKey) = uint8(1);
            options.attributes("format") = "The function_handle object is saved in a .mat file, which is then written as binary data to this variable.";

            var = NetCDFRealVariable(self,name=name,dimensions=dim,attributes=options.attributes,type=class(binaryData));
            self.addRealVariablePrimitive(var);
            var.value = binaryData;
        end


        function varargout = readVariables(self,variableNames)
            % readVariables Read one or more variables from file.
            %
            % ```matlab
            % [x,y] = ncfile.readVariables('x','y');
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: varargout = readVariables(variableNames)
            % - Parameter variableNames: (cell|string) Variable names or paths to read.
            % - Returns varargout: (cell) Variable data returned in the same order as variableNames.
            arguments
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Repeating)
                variableNames char
            end
            varargout = cell(size(variableNames));
            for iArg=1:length(variableNames)
                varargout{iArg} = self.variableWithName(variableNames{iArg}).value;
            end
        end

        function varargout = readVariablesAtIndexAlongDimension(self,dimName,index,variableNames)
            % readVariablesAtIndexAlongDimension Read variables at a single index along a named dimension.
            %
            % ```matlab
            % [u,v] = ncfile.readVariablesAtIndexAlongDimension('t',100,'u','v');
            % ```
            %
            % - Topic: Working with variables
            % - Declaration: varargout = readVariablesAtIndexAlongDimension(dimName,index,variableNames)
            % - Parameter dimName: (char|string) Dimension name to index.
            % - Parameter index: (double) 1-based index along dimName.
            % - Parameter variableNames: (cell|string) Variable names or paths to read.
            % - Returns varargout: (cell) Variable data slices returned in the same order as variableNames.
            arguments
                self NetCDFGroup {mustBeNonempty}
                dimName char {mustBeNonempty}
                index  (1,1) double {mustBePositive} = 1
            end
            arguments (Repeating)
                variableNames char
            end
            varargout = cell(size(variableNames));
            for iArg=1:length(variableNames)
                varargout{iArg} = self.variableWithName(variableNames{iArg}).valueAlongDimensionAtIndex(dimName,index);
            end
        end

        function varargout = variableWithName(self,variableName)
            % variableWithName Retrieve a variable by name or path.
            %
            % ```matlab
            % var = ncfile.variableWithName('x');
            % ```
            %
            % var will be either NetCDFRealVariable or
            % NetCDFComplexVariable.
            %
            % - Topic: Working with variables
            % - Declaration: var = variableWithName(,variableName)
            % - Parameter variableName: (char|string) Variable name or group-qualified path.
            % - Returns var: (NetCDFVariable) Matching variable object.
            arguments
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Repeating)
                variableName char
            end
            varargout = cell(size(variableName));
            
            % The logical here is a little complicated.
            % If the user provides an exact path, but the variable isn't
            % there, then we error.
            for iArg=1:length(variableName)
                variablePath = split(variableName{iArg},"/");
                if length(variablePath) > 1
                    grp = self;
                    for iGroup=1:(length(variablePath)-1)
                        grp = grp.groupWithName(variablePath(iGroup));
                    end
                    varargout{iArg} = grp.variableWithName(variablePath(end));
                else
                    if isKey(self.complexVariableNameMap,variableName{iArg})
                        varargout{iArg} = self.complexVariableNameMap(variableName{iArg});
                    elseif isKey(self.realVariableNameMap,variableName{iArg})
                        varargout{iArg} = self.realVariableNameMap(variableName{iArg});
                    else
                        varPaths = self.variablePathsWithName(variableName{iArg});
                        numPaths = length(varPaths);
                        if numPaths == 0
                            error('Unable to find a variable with the name %s',variableName{iArg});
                        elseif numPaths == 1
                            varargout{iArg} = self.variableWithName(varPaths);
                        else
                            error('Found more than one variable with the name %s',variableName{iArg});
                        end
                    end
                end
            end
        end

        function varPaths = variablePathsWithName(self,variableName)
            % variablePathsWithName Return all group-qualified paths matching a variable name.
            %
            % ```matlab
            % varPaths = ncfile.variablePathsWithName('x');
            % ```
            %
            % 
            %
            % - Topic: Working with variables
            % - Declaration: varPaths = variablePathsWithName(grp,variableName)
            % - Parameter variableName: (char|string) Variable name.
            % - Returns varPaths: (string) Matching variable paths.
            arguments
                self NetCDFGroup {mustBeNonempty}
                variableName char
            end
            varPaths = string.empty(0,0);
            if isKey(self.complexVariableNameMap,variableName)
                varPaths = cat(1,varPaths,self.complexVariableNameMap(variableName).namePath);
            elseif isKey(self.realVariableNameMap,variableName)
                varPaths = cat(1,varPaths,self.realVariableNameMap(variableName).namePath);
            end
            for iGroup=1:length(self.groups)
                varPaths = cat(1,varPaths,self.groups(iGroup).variablePathsWithName(variableName));
            end
        end

        function bool = hasVariableWithName(self,variablePath)
            % hasVariableWithName Determine whether a variable exists in this group.
            %
            % ```matlab
            % bool = ncfile.hasVariableWithName('t');
            % ```
            %
            % will return true if there exists one or more variables with
            % that name.
            %
            % ```matlab
            % bool = ncfile.hasVariableWithName('wave-vortex/t');
            % ```
            % 
            % will return true only if exactly that variable path exists.
            %

            % - Topic: Working with variables
            % - Declaration: tf = hasVariableWithName(grp,variablePath)
            % - Parameter variablePath: (char|string) Variable name or group-qualified path.
            % - Returns tf: (logical) True if the variable exists.
            arguments
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Repeating)
                variablePath char
            end
            bool = zeros(length(variablePath),1);

            for iArg=1:length(variablePath)
                variablePathComponents = split(variablePath{iArg},"/");
                if length(variablePathComponents) > 1
                    bool(iArg) = any(self.variablePathsWithName(variablePathComponents{end}) == variablePath{iArg});
                else
                    bool(iArg) = any(~isempty(self.variablePathsWithName(variablePathComponents{end})));
                end
            end
            bool = logical(bool);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Groups
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bool = hasGroupWithName(self,groupName)
            % hasGroupWithName Determine whether a child group exists in this group.
            %
            % - Topic: Working with groups
            % - Declaration: tf = hasGroupWithName(grp,groupName)
            % - Parameter groupName: (char|string) Child group name.
            % - Returns tf: (logical) True if the group exists.
            arguments
                self NetCDFGroup {mustBeNonempty}
            end
            arguments (Repeating)
                groupName char
            end
            bool = zeros(length(groupName),1);

            for iArg=1:length(groupName)
                bool(iArg) = isKey(self.groupNameMap,groupName{iArg});
            end
            bool = logical(bool);
        end

        function grp = groupWithName(self,name)
            % groupWithName Retrieve a child group by name.
            %
            % - Topic: Working with groups
            % - Declaration: grpOut = groupWithName(grp,name)
            % - Parameter name: (char|string) Child group name.
            % - Returns grpOut: (NetCDFGroup) Matching child group.
            if ~isKey(self.groupNameMap,name)
                error('A group with the name %s does not exist!',name);
            end
            grp = self.groupNameMap(name);
        end

        function grp = addGroup(self, name)
            % addGroup Create a new child group.
            %
            % - Topic: Working with groups
            % - Declaration: grpOut = addGroup(grp,name)
            % - Parameter name: (char|string) New child group name.
            % - Returns grpOut: (NetCDFGroup) Created child group.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeText}
            end
            if isKey(self.groupNameMap,name)
                error('A group with that name already exists.');
            end
            grp = NetCDFGroup(parentGroup=self,name=name);
            self.addGroupPrimitive(grp);
        end

        function addDuplicateGroup(self,sourceGroup,options)
            % addDuplicateGroup Duplicate an existing group into this group.
            %
            % - Topic: Working with groups
            % - Declaration: addDuplicateGroup(grp,sourceGroup,options)
            % - Parameter sourceGroup: (NetCDFGroup) Group to duplicate.
            % - Parameter options.copyVariables: (logical) Copy variables. Default true.
            % - Parameter options.copyDimensions: (logical) Copy dimensions. Default true.
            % - Parameter options.copyAttributes: (logical) Copy attributes. Default true.
            arguments
                self NetCDFGroup
                sourceGroup NetCDFGroup
                options.shouldAddToSelf = false
                options.indexRange = dictionary(string.empty,cell.empty)
            end
            if options.shouldAddToSelf
                grp = self;
            else
                grp = self.addGroup(sourceGroup.name);
            end
            for iDim=1:length(sourceGroup.dimensions)
                sourceDim = sourceGroup.dimensions(iDim);
                if ~grp.hasDimensionWithName(sourceDim.name)
                    if sourceDim.isMutable
                        dim = NetCDFDimension(grp,name=sourceDim.name,nPoints=Inf);
                    else
                        dim = NetCDFDimension(grp,name=sourceDim.name,nPoints=sourceDim.nPoints);
                    end
                    dim.isMutable = sourceDim.isMutable;
                    grp.addDimensionPrimitive(dim);
                end
            end

            for iVar=1:length(sourceGroup.realVariables)
                sourceVar = sourceGroup.realVariables(iVar);
                if ~grp.hasVariableWithName(sourceVar.name)
                    if sourceVar.attributes.Count
                        attributesCopy = containers.Map(sourceVar.attributes.keys,sourceVar.attributes.values);
                    else
                        attributesCopy = containers.Map();
                    end
                    dimNames = {sourceVar.dimensions.name};
                    var = NetCDFRealVariable(grp,name=sourceVar.name,dimensions=grp.dimensionWithName(dimNames{:}),attributes=attributesCopy,type=sourceVar.type);
                    grp.addRealVariablePrimitive(var);
                    mutableDimID = find([var.dimensions.isMutable],1,'first');
                    if isempty(mutableDimID)
                        if attributesCopy.isKey(NetCDFVariable.GLNetCDFSchemaIsFunctionHandleTypeKey) && attributesCopy(NetCDFVariable.GLNetCDFSchemaIsFunctionHandleTypeKey)
                            data = netcdf.getVar(sourceGroup.id,sourceVar.id);
                            netcdf.putVar(grp.id, var.id, data);
                        else
                            var.value = sourceVar.value;
                        end
                    else
                        mutableDim = var.dimensions(mutableDimID);
                        if options.indexRange.isKey(sourceDim.namePath)
                            indexRange = options.indexRange{sourceDim.namePath};
                        else
                            indexRange = 1:mutableDim.nPoints;
                        end
                        for index=1:length(indexRange)
                            data = sourceVar.valueAlongDimensionAtIndex(mutableDim.name,indexRange(index));
                            var.setValueAlongDimensionAtIndex(data,mutableDim.name,index)
                        end
                    end
                end
            end

            for iVar=1:length(sourceGroup.complexVariables)
                sourceVar = sourceGroup.complexVariables(iVar);
                if ~grp.hasVariableWithName(sourceVar.name)
                    if sourceVar.attributes.Count
                        attributesCopy = containers.Map(sourceVar.attributes.keys,sourceVar.attributes.values);
                    else
                        attributesCopy = containers.Map();
                    end
                    dimNames = {sourceVar.dimensions.name};
                    var = NetCDFComplexVariable(group=grp,name=sourceVar.name,dimensions=grp.dimensionWithName(dimNames{:}),attributes=attributesCopy,type=sourceVar.type);
                    grp.addComplexVariablePrimitive(var);

                    mutableDimID = find([var.dimensions.isMutable],1,'first');
                    if isempty(mutableDimID)
                        var.value = sourceVar.value;
                    else
                        mutableDim = var.dimensions(mutableDimID);
                        for index=1:mutableDim.nPoints
                            data = sourceVar.valueAlongDimensionAtIndex(mutableDim.name,index);
                            var.setValueAlongDimensionAtIndex(data,mutableDim.name,index)
                        end
                    end
                end
            end

            for attName = string(sourceGroup.attributes.keys)
                grp.addAttribute(attName,sourceGroup.attributes(attName));
            end

            for iGroup=1:length(sourceGroup.groups)
                group = sourceGroup.groups(iGroup);
                grp.addDuplicateGroup(group,indexRange=options.indexRange);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Attributes
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function addAttribute(self,name,data)
            % addAttribute Add or replace a global attribute on this group.
            %
            % - Topic: Working with global attributes
            % - Declaration: addAttribute(grp,name,data)
            % - Parameter name: (char|string) Attribute name.
            % - Parameter data: (:) Attribute value.
            netcdf.putAtt(self.id,netcdf.getConstant('NC_GLOBAL'), name, data);
            self.attributes(name) = data;
        end

        function bool = hasAttributeWithName(self,name)
            % hasAttributeWithName Determine whether a global attribute exists on this group.
            %
            % - Topic: Working with global attributes
            % - Declaration: tf = hasAttributeWithName(grp,name)
            % - Parameter name: (char|string) Attribute name.
            % - Returns tf: (logical) True if the attribute exists.
            bool = isKey(self.attributes,name);
        end

        function dump(self,options)
            % dump Display group contents to the command window.
            %
            % - Topic: Accessing group properties
            % - Declaration: dump(grp,options)
            % - Parameter options.indent: (double) Number of spaces to indent nested output. Default 0.
            % - Parameter options.maxArrayElements: (double) Max elements to display for array-valued attributes. Default 10.
            arguments
                self NetCDFGroup
                options.indentLevel = 0
            end
            indent0 = sprintf(repmat('   ',1,options.indentLevel));
            indent1 = sprintf(repmat('   ',1,options.indentLevel+1));
            indent2 = sprintf(repmat('   ',1,options.indentLevel+2));
            indent3 = sprintf(repmat('   ',1,options.indentLevel+3));

            if isempty(self.parentGroup) && isa(self,'NetCDFFile')
                [~,filename,~] = fileparts(self.path);
                fprintf('%snetcdf: %s {\n',indent0,filename);
            else
                fprintf('\n%sgroup: %s {\n',indent0,self.name);
            end

            if ~isempty(self.dimensions)
                fprintf('%sdimensions: \n',indent1);
                for iDim=1:length(self.dimensions)
                    dimension = self.dimensions(iDim);
                    fprintf('%s%s = %d\n',indent2,dimension.name,dimension.nPoints);
                end
            end

            if ~isempty(self.realVariables)
                fprintf('\n%svariables: \n',indent1);
                for iVar=1:length(self.realVariables)
                    variable = self.realVariables(iVar);
                    fprintf('%s%s %s(',indent2,variable.type,variable.name);
                    for iDim=1:length(variable.dimensions)
                        if iDim==length(variable.dimensions)
                            fprintf('%s',variable.dimensions(iDim).name);
                        else
                            fprintf('%s,',variable.dimensions(iDim).name);
                        end
                    end
                    fprintf(')\n');
                    for attName = string(variable.attributes.keys)
                        if isa(variable.attributes(attName),'char') || isa(variable.attributes(attName),'string')
                            fprintf('%s%s = \"%s\"\n',indent3,attName,variable.attributes(attName));
                        elseif isa(variable.attributes(attName),'double') || isa(variable.attributes(attName),'single')
                            fprintf('%s%s = %f\n',indent3,attName,variable.attributes(attName));
                        else
                            fprintf('%s%s = %d\n',indent3,attName,variable.attributes(attName));
                        end
                    end
                end
            end

            if ~isempty(self.complexVariables)
                % fprintf('\n%scomplex variables: \n',indent1);
                for iVar=1:length(self.complexVariables)
                    variable = self.complexVariables(iVar);
                    fprintf('%scomplex %s %s(',indent2,variable.type,variable.name);
                    for iDim=1:length(variable.dimensions)
                        if iDim==length(variable.dimensions)
                            fprintf('%s',variable.dimensions(iDim).name);
                        else
                            fprintf('%s,',variable.dimensions(iDim).name);
                        end
                    end
                    fprintf(')\n');
                    for attName = string(variable.attributes.keys)
                        if isa(variable.attributes(attName),'char')
                            fprintf('%s%s = \"%s\"\n',indent3,attName,variable.attributes(attName));
                        elseif isa(variable.attributes(attName),'string')
                            fprintf('%s%s = \"%s\"\n',indent3,attName,strjoin(variable.attributes(attName),", "));
                        elseif isa(variable.attributes(attName),'double') || isa(variable.attributes(attName),'single')
                            fprintf('%s%s = %f\n',indent3,attName,variable.attributes(attName));
                        else
                            fprintf('%s%s = %d\n',indent3,attName,variable.attributes(attName));
                        end
                    end
                end
            end

            if self.attributes.Count > 0
                fprintf('\n%sattributes: \n',indent1);
                for attName = string(self.attributes.keys)
                    if isa(self.attributes(attName),'char')
                        fprintf('%s%s = \"%s\"\n',indent2,attName,self.attributes(attName));
                    elseif isa(self.attributes(attName),'string')
                        fprintf('%s%s = \"%s\"\n',indent2,attName,strjoin(self.attributes(attName),", "));
                    elseif isa(self.attributes(attName),'double') || isa(self.attributes(attName),'single')
                        fprintf('%s%s = %f\n',indent2,attName,self.attributes(attName));
                    else
                        fprintf('%s%s = %d\n',indent2,attName,self.attributes(attName));
                    end
                end
            end

            if ~isempty(self.groups)
                for iGroup=1:length(self.groups)
                    group = self.groups(iGroup);
                    group.dump(indentLevel=options.indentLevel+1);
                end
            end

            if isempty(self.parentGroup) && isa(self,'NetCDFFile')
                fprintf('%s}\n',indent0);
            else
                fprintf('%s} // group %s\n',indent0,self.name);
            end
        end
    end

    methods (Static, Hidden)
        function [binaryData, warningMessage, warningIdentifier] = serializeFunctionHandle(variableName, value)
            arguments
                variableName {mustBeText}
                value function_handle
            end
            warningState = warning;
            tmpfile = strcat(tempname, '.mat');
            cleanupWarningState = onCleanup(@() warning(warningState));
            cleanupTmpFile = onCleanup(@() NetCDFGroup.deleteIfPresent(tmpfile));

            warning('off', 'all');
            lastwarn('');
            try
                tmpstruct.(char(variableName)) = value;
                save(tmpfile, "-struct", "tmpstruct");
                roundTrip = load(tmpfile, char(variableName));
                if ~isfield(roundTrip, char(variableName)) || ~isa(roundTrip.(char(variableName)), 'function_handle')
                    error('NetCDFGroup:InvalidFunctionHandleRoundTrip', 'Loading the MAT-file did not reconstruct a function_handle.');
                end
                [warningMessage, warningIdentifier] = lastwarn;

                fileID = fopen(tmpfile, 'rb');
                if fileID == -1
                    error('NetCDFGroup:UnableToReadSerializedFunctionHandle', 'Unable to read the temporary MAT-file for the serialized function handle.');
                end
                cleanupFileID = onCleanup(@() fclose(fileID));
                binaryData = fread(fileID, '*uint8');
            catch ME
                error('NetCDFGroup:UnserializableFunctionHandle', '%s', NetCDFGroup.functionHandleErrorMessage(value, ME.message));
            end
        end

        function message = functionHandleErrorMessage(value, originalError)
            arguments
                value function_handle
                originalError {mustBeText}
            end
            message = strcat("MATLAB could not serialize the function handle ", NetCDFGroup.functionHandleSummary(value), " to a MAT-file. ", NetCDFGroup.functionHandleGuidance(value), " Original error: ", string(originalError));
        end

        function message = functionHandleWarningMessage(name, value, warningIdentifier, warningMessage)
            arguments
                name {mustBeText}
                value function_handle
                warningIdentifier {mustBeText}
                warningMessage {mustBeText}
            end
            message = strcat("MATLAB warned while serializing function handle ", string(name), " = ", NetCDFGroup.functionHandleSummary(value), ". ", NetCDFGroup.warningSummaryText(warningIdentifier, warningMessage), " ", NetCDFGroup.functionHandleGuidance(value));
        end
    end

    methods (Static, Access = private)
        function summary = functionHandleSummary(value)
            arguments
                value function_handle
            end
            functionText = string(func2str(value));
            workspaceNames = NetCDFGroup.capturedWorkspaceVariableNames(value);
            if isempty(workspaceNames)
                summary = strcat("`", functionText, "` with no captured workspace variables detected.");
            else
                summary = strcat("`", functionText, "` capturing workspace variables: ", strjoin(workspaceNames, ", "), ".");
            end
        end

        function guidance = functionHandleGuidance(value)
            arguments
                value function_handle
            end
            info = functions(value);
            guidance = "Prefer named functions or anonymous functions that capture only plain numeric, logical, string, char, struct, or cell data. Keep any referenced functions, package folders, and class definitions available on the MATLAB path when the file is read back. Avoid capturing runtime resources such as onCleanup objects, file IDs, timers, graphics handles, Java/.NET objects, or transient handle objects if the handle must survive a restart.";
            if isfield(info, 'function') && strlength(string(info.function)) > 0
                guidance = strcat(guidance, " Reported function name: ", string(info.function), ".");
            end
        end

        function names = capturedWorkspaceVariableNames(value)
            arguments
                value function_handle
            end
            info = functions(value);
            names = strings(0, 1);
            if ~isfield(info, 'workspace') || isempty(info.workspace)
                return
            end
            workspace = info.workspace;
            if iscell(workspace)
                for iWorkspace = 1:numel(workspace)
                    if isstruct(workspace{iWorkspace})
                        names = union(names, string(fieldnames(workspace{iWorkspace})), 'stable');
                    end
                end
            elseif isstruct(workspace)
                names = union(names, string(fieldnames(workspace)), 'stable');
            end
        end

        function text = warningSummaryText(warningIdentifier, warningMessage)
            arguments
                warningIdentifier {mustBeText}
                warningMessage {mustBeText}
            end
            if strlength(string(warningIdentifier)) > 0
                text = strcat("Warning ", string(warningIdentifier), ": ", string(warningMessage), ".");
            else
                text = strcat("Warning: ", string(warningMessage), ".");
            end
        end

        function deleteIfPresent(path)
            arguments
                path {mustBeText}
            end
            if isfile(path)
                delete(path);
            end
        end
    end

    methods (Access=protected)

        function initializeGroupFromFile(self,parentGroup,id)
            % initializeGroupFromFile Populate this group from an open NetCDF file/group ID.
            %
            % - Topic: Initializing
            % - Declaration: initializeGroupFromFile(grp,parentGroup,id)
            % - Parameter parentGroup: (NetCDFGroup) Parent group (or []).
            % - Parameter id: (double) NetCDF group ID.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                parentGroup
                id (1,1) double {mustBeNonnegative,mustBeInteger}
            end
            self.id = id;
            if ~isempty(parentGroup)
                self.parentGroup = parentGroup;
            end
            self.name = netcdf.inqGrpName(self.id);
            if ~isempty(self.parentGroup)
                if isequal(self.parentGroup.groupPath,"")
                    self.groupPath = self.name;
                else
                    self.groupPath = self.parentGroup.groupPath + "/" + self.name;
                end
            end

            % Fetch all dimensions and convert to NetCDFDimension objects
            if ~isempty(self.parentGroup) && ~isempty(self.parentGroup.dimensions)
                self.addDimensionPrimitive(self.parentGroup.dimensions);
            end
            dimensionIDs = netcdf.inqDimIDs(self.id);
            for iDim=1:length(dimensionIDs)
                self.addDimensionPrimitive(NetCDFDimension(self,id=dimensionIDs(iDim)));
            end

            % Fetch all variables and convert to NetCDFRealVariable objects
            variableIDs = netcdf.inqVarIDs(self.id);
            for iVar=1:length(variableIDs)
                self.addRealVariablePrimitive(NetCDFRealVariable(self,id=variableIDs(iVar)));
            end

            % Grab all the global attributes
            [~,~,ngatts,~] = netcdf.inq(self.id);
            for iAtt=0:(ngatts-1)
                gattname = netcdf.inqAttName(self.id,netcdf.getConstant('NC_GLOBAL'),iAtt);
                self.attributes(gattname) = netcdf.getAtt(self.id,netcdf.getConstant('NC_GLOBAL'),gattname);
            end

            % Now check to see if any variables form a complex variable
            self.addComplexVariablePrimitive(NetCDFComplexVariable.complexVariablesFromVariables(self.realVariables));
            if ~isempty(self.complexVariables)
                self.removeVariablePrimitive([self.complexVariables.realp self.complexVariables.imagp]);
            end

            groupIDs = netcdf.inqGrps(self.id);
            for iGrp=1:length(groupIDs)
                self.addGroupPrimitive(NetCDFGroup(parentGroup=self,id=groupIDs(iGrp)));
            end
        end

        function initGroup(self,parentGroup, name)
            % initGroup Initialize a new group in an open NetCDF file.
            %
            % - Topic: Initializing
            % - Declaration: initGroup(grp,parentGroup,name)
            % - Parameter parentGroup: (NetCDFGroup) Parent group (or []).
            % - Parameter name: (char|string) Group name.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty} 
                parentGroup (1,1) NetCDFGroup {mustBeNonempty} 
                name string {mustBeText}
            end
            self.parentGroup = parentGroup;
            self.name = name;
            self.id = netcdf.defGrp(self.parentGroup.id,self.name);
            if ~isempty(self.parentGroup.dimensions)
                self.addDimensionPrimitive(self.parentGroup.dimensions);
            end
            if ~isempty(self.parentGroup)
                if isequal(self.parentGroup.groupPath,"")
                    self.groupPath = self.name;
                else
                    self.groupPath = self.parentGroup.groupPath + "/" + self.name;
                end
            end
        end

        function addDimensionPrimitive(self,dimension)
            % addDimensionPrimitive Register a NetCDFDimension with this group.
            %
            % - Topic: Working with dimensions
            % - Declaration: addDimensionPrimitive(grp,dimension)
            % - Parameter dimension: (NetCDFDimension) Dimension to register.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                dimension (:,1) NetCDFDimension {mustBeNonempty}
            end
            for iDim=1:length(dimension)
                self.dimensions(end+1) = dimension(iDim);
                self.dimensionIDMap(dimension(iDim).id) = dimension(iDim);
                self.dimensionNameMap(dimension(iDim).name) = dimension(iDim);
            end
            self.dimensions = reshape(self.dimensions,[],1);
        end

        function addRealVariablePrimitive(self,variable)
            % addRealVariablePrimitive Register a NetCDFRealVariable with this group.
            %
            % - Topic: Working with variables
            % - Declaration: addRealVariablePrimitive(grp,variable)
            % - Parameter variable: (NetCDFRealVariable) Variable to register.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFRealVariable {mustBeNonempty}
            end
            for iVar=1:length(variable)
                self.realVariables(end+1) = variable(iVar);
                self.realVariableIDMap(variable(iVar).id) = variable(iVar);
                self.realVariableNameMap(variable(iVar).name) = variable(iVar);
            end
            self.realVariables = reshape(self.realVariables,[],1);
        end

        function removeVariablePrimitive(self,variable)
            % removeVariablePrimitive Remove a variable from this group bookkeeping.
            %
            % - Topic: Working with variables
            % - Declaration: removeVariablePrimitive(grp,variable)
            % - Parameter variable: (NetCDFVariable) Variable to remove.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFRealVariable
            end
            for iVar=1:length(variable)
                self.realVariables(self.realVariables==variable(iVar)) = [];
                self.realVariableIDMap(variable(iVar).id) = [];
                self.realVariableNameMap(variable(iVar).name) = [];
            end
        end

        function addComplexVariablePrimitive(self,variable)
            % addComplexVariablePrimitive Register a NetCDFComplexVariable with this group.
            %
            % - Topic: Working with variables
            % - Declaration: addComplexVariablePrimitive(grp,variable)
            % - Parameter variable: (NetCDFComplexVariable) Variable to register.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                variable (:,1) NetCDFComplexVariable
            end
            for iVar=1:length(variable)
                self.complexVariables(end+1) = variable(iVar);
                self.complexVariableNameMap(variable(iVar).name) = variable(iVar);
            end
            self.complexVariables = reshape(self.complexVariables,[],1);
        end

        function addGroupPrimitive(self,group)
            % addGroupPrimitive Register a child group with this group.
            %
            % - Topic: Working with groups
            % - Declaration: addGroupPrimitive(grp,group)
            % - Parameter group: (NetCDFGroup) Child group to register.
            arguments
                self (1,1) NetCDFGroup {mustBeNonempty}
                group (1,1) NetCDFGroup {mustBeNonempty}
            end
            self.groups(end+1) = group;
            self.groups = reshape(self.groups,[],1);
            self.groupIDMap(group.id) = group;
            self.groupNameMap(group.name) = group;
        end
    end
end
