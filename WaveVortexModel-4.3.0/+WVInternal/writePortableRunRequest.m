function writePortableRunRequest(requestPath,modelFiles,configuration,failureStage)
% Validate, encode, and transactionally install one portable run request.
arguments
    requestPath (1,1) string
    modelFiles string
    configuration (1,1) struct
    failureStage (1,1) string
end

[document,validated] = validateRequest(requestPath,modelFiles,configuration);
injectFailure(failureStage,"validation");
encoded = encodeDocument(document,validated.destinationIdentifiers);
injectFailure(failureStage,"encoding");

requestPath = validated.requestPath;
temporaryPath = string(tempname(fileparts(requestPath))) + ".json";
cleanup = onCleanup(@()deleteIfPresent(temporaryPath));
injectFailure(failureStage,"temporary-file");

bytes = unicode2native(char(encoded),"UTF-8");
fileIdentifier = fopen(temporaryPath,"wb");
if fileIdentifier < 0
    error("WaveVortexModel:PortableRunRequestTemporaryFile","Unable to create the temporary request file %s.",temporaryPath);
end
try
    written = fwrite(fileIdentifier,bytes,"uint8");
    closeStatus = fclose(fileIdentifier);
    fileIdentifier = -1;
catch exception
    if fileIdentifier >= 0
        fclose(fileIdentifier);
    end
    rethrow(exception)
end
if written ~= numel(bytes) || closeStatus ~= 0
    error("WaveVortexModel:PortableRunRequestTemporaryFile","Unable to flush and close the temporary request file %s.",temporaryPath);
end

observedBytes = readBytes(temporaryPath);
if ~isequal(observedBytes(:),uint8(bytes(:)))
    error("WaveVortexModel:PortableRunRequestValidation","The temporary request changed while it was being validated.");
end
validateEncodedDocument(native2unicode(observedBytes',"UTF-8"),document);

if failureStage == "cleanup"
    deleteIfPresent(temporaryPath);
    injectFailure(failureStage,"cleanup");
end
injectFailure(failureStage,"replacement");
atomicReplace(temporaryPath,requestPath);
clear cleanup
end

function [document,validated] = validateRequest(requestPath,modelFiles,configuration)
requiredFields = ["schemaVersion","method","finalTime","initialStep","cfl", ...
    "timeStepConstraint","maximumStep","relativeTolerance", ...
    "absoluteToleranceScale","outputPolicy","destinations", ...
    "fftProvider","threads","reportPath"];
observedFields = string(fieldnames(configuration));
if ~isequal(sort(observedFields),sort(requiredFields(:)))
    error("WaveVortexModel:PortableRunRequestOptions","The portable run-request configuration contains unknown or missing options.");
end
if ~isvector(modelFiles)
    error("WaveVortexModel:PortableRunRequestModelFiles","modelFiles must be an ordered vector of paths.");
end
modelFiles = reshape(modelFiles,1,[]);
if any(ismissing(modelFiles) | strlength(modelFiles) == 0)
    error("WaveVortexModel:PortableRunRequestModelFiles","Every modelFiles entry must be a nonempty path.");
end

requestPath = canonicalPath(requestPath);
requestFolder = string(fileparts(requestPath));
if ~isfolder(requestFolder)
    error("WaveVortexModel:PortableRunRequestPath","The request destination folder does not exist.");
end
resolvedModelFiles = strings(size(modelFiles));
metadata = cell(size(modelFiles));
for iFile = 1:numel(modelFiles)
    resolvedModelFiles(iFile) = resolveDocumentPath(requestFolder,modelFiles(iFile));
    if ~isfile(resolvedModelFiles(iFile))
        error("WaveVortexModel:PortableRunRequestMissingFile","Referenced model file does not exist: %s",modelFiles(iFile));
    end
    if any(resolvedModelFiles(1:iFile-1) == resolvedModelFiles(iFile))
        error("WaveVortexModel:PortableRunRequestDuplicateFile","modelFiles entries must not alias each other.");
    end
    metadata{iFile} = inspectModelFile(resolvedModelFiles(iFile));
end
referenceSignature = metadata{1}.bundleSignature;
for iFile = 2:numel(metadata)
    if ~isequaln(metadata{iFile}.bundleSignature,referenceSignature)
        error("WaveVortexModel:PortableRunRequestInconsistentBundle","The source files do not describe one compatible model bundle.");
    end
end
fileIdentifiers = string(cellfun(@(entry)entry.fileIdentifier,metadata,UniformOutput=false));
if numel(unique(fileIdentifiers)) ~= numel(fileIdentifiers)
    error("WaveVortexModel:PortableRunRequestDuplicateIdentifier","Source files must have unique stable output-file identifiers.");
end

schemaVersion = configuration.schemaVersion;
method = string(configuration.method);
integration = validateIntegration(schemaVersion,method,configuration);
destinations = validateDestinations(configuration.destinations,string(configuration.outputPolicy), ...
    fileIdentifiers,metadata,requestFolder,requestPath,resolvedModelFiles);

fftProvider = string(configuration.fftProvider);
threads = double(configuration.threads);
if schemaVersion == 1
    if strlength(fftProvider) == 0
        fftProvider = "reference";
    end
    if isnan(threads)
        threads = 1;
    end
end
if ~isscalar(threads) || (~isnan(threads) && (~isfinite(threads) || threads <= 0 || fix(threads) ~= threads))
    error("WaveVortexModel:PortableRunRequestExecution","threads must be a positive integer.");
end
if ~ismember(fftProvider,["","native-fftw","reference"])
    error("WaveVortexModel:PortableRunRequestExecution","fftProvider must be native-fftw or reference.");
end
if fftProvider == "reference" && ~isnan(threads) && threads ~= 1
    error("WaveVortexModel:PortableRunRequestExecution","The reference FFT provider requires one thread.");
end

reportPath = string(configuration.reportPath);
if ismissing(reportPath)
    error("WaveVortexModel:PortableRunRequestReport","reportPath must not be missing.");
end
if strlength(reportPath) == 0
    [~,requestName] = fileparts(requestPath);
    reportPath = requestName + "-report.json";
end
resolvedReport = resolveDocumentPath(requestFolder,reportPath);
reservedPaths = [requestPath resolvedModelFiles destinations.resolvedPaths];
if any(resolvedReport == reservedPaths)
    error("WaveVortexModel:PortableRunRequestAlias","The report path must not alias the request, a model file, or an output destination.");
end

document = struct;
document.schemaIdentifier = "wave-vortex-run-request-v" + string(schemaVersion);
document.schemaVersion = schemaVersion;
document.modelFiles = modelFiles;
document.integration = integration;
document.outputPolicy = string(configuration.outputPolicy);
document.destinations = destinations.map;
document.fftProvider = fftProvider;
document.threads = threads;
document.report = reportPath;
validated = struct("requestPath",requestPath,"destinationIdentifiers",destinations.identifiers);
end

function integration = validateIntegration(schemaVersion,method,configuration)
if ~isscalar(schemaVersion) || ~ismember(schemaVersion,[1 2]) || fix(schemaVersion) ~= schemaVersion
    error("WaveVortexModel:PortableRunRequestSchema","schemaVersion must be exactly 1 or 2.");
end
if ~isscalar(configuration.finalTime) || ~isfinite(configuration.finalTime)
    error("WaveVortexModel:PortableRunRequestIntegration","finalTime must be finite.");
end
allowedV1 = ["fixed-rk4","adaptive-rk23"];
allowedV2 = [allowedV1,"adaptive-rk45","adaptive-rk78"];
if schemaVersion == 1 && ~ismember(method,allowedV1)
    error("WaveVortexModel:PortableRunRequestIntegration","The integration method is unavailable in schema version %d.",schemaVersion);
end
if schemaVersion == 2 && ~ismember(method,["" allowedV2])
    error("WaveVortexModel:PortableRunRequestIntegration","The integration method is unavailable in schema version %d.",schemaVersion);
end

initialStep = optionalPositive(configuration.initialStep,"initialStep");
cfl = optionalPositive(configuration.cfl,"cfl");
maximumStep = optionalPositive(configuration.maximumStep,"maximumStep");
relativeTolerance = optionalPositive(configuration.relativeTolerance,"relativeTolerance");
absoluteToleranceScale = optionalPositive(configuration.absoluteToleranceScale,"absoluteToleranceScale");
constraint = string(configuration.timeStepConstraint);
hasConstraint = strlength(constraint) > 0;

integration = struct("finalTime",double(configuration.finalTime));
if strlength(method) > 0
    integration = struct;
    integration.method = method;
    integration.finalTime = double(configuration.finalTime);
end
effectiveMethod = method;
if strlength(effectiveMethod) == 0
    effectiveMethod = "adaptive-rk78";
end
if effectiveMethod == "fixed-rk4"
    if ~isnan(maximumStep) || ~isnan(relativeTolerance) || ~isnan(absoluteToleranceScale)
        error("WaveVortexModel:PortableRunRequestIntegration","Adaptive tolerances and maximumStep are not valid for fixed-rk4.");
    end
    if schemaVersion == 1
        if isnan(initialStep) || ~isnan(cfl) || hasConstraint
            error("WaveVortexModel:PortableRunRequestIntegration","Schema-v1 fixed-rk4 requires initialStep and rejects CFL controls.");
        end
        integration.initialStep = initialStep;
        return
    end
    if isnan(initialStep) == isnan(cfl)
        error("WaveVortexModel:PortableRunRequestIntegration","Schema-v2 fixed-rk4 requires exactly one of initialStep or cfl.");
    end
    if ~isnan(initialStep)
        if hasConstraint
            error("WaveVortexModel:PortableRunRequestIntegration","timeStepConstraint requires cfl.");
        end
        integration.initialStep = initialStep;
        return
    end
    if ~ismember(constraint,["advective","oscillatory","min"])
        error("WaveVortexModel:PortableRunRequestIntegration","CFL-selected fixed-rk4 requires timeStepConstraint=advective, oscillatory, or min.");
    end
    integration.cfl = cfl;
    integration.timeStepConstraint = constraint;
    return
end

if ~isnan(cfl) || hasConstraint
    error("WaveVortexModel:PortableRunRequestIntegration","CFL controls are valid only for fixed-rk4.");
end
if schemaVersion == 1 && any(isnan([initialStep maximumStep relativeTolerance absoluteToleranceScale]))
    error("WaveVortexModel:PortableRunRequestIntegration","Adaptive methods require initialStep, maximumStep, relativeTolerance, and absoluteToleranceScale.");
end
if ~isnan(initialStep)
    integration.initialStep = initialStep;
end
if ~isnan(maximumStep)
    integration.maximumStep = maximumStep;
end
if ~isnan(relativeTolerance)
    integration.relativeTolerance = relativeTolerance;
end
if ~isnan(absoluteToleranceScale)
    integration.absoluteToleranceScale = absoluteToleranceScale;
end
end

function value = optionalPositive(value,name)
if ~isscalar(value) || ~isreal(value) || (~isnan(value) && (~isfinite(value) || value <= 0))
    error("WaveVortexModel:PortableRunRequestIntegration","%s must be finite and positive when supplied.",name);
end
value = double(value);
end

function result = validateDestinations(destinationMap,policy,fileIdentifiers,sourceMetadata,requestFolder,requestPath,resolvedModelFiles)
if ~isa(destinationMap,"dictionary")
    error("WaveVortexModel:PortableRunRequestDestinations","destinations must be a string-to-string dictionary.");
end
if ~ismember(policy,["create","replace","append"])
    error("WaveVortexModel:PortableRunRequestDestinations","outputPolicy must be create, replace, or append.");
end
identifiers = string(destinationMap.keys);
identifiers = reshape(identifiers,1,[]);
if (policy == "create" || policy == "replace" || ~isempty(identifiers)) && ...
        ~isequal(sort(identifiers),sort(fileIdentifiers))
    error("WaveVortexModel:PortableRunRequestDestinations", ...
        "Output destination remapping requires exactly one path for every stable file identifier: %s", ...
        strjoin(sort(fileIdentifiers),", "));
end
values = strings(size(identifiers));
resolved = strings(size(identifiers));
for iDestination = 1:numel(identifiers)
    values(iDestination) = string(destinationMap(identifiers(iDestination)));
    if ismissing(values(iDestination)) || strlength(values(iDestination)) == 0
        error("WaveVortexModel:PortableRunRequestDestinations","Every output destination must be a nonempty path.");
    end
    resolved(iDestination) = resolveDocumentPath(requestFolder,values(iDestination));
    if any(resolved(1:iDestination-1) == resolved(iDestination))
        error("WaveVortexModel:PortableRunRequestAlias","Output destinations must not alias each other.");
    end
    if resolved(iDestination) == requestPath
        error("WaveVortexModel:PortableRunRequestAlias","An output destination must not alias the request document.");
    end
    if policy ~= "append" && any(resolved(iDestination) == resolvedModelFiles)
        error("WaveVortexModel:PortableRunRequestAlias","Create and replace destinations must not alias a source model file.");
    end
    if policy == "create" && isfile(resolved(iDestination))
        error("WaveVortexModel:PortableRunRequestDestinations","Create policy requires every output destination to be absent.");
    elseif policy ~= "create" && ~isfile(resolved(iDestination))
        error("WaveVortexModel:PortableRunRequestDestinations","Replace and append policies require every mapped output destination to exist.");
    end
    if policy == "append"
        sourceIndex = find(fileIdentifiers == identifiers(iDestination),1);
        destinationMetadata = inspectModelFile(resolved(iDestination));
        if destinationMetadata.fileIdentifier ~= identifiers(iDestination) || ...
                ~isequaln(destinationMetadata.bundleSignature,sourceMetadata{sourceIndex}.bundleSignature)
            error("WaveVortexModel:PortableRunRequestInconsistentBundle", ...
                "An append destination is incompatible with its source output-file identifier.");
        end
    end
end
[identifiers,order] = sort(identifiers);
values = values(order);
if isempty(identifiers)
    sortedMap = configureDictionary("string","string");
else
    sortedMap = configureDictionary("string","string");
    sortedMap(identifiers) = values;
end
result = struct("map",sortedMap, ...
    "identifiers",identifiers,"resolvedPaths",resolved(order));
end

function metadata = inspectModelFile(path)
try
    information = ncinfo(path);
catch exception
    throwAsCaller(addCause(MException("WaveVortexModel:PortableRunRequestNetCDF", ...
        "Unable to inspect NetCDF metadata for %s.",path),exception));
end
annotatedClass = textAttribute(information.Attributes,"AnnotatedClass","");
wvTransform = textAttribute(information.Attributes,"WVTransform","");
if strlength(annotatedClass) == 0 && strlength(wvTransform) == 0
    error("WaveVortexModel:PortableRunRequestContract","The NetCDF root has no transform identity.");
end
if strlength(annotatedClass) > 0 && strlength(wvTransform) > 0 && annotatedClass ~= wvTransform
    error("WaveVortexModel:PortableRunRequestContract","WVTransform and AnnotatedClass root metadata disagree.");
end
transformClass = annotatedClass;
if strlength(transformClass) == 0
    transformClass = wvTransform;
end
modelVersion = textAttribute(information.Attributes,"model_version","");
majorVersion = regexp(modelVersion,"^4(?:\.|$)","once");
if isempty(majorVersion)
    error("WaveVortexModel:PortableRunRequestContract","The structural portable profile requires a WaveVortexModel 4.x file.");
end

dynamicsMode = numericAttribute(information.Attributes,"WVModelIsDynamicsLinear",0);
switch transformClass
    case "WVTransformConstantStratification"
        bundleSignature = inspectConstantStratificationMetadata(path,information,transformClass,modelVersion,dynamicsMode);
    case "WVTransformBarotropicQG"
        bundleSignature = inspectBarotropicQGMetadata(path,information,transformClass,modelVersion,dynamicsMode);
    otherwise
        error("WaveVortexModel:PortableRunRequestContract", ...
            "The portable runtime does not support transform %s.",transformClass);
end

fileIdentifier = textAttribute(information.Attributes,"portableFileIdentifier","");
if strlength(fileIdentifier) == 0
    fileIdentifier = regexprep(path,"[^A-Za-z0-9_.-]","-");
end
if strlength(fileIdentifier) == 0
    fileIdentifier = "unnamed";
end
metadata = struct("fileIdentifier",fileIdentifier,"bundleSignature",bundleSignature);
end

function signature = inspectConstantStratificationMetadata(path,information,transformClass,modelVersion,dynamicsMode)
supportedGroups = false(size(information.Groups));
completeRestartGroups = 0;
for iGroup = 1:numel(information.Groups)
    group = information.Groups(iGroup);
    if textAttribute(group.Attributes,"AnnotatedClass","") ~= "WVModelOutputGroupEvenlySpaced"
        continue
    end
    supportedGroups(iGroup) = true;
    validatePortableGroupContracts(group);
    variableNames = string({group.Variables.Name});
    completeFamilies = true;
    for family = ["Ap","Am","A0"]
        hasPlain = any(variableNames == family);
        hasReal = any(variableNames == family + "_real");
        hasImaginary = any(variableNames == family + "_imag");
        if hasPlain && (hasReal || hasImaginary) || xor(hasReal,hasImaginary)
            error("WaveVortexModel:PortableRunRequestContract","A coefficient family has ambiguous or incomplete storage in %s.",group.Name);
        end
        completeFamilies = completeFamilies && (hasPlain || (hasReal && hasImaginary));
    end
    completeRestartGroups = completeRestartGroups + completeFamilies;
end
if ~any(supportedGroups) || completeRestartGroups ~= 1
    error("WaveVortexModel:PortableRunRequestContract","Each source file must contain supported output and exactly one complete coefficient stream.");
end

dimensionNames = ["x","y","z"];
dimensionLengths = zeros(1,3);
for iDimension = 1:3
    index = find(string({information.Dimensions.Name}) == dimensionNames(iDimension),1);
    if isempty(index)
        error("WaveVortexModel:PortableRunRequestContract","The NetCDF root is missing dimension %s.",dimensionNames(iDimension));
    end
    dimensionLengths(iDimension) = information.Dimensions(index).Length;
end
scalarNames = ["Lx","Ly","Lz","N0","g","rho0","planetaryRadius","rotationRate","latitude","isHydrostatic","shouldAntialias"];
scalarValues = zeros(size(scalarNames));
for iScalar = 1:numel(scalarNames)
    variableIndex = find(string({information.Variables.Name}) == scalarNames(iScalar),1);
    if isempty(variableIndex) || prod(double(information.Variables(variableIndex).Size)) ~= 1
        error("WaveVortexModel:PortableRunRequestContract","The NetCDF root is missing scalar configuration variable %s.",scalarNames(iScalar));
    end
    scalarValues(iScalar) = double(ncread(path,scalarNames(iScalar)));
    if ~isfinite(scalarValues(iScalar))
        error("WaveVortexModel:PortableRunRequestContract","Configuration variable %s must be finite.",scalarNames(iScalar));
    end
end
signature = struct("transformClass",transformClass,"modelVersion",modelVersion, ...
    "dimensionLengths",dimensionLengths,"scalarValues",scalarValues,"dynamicsMode",dynamicsMode);
end

function signature = inspectBarotropicQGMetadata(path,information,transformClass,modelVersion,dynamicsMode)
dimensionNames = ["x","y"];
dimensionLengths = zeros(1,2);
coordinateValues = cell(1,2);
for iDimension = 1:2
    dimensionName = dimensionNames(iDimension);
    dimensionLengths(iDimension) = rootDimensionLength(information,dimensionName);
    if dimensionLengths(iDimension) < 2
        error("WaveVortexModel:PortableRunRequestContract", ...
            "Barotropic QG dimension %s must contain at least two points.",dimensionName);
    end
    coordinate = rootVariable(information,dimensionName);
    coordinateDimensions = arrayfun(@dimensionLeafName,coordinate.Dimensions);
    if string(coordinate.Datatype) ~= "double" || ...
            ~isequal(coordinateDimensions,dimensionName) || ...
            prod(double(coordinate.Size)) ~= dimensionLengths(iDimension)
        error("WaveVortexModel:PortableRunRequestContract", ...
            "Barotropic QG coordinate %s must be a double vector on its matching root dimension.",dimensionName);
    end
    values = reshape(double(ncread(path,dimensionName)),1,[]);
    if any(~isfinite(values)) || any(diff(values) <= 0)
        error("WaveVortexModel:PortableRunRequestContract", ...
            "Barotropic QG coordinate %s must be finite and strictly increasing.",dimensionName);
    end
    coordinateValues{iDimension} = values;
end

scalarNames = ["Lx","Ly","h","j","g","planetaryRadius","rotationRate","latitude","shouldAntialias"];
scalarValues = zeros(size(scalarNames));
for iScalar = 1:numel(scalarNames)
    scalarValues(iScalar) = readRootScalar(path,information,scalarNames(iScalar));
end
positiveScalars = scalarValues([1 2 3 5 6 7]);
latitude = scalarValues(8);
if any(positiveScalars <= 0) || abs(latitude) < 5 || abs(latitude) > 85
    error("WaveVortexModel:PortableRunRequestContract", ...
        "Barotropic QG lengths and physical constants must be positive, and absolute latitude must lie in [5,85].");
end
if ~ismember(scalarValues(4),[0 1])
    error("WaveVortexModel:PortableRunRequestContract","Barotropic QG mode j must be exactly 0 or 1.");
end
if ~ismember(scalarValues(9),[0 1])
    error("WaveVortexModel:PortableRunRequestContract","Barotropic QG shouldAntialias must be logical scalar metadata.");
end
t0 = readRootScalar(path,information,"t0");

compactLength = rootDimensionLength(information,"kl");
if compactLength < 1
    error("WaveVortexModel:PortableRunRequestContract","Barotropic QG compact kl storage must be nonempty.");
end
supportedGroups = false(size(information.Groups));
completeRestartGroups = 0;
for iGroup = 1:numel(information.Groups)
    group = information.Groups(iGroup);
    if textAttribute(group.Attributes,"AnnotatedClass","") ~= "WVModelOutputGroupEvenlySpaced"
        continue
    end
    supportedGroups(iGroup) = true;
    validatePortableGroupContracts(group);
    completeRestartGroups = completeRestartGroups + validateBarotropicQGRestartGroup(group,compactLength);
end
if ~any(supportedGroups) || completeRestartGroups ~= 1
    error("WaveVortexModel:PortableRunRequestContract", ...
        "Each Barotropic QG source file must contain supported output and exactly one declared compact WVCoefficients A0 stream.");
end

signature = struct("transformClass",transformClass,"modelVersion",modelVersion, ...
    "dimensionLengths",dimensionLengths,"coordinateValues",{coordinateValues}, ...
    "scalarValues",scalarValues,"t0",t0,"compactLength",compactLength,"dynamicsMode",dynamicsMode);
end

function isComplete = validateBarotropicQGRestartGroup(group,compactLength)
isComplete = false;
variableNames = string({group.Variables.Name});
hasPlain = any(variableNames == "A0");
hasReal = any(variableNames == "A0_real");
hasImaginary = any(variableNames == "A0_imag");
if ~hasPlain && ~hasReal && ~hasImaginary
    return
end
coefficientObservers = declaredObserverCount(group,"WVCoefficients");
if coefficientObservers == 0
    return
elseif coefficientObservers ~= 1
    error("WaveVortexModel:PortableRunRequestContract", ...
        "A Barotropic QG output group declares an ambiguous WVCoefficients observer contract in %s.",group.Name);
end
if hasPlain && (hasReal || hasImaginary) || xor(hasReal,hasImaginary)
    error("WaveVortexModel:PortableRunRequestContract", ...
        "Barotropic QG A0 has ambiguous or incomplete compact storage in %s.",group.Name);
end
for forbidden = ["Ap","Ap_real","Ap_imag","Am","Am_real","Am_imag"]
    if any(variableNames == forbidden)
        error("WaveVortexModel:PortableRunRequestContract", ...
            "Barotropic QG compact state must not contain dummy %s storage in %s.",forbidden,group.Name);
    end
end
if hasPlain
    validateCompactVariable(group,"A0",compactLength,false);
else
    validateCompactVariable(group,"A0_real",compactLength,true);
    validateCompactVariable(group,"A0_imag",compactLength,true);
    realVariable = group.Variables(find(variableNames == "A0_real",1));
    imaginaryVariable = group.Variables(find(variableNames == "A0_imag",1));
    if ~isequal(realVariable.Size,imaginaryVariable.Size) || ...
            ~isequal(arrayfun(@dimensionLeafName,realVariable.Dimensions), ...
            arrayfun(@dimensionLeafName,imaginaryVariable.Dimensions))
        error("WaveVortexModel:PortableRunRequestContract", ...
            "Barotropic QG compact A0 components must have identical dimensions in %s.",group.Name);
    end
    validateComplexMarker(realVariable,"isRealPart",1,group.Name);
    validateComplexMarker(realVariable,"isImaginaryPart",0,group.Name);
    validateComplexMarker(imaginaryVariable,"isRealPart",0,group.Name);
    validateComplexMarker(imaginaryVariable,"isImaginaryPart",1,group.Name);
end
isComplete = true;
end

function validateCompactVariable(group,name,compactLength,allowTimeSeries)
variableNames = string({group.Variables.Name});
variable = group.Variables(find(variableNames == name,1));
dimensionNames = arrayfun(@dimensionLeafName,variable.Dimensions);
if allowTimeSeries
    validDimensions = isequal(dimensionNames,"kl") || ...
        (isequal(dimensionNames,["kl" "t"]) && variable.Size(2) > 0);
else
    validDimensions = isequal(dimensionNames,"kl");
end
if string(variable.Datatype) ~= "double" || ~validDimensions || ...
        isempty(variable.Size) || variable.Size(1) ~= compactLength
    error("WaveVortexModel:PortableRunRequestContract", ...
        "Barotropic QG compact variable %s has an incompatible type or shape in %s.",name,group.Name);
end
if allowTimeSeries
    validateComplexMarker(variable,"isComplex",1,group.Name);
end
end

function validateComplexMarker(variable,name,expectedValue,groupName)
value = numericAttribute(variable.Attributes,name,NaN);
if ~isequal(value,double(expectedValue))
    error("WaveVortexModel:PortableRunRequestContract", ...
        "Barotropic QG compact variable %s has an invalid %s marker in %s.",variable.Name,name,groupName);
end
end

function count = declaredObserverCount(group,className)
count = 0;
for iGroup = 1:numel(group.Groups)
    child = group.Groups(iGroup);
    if dimensionLeafName(child) == "observingSystems"
        count = countAnnotatedClass(child,className);
        return
    end
end
end

function count = countAnnotatedClass(group,className)
count = double(textAttribute(group.Attributes,"AnnotatedClass","") == className);
for iGroup = 1:numel(group.Groups)
    count = count + countAnnotatedClass(group.Groups(iGroup),className);
end
end

function lengthValue = rootDimensionLength(information,name)
index = find(string({information.Dimensions.Name}) == name,1);
if isempty(index)
    error("WaveVortexModel:PortableRunRequestContract","The NetCDF root is missing dimension %s.",name);
end
lengthValue = double(information.Dimensions(index).Length);
end

function variable = rootVariable(information,name)
index = find(string({information.Variables.Name}) == name,1);
if isempty(index)
    error("WaveVortexModel:PortableRunRequestContract","The NetCDF root is missing variable %s.",name);
end
variable = information.Variables(index);
end

function value = readRootScalar(path,information,name)
variable = rootVariable(information,name);
if prod(double(variable.Size)) ~= 1
    error("WaveVortexModel:PortableRunRequestContract","The NetCDF root is missing scalar configuration variable %s.",name);
end
value = double(ncread(path,name));
if ~isscalar(value) || ~isfinite(value)
    error("WaveVortexModel:PortableRunRequestContract","Configuration variable %s must be a finite scalar.",name);
end
end

function name = dimensionLeafName(metadata)
name = string(metadata.Name);
parts = split(name,"/");
name = parts(end);
end

function validatePortableGroupContracts(group)
scheduleType = textAttribute(group.Attributes,"portableScheduleTypeIdentifier","");
if strlength(scheduleType) > 0
    if numericAttribute(group.Attributes,"portableScheduleContractVersion",NaN) ~= 1 || ...
            strlength(textAttribute(group.Attributes,"portableScheduleConfiguration","")) == 0
        error("WaveVortexModel:PortableRunRequestContract","An output group declares an unsupported portable schedule contract.");
    end
end
for attributeName = ["portableObservationSchemaVersion","portableObserverContractVersion"]
    version = numericAttribute(group.Attributes,attributeName,NaN);
    if ~isnan(version) && version ~= 1
        error("WaveVortexModel:PortableRunRequestContract","An output group declares an unsupported portable observer contract.");
    end
end
end

function value = textAttribute(attributes,name,defaultValue)
if isempty(attributes)
    value = string(defaultValue);
    return
end
index = find(string({attributes.Name}) == name,1);
if isempty(index)
    value = string(defaultValue);
else
    value = string(attributes(index).Value);
    if ~isscalar(value)
        error("WaveVortexModel:PortableRunRequestContract","Attribute %s must be scalar text.",name);
    end
end
end

function value = numericAttribute(attributes,name,defaultValue)
if isempty(attributes)
    value = defaultValue;
    return
end
index = find(string({attributes.Name}) == name,1);
if isempty(index)
    value = defaultValue;
else
    value = double(attributes(index).Value);
    if ~isscalar(value)
        error("WaveVortexModel:PortableRunRequestContract","Attribute %s must be scalar numeric metadata.",name);
    end
end
end

function encoded = encodeDocument(document,destinationIdentifiers)
lines = strings(0,1);
append("{");
append("  " + jsonString("schemaIdentifier") + ": " + jsonString(document.schemaIdentifier) + ",");
append("  " + jsonString("schemaVersion") + ": " + integerText(document.schemaVersion) + ",");
modelValues = arrayfun(@jsonString,document.modelFiles);
append("  " + jsonString("modelFiles") + ": [" + strjoin(modelValues,", ") + "],");
append("  " + jsonString("integration") + ": {");
integrationFields = string(fieldnames(document.integration));
for iField = 1:numel(integrationFields)
    field = integrationFields(iField);
    value = document.integration.(field);
    if isstring(value)
        valueText = jsonString(value);
    else
        valueText = numberText(value);
    end
    suffix = conditional(iField < numel(integrationFields),",","");
    append("    " + jsonString(field) + ": " + valueText + suffix);
end
append("  },");
append("  " + jsonString("output") + ": {");
append("    " + jsonString("policy") + ": " + jsonString(document.outputPolicy) + ",");
if isempty(destinationIdentifiers)
    append("    " + jsonString("destinations") + ": {}");
else
    append("    " + jsonString("destinations") + ": {");
    for iDestination = 1:numel(destinationIdentifiers)
        identifier = destinationIdentifiers(iDestination);
        suffix = conditional(iDestination < numel(destinationIdentifiers),",","");
        append("      " + jsonString(identifier) + ": " + jsonString(document.destinations(identifier)) + suffix);
    end
    append("    }");
end
append("  },");
append("  " + jsonString("execution") + ": {");
if strlength(document.fftProvider) > 0
    suffix = conditional(~isnan(document.threads),",","");
    append("    " + jsonString("fftProvider") + ": " + jsonString(document.fftProvider) + suffix);
end
if ~isnan(document.threads)
    append("    " + jsonString("threads") + ": " + integerText(document.threads));
end
append("  },");
append("  " + jsonString("report") + ": " + jsonString(document.report));
append("}");
encoded = strjoin(lines,newline) + newline;

    function append(line)
        lines(end+1,1) = line;
    end
end

function encoded = jsonString(value)
characters = char(string(value));
pieces = strings(1,numel(characters));
for iCharacter = 1:numel(characters)
    code = double(characters(iCharacter));
    switch code
        case 8
            pieces(iCharacter) = string(char([92 98]));
        case 9
            pieces(iCharacter) = string(char([92 116]));
        case 10
            pieces(iCharacter) = string(char([92 110]));
        case 12
            pieces(iCharacter) = string(char([92 102]));
        case 13
            pieces(iCharacter) = string(char([92 114]));
        case 34
            pieces(iCharacter) = string(char([92 34]));
        case 92
            pieces(iCharacter) = string(char([92 92]));
        otherwise
            if code < 32
                pieces(iCharacter) = string(char([92 117])) + lower(string(dec2hex(code,4)));
            else
                pieces(iCharacter) = string(characters(iCharacter));
            end
    end
end
quote = string(char(34));
encoded = quote + join(pieces,"") + quote;
end

function value = numberText(value)
value = double(value);
if ~isscalar(value) || ~isfinite(value)
    error("WaveVortexModel:PortableRunRequestEncoding","JSON numbers must be finite scalars.");
end
if value == 0
    value = 0;
end
value = string(sprintf("%.17g",value));
end

function value = integerText(value)
value = string(sprintf("%.0f",double(value)));
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function validateEncodedDocument(encoded,document)
try
    decoded = jsondecode(char(encoded));
catch exception
    throwAsCaller(addCause(MException("WaveVortexModel:PortableRunRequestValidation", ...
        "The encoded request is not valid JSON."),exception));
end
expectedFields = ["schemaIdentifier","schemaVersion","modelFiles","integration","output","execution","report"];
if ~isequal(string(fieldnames(decoded)),expectedFields(:)) || ...
        string(decoded.schemaIdentifier) ~= document.schemaIdentifier || ...
        double(decoded.schemaVersion) ~= document.schemaVersion || ...
        ~isequal(reshape(string(decoded.modelFiles),1,[]),document.modelFiles)
    error("WaveVortexModel:PortableRunRequestValidation","The encoded request did not pass schema-specific round-trip validation.");
end
end

function path = resolveDocumentPath(requestFolder,path)
path = string(path);
if ~isAbsolutePath(path)
    path = fullfile(requestFolder,path);
end
path = canonicalPath(path);
end

function tf = isAbsolutePath(path)
tf = javaObject("java.io.File",char(path)).isAbsolute();
end

function path = canonicalPath(path)
path = string(javaObject("java.io.File",char(path)).getCanonicalPath());
end

function bytes = readBytes(path)
fileIdentifier = fopen(path,"rb");
if fileIdentifier < 0
    error("WaveVortexModel:PortableRunRequestTemporaryFile","Unable to reopen the temporary request file.");
end
cleanup = onCleanup(@()fclose(fileIdentifier));
bytes = fread(fileIdentifier,Inf,"*uint8");
clear cleanup
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function atomicReplace(source,destination)
sourcePath = javaObject("java.io.File",char(source)).toPath();
destinationPath = javaObject("java.io.File",char(destination)).toPath();
atomicMove = javaMethod("valueOf","java.nio.file.StandardCopyOption","ATOMIC_MOVE");
replaceExisting = javaMethod("valueOf","java.nio.file.StandardCopyOption","REPLACE_EXISTING");
options = javaArray("java.nio.file.CopyOption",2);
options(1) = atomicMove;
options(2) = replaceExisting;
try
    javaMethod("move","java.nio.file.Files",sourcePath,destinationPath,options);
catch exception
    if ~contains(string(exception.message),"AtomicMoveNotSupportedException")
        rethrow(exception)
    end
    fallback = javaArray("java.nio.file.CopyOption",1);
    fallback(1) = replaceExisting;
    javaMethod("move","java.nio.file.Files",sourcePath,destinationPath,fallback);
end
end

function injectFailure(failureStage,stage)
if failureStage == stage
    error("WaveVortexModel:PortableRunRequestInjectedFailure","Injected %s failure.",stage);
end
end
