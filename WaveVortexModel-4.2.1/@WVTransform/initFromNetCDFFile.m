function initFromNetCDFFile(wvt,ncfile,options)
% initialize the flow from a NetCDF file
%
% Clears variables Ap,Am,A0 and then sets them the values found in the file
% at the requested time.
% 
% This is intended to be used in conjunction with
% [`waveVortexTransformFromFile`](/classes/transforms/wvtransform/wavevortextransformfromfile.html)
% e.g.,
%
% ```matlab
% [wvt,ncfile] = WVTransform.waveVortexTransformFromFile('cyprus-eddy.nc');
% t = ncfile.readVariables('t');
% for iTime=1:length(t)
%     wvt.initFromNetCDFFile(ncfile,iTime=iTime)
%     // some analysis
% end
% ```
%
% Note that this method only lightly checks that you are reading from a
% file that is compatible with this transform! So be careful.
%
% See also the users guide for [reading and writing to
% file](/users-guide/reading-and-writing-to-file.html).
% 
% - Topic: Initial conditions
% - Declaration: initWithFile(self,path,options)
% - Parameter ncfile: a NetCDF file object
% - Parameter iTime: (optional) time index to initialize from (default 1)
arguments
    wvt WVTransform {mustBeNonempty}
    ncfile NetCDFFile {mustBeNonempty}
    options.iTime (1,1) double {mustBePositive} = 1
    options.shouldDisplayInit (1,1) = 0
end

wvt.t0 = ncfile.readVariables('t0');

coefficientNames = {};
if wvt.hasWaveComponent
    coefficientNames = {'Ap','Am'};
end
if wvt.hasPVComponent
    coefficientNames{end+1} = 'A0';
end
stateVariableNames = {};
if ~isempty(coefficientNames) && all(ncfile.hasVariableWithName(coefficientNames{:}))
    stateVariableNames = coefficientNames;
elseif all(ncfile.hasVariableWithName('u','v','eta'))
    stateVariableNames = {'u','v','eta'};
elseif ncfile.hasVariableWithName('A0')
    stateVariableNames = {'A0'};
end
stateGroup = ncfile;
if ~isempty(stateVariableNames)
    stateGroup = commonVariableGroup(ncfile,stateVariableNames);
end

hasTimeDimension = 0;
if stateGroup.hasDimensionWithName('t') && stateGroup.hasVariableWithName('t')
    hasTimeDimension = 1;
    tDim = stateGroup.readVariables('t');
    if isinf(options.iTime)
        iTime = length(tDim);
    elseif options.iTime > length(tDim)
        error('Index out of bounds! There are %d time points in this file, you requested %d.',length(tDim),options.iTime);
    else
        iTime = options.iTime;
    end
    wvt.t = tDim(iTime);
else
    wvt.t = stateGroup.readVariables('t');
end

if all(ismember({'Ap','Am','A0'},stateVariableNames))
    if hasTimeDimension == 1
        [wvt.A0,wvt.Ap,wvt.Am] = stateGroup.readVariablesAtIndexAlongDimension('t',iTime,'A0','Ap','Am');
    else
        [wvt.A0,wvt.Ap,wvt.Am] = stateGroup.readVariables('A0','Ap','Am');
    end
    if options.shouldDisplayInit == 1
        fprintf('%s initialized from Ap, Am, A0.\n',ncfile.attributes('WVTransform'));
    end
elseif all(ismember({'u','v','eta'},stateVariableNames))
    if hasTimeDimension == 1
        [u,v,eta] = stateGroup.readVariablesAtIndexAlongDimension('t',iTime,'u','v','eta');
    else
        [u,v,eta] = stateGroup.readVariables('u','v','eta');
    end
    [wvt.Ap,wvt.Am,wvt.A0] = wvt.transformUVEtaToWaveVortex(u,v,eta);
    if options.shouldDisplayInit == 1
        fprintf('%s initialized from u, u, eta.\n',ncfile.attributes('WVTransform'));
    end
elseif isequal(stateVariableNames,{'A0'})
    if hasTimeDimension == 1
        wvt.A0 = reshape(stateGroup.readVariablesAtIndexAlongDimension('t',iTime,'A0'),wvt.spectralMatrixSize);
    else
        wvt.A0 = stateGroup.readVariables('A0');
    end
    if options.shouldDisplayInit == 1
        fprintf('%s initialized from A0.\n',ncfile.attributes('AnnotatedClass'));
    end
else
    warning('%s initialized without data.\n',ncfile.attributes('AnnotatedClass'));
end

end

function group = commonVariableGroup(ncfile,variableNames)
variables = cell(1,length(variableNames));
for iVariable = 1:length(variableNames)
    paths = ncfile.variablePathsWithName(variableNames{iVariable});
    if length(paths) ~= 1
        error('A restart file must contain exactly one %s state variable, but %d were found.',variableNames{iVariable},length(paths));
    end
    variables{iVariable} = ncfile.variableWithName(char(paths));
end
group = variables{1}.group;
for iVariable = 2:length(variables)
    if variables{iVariable}.group ~= group
        error('The restart state variables must be stored together in one NetCDF group.');
    end
end
if ~group.hasVariableWithName('t')
    error('The NetCDF group containing the restart state does not contain time.');
end
end
