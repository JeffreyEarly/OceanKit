function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
% Initialize a WVTransform instance from an existing file
%
% A WVTransform instance can be recreated from a NetCDF file and .mat
% sidecar file if the default variables were save to file. For example,
%
% ```matlab
% wvt = WVTransform.waveVortexTransformFromFile('cyprus-eddy.nc',iTime=Inf);
% ```
%
% will create a WVTransform instance populated with values from the last
% time point of the file. Note that this is a static function---it is a
% function defined on the class, not an instance variable---so requires we
% prepend `WVTransform.` The result of this function call is an instance
% variable.
%
% Restoring only the transform closes the NetCDF file before returning. If
% you intend to read more than one time point from the save file, request
% the second output and hold onto that read-only NetCDFFile instance, then call
% [`initFromNetCDFFile`](/classes/transforms/wvtransform/initfromnetcdffile.html). This
% avoids the relatively expensive operation recreating the WVTransform, and
% simply reads the appropriate data from file. The caller owns the returned
% NetCDFFile and must close it. Set `shouldReadOnly=false` only when writable
% access is required.
%
% See also the users guide for [reading and writing to
% file](/users-guide/reading-and-writing-to-file.html).
%
% - Topic: Initialization
% - Declaration: [wvt,ncfile] = WVTransform.waveVortexTransformFromFile(path,options)
% - Parameter path: path to a NetCDF file
% - Parameter iTime: (optional) time index to initialize from (default 1).
% - Parameter shouldReadOnly: (optional) open the returned NetCDFFile read-only (default true).
% - Returns wvt: an instance of a WVTransform subclass
% - Returns ncfile: a caller-owned NetCDFFile instance pointing to the file
arguments (Input)
    path char {mustBeFile}
    options.iTime (1,1) double {mustBePositive} = 1
    options.shouldReadOnly logical = true
end
arguments (Output)
    wvt WVTransform
    ncfile NetCDFFile
end
wvtClassName = transformClassNameFromFile(path);
[wvt,ncfile] = feval(strcat(wvtClassName,'.waveVortexTransformFromFile'),path,'iTime',options.iTime,'shouldReadOnly',options.shouldReadOnly);

if nargout < 2
    ncfile.close();
end

% totalForcingGroups = ncfile.attributes('TotalForcingGroups');
% for iForce=1:totalForcingGroups
%     forceGroup = ncfile.groupWithName("forcing-"+iForce);
%     forcingClassName = forceGroup.attributes('WVForcing');
%     force = feval(strcat(forcingClassName,'.forcingFromFile'),forceGroup,wvt);
%     self.addForcing(force);
% end

end

function wvtClassName = transformClassNameFromFile(path)
ncfile = NetCDFFile(path,shouldReadOnly=true);
cleanup = onCleanup(@()closeNetCDFFileIfOpen(ncfile));

if ~isKey(ncfile.attributes,'WVTransform')
    error("WVTransform:MissingTransformClass","Unable to find the WVTransform attribute in %s. This file may not have been created by WaveVortexModel.",path);
end

wvtClassName = ncfile.attributes('WVTransform');
if ~(ischar(wvtClassName) || (isstring(wvtClassName) && isscalar(wvtClassName)))
    error("WVTransform:InvalidTransformClass","The WVTransform attribute in %s must contain a class name.",path);
end
wvtClassName = char(wvtClassName);
if exist(wvtClassName,'class') ~= 8
    error("WVTransform:InvalidTransformClass","The WVTransform attribute in %s names unavailable class %s.",path,wvtClassName);
end
end

function closeNetCDFFileIfOpen(ncfile)
if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
    ncfile.close();
end
end
