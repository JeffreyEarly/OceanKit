function ncfile = writeToFile(self, path, properties, options)
% Write this bootstrap ensemble to a NetCDF restart file.
%
% `writeToFile` persists the canonical bootstrap restart state. If the
% lazy `bootstrapKappa` or `bootstrapCoherence` vectors have already been
% materialized, they are written as optional cached diagnostics as well.
%
% - Topic: Write to file
% - Declaration: ncfile = writeToFile(self,path,properties,options)
% - Parameter path: path to write file
% - Parameter properties: optional extra property names to persist
% - Parameter shouldOverwriteExisting: (optional) boolean indicating whether or not to overwrite an existing file at the path. Default `false`.
% - Parameter shouldAddRequiredProperties: (optional) boolean indicating whether or not to include the canonical required restart properties. Default `true`.
% - Parameter attributes: (optional) dictionary containing additional NetCDF attributes to write
% - Returns ncfile: open `NetCDFFile` handle for the written restart file
arguments (Input)
    self (1,1) GriddedStreamfunctionBootstrap
    path {mustBeTextScalar}
end
arguments (Input, Repeating)
    properties char
end
arguments (Input)
    options.shouldOverwriteExisting logical = false
    options.shouldAddRequiredProperties logical = true
    options.attributes = configureDictionary("string", "string")
end
arguments (Output)
    ncfile NetCDFFile
end

path = char(path);
propertyNames = properties;
if self.nBootstraps > 0 && numel(self.bootstrapKappaCache) == self.nBootstraps
    propertyNames{end + 1} = 'bootstrapKappaCache'; %#ok<AGROW>
end
if self.nBootstraps > 0 && numel(self.bootstrapCoherenceCache) == self.nBootstraps
    propertyNames{end + 1} = 'bootstrapCoherenceCache'; %#ok<AGROW>
end

ncfile = writeToFile@CAAnnotatedClass(self, path, propertyNames{:}, ...
    shouldOverwriteExisting=options.shouldOverwriteExisting, ...
    shouldAddRequiredProperties=options.shouldAddRequiredProperties, ...
    attributes=options.attributes);
end
