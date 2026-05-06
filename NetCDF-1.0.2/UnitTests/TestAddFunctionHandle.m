testPath = fullfile(fileparts(mfilename('fullpath')), 'test_function_handle.nc');
cleanup = onCleanup(@() deleteIfPresent(testPath));
deleteIfPresent(testPath);

%%
m = 1; b = 2;
f = @(x) m*x + b;

ncfile = NetCDFFile(testPath,shouldOverwriteExisting=1);
ncfile.addFunctionHandle('f',f);

%%
ncfile.close();
clear ncfile f m b
ncfile = NetCDFFile(testPath);
f = ncfile.readVariables('f');
ncfile.close();

%%
cleanupToken = onCleanup(@() []);
f = @(x) cleanupToken;
ncfile = NetCDFFile(testPath,shouldOverwriteExisting=1);
ncfile.addFunctionHandle('f',f);
ncfile.close();
clear ncfile f
ncfile = NetCDFFile(testPath);
f = ncfile.readVariables('f');
assert(isa(f, 'function_handle'), 'Warning-only function handles should still round-trip through NetCDF persistence.');
ncfile.close();

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
