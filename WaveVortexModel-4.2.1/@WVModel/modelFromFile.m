function model = modelFromFile(path)
% Initialize a model from an existing file
%
% A WVModel will be initialized from the specified path. The model will
% have this file designated as its only output file. All groups in that
% file are restored, but other files previously written by the same model
% are not reconstructed. The file must contain one complete coefficient
% stream. Linear versus nonlinear dynamics is restored when the model
% metadata is present; older files retain the nonlinear default.
%
% Runtime integrator objects are not persisted. Configure the desired
% fixed or adaptive integrator before continuing the restored model.
%
% - Topic: Initialization
% - Declaration: model = modelFromFile(path)
% - Parameter path: path to a NetCDF file
% - Returns model: restored `WVModel` instance
    arguments
        path char {mustBeFile}
    end

    [wvt,reader] = WVTransform.waveVortexTransformFromFile(path,iTime=Inf,shouldReadOnly=true);
    readerCleanup = onCleanup(@()closeNetCDFFileIfOpen(reader));
    isDynamicsLinear = false;
    if isKey(reader.attributes,'WVModelIsDynamicsLinear')
        isDynamicsLinear = logical(reader.attributes('WVModelIsDynamicsLinear'));
    end
    model = WVModel(wvt,shouldUseLinearDynamics=isDynamicsLinear);
    reader.close();
    clear readerCleanup
    ncfile = NetCDFFile(path,shouldReadOnly=false);
    try
        outputFile = WVModelOutputFile.modelOutputFileFromFile(ncfile,model);
        model.addOutputFile(outputFile);
    catch exception
        if ~isempty(ncfile.id)
            ncfile.close();
        end
        rethrow(exception)
    end
end


function closeNetCDFFileIfOpen(ncfile)
    if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
        ncfile.close();
    end
end
