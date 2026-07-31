function model = modelFromFile(path)
% Initialize a model from an existing file
%
% A WVModel will be initialized from the specified path. The model will
% have this file designated as its outputFile. Integrating the model will
% thus write to this file.
%
% - Topic: Initialization
% - Declaration: model = modelFromFile(path)
% - Parameter path: path to a NetCDF file
    arguments
        path char {mustBeFile}
    end

    wvt = WVTransform.waveVortexTransformFromFile(path,iTime=Inf,shouldReadOnly=true);
    model = WVModel(wvt);
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
