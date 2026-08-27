function writePortableRunRequestForTesting(path,modelFiles,configuration,failureStage)
% Exercise portable-request writing with an injected transactional failure.
%
% - Topic: Write model output
% - Declaration: WVModel.writePortableRunRequestForTesting(path,modelFiles,configuration,failureStage)
% - Developer: true

arguments
    path (1,1) string {mustBeNonzeroLengthText}
    modelFiles string {mustBeNonempty}
    configuration (1,1) struct
    failureStage (1,1) string {mustBeMember(failureStage,["","validation","encoding","temporary-file","replacement","cleanup"])}
end
WVInternal.writePortableRunRequest(path,modelFiles,configuration,failureStage);
end
