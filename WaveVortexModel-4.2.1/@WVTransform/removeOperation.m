function removeOperation(self,operation)
% Remove the exact registered operation and its cached outputs.
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    operation (1,1) WVOperation {mustBeNonempty}
end

if ~isKey(self.operationNameMap,operation.name) || self.operationNameMap{operation.name} ~= operation
    error("The requested operation is not registered with this transform.")
end
for iOutput = 1:operation.nVarOut
    outputName = operation.outputVariables(iOutput).name;
    if ~isKey(self.operationVariableNameMap,outputName) || self.operationVariableNameMap(outputName).modelOp ~= operation
        error("The registered output map for operation '%s' is inconsistent.",operation.name)
    end
end

self.removePropertyAnnotation(operation.outputVariables);
for iOutput = 1:operation.nVarOut
    outputName = operation.outputVariables(iOutput).name;
    self.operationVariableNameMap(outputName) = [];
    self.removeFromVariableCache(outputName);
end
self.operationNameMap(operation.name) = [];
end
