function varargout = variableWithName(self,variableNames)
% Compute or retrieve one or more registered transform variables.
%
% - Topic: State variables
% - Declaration: varargout = variableWithName(variableNames)
% - Parameter variableNames: names of registered state variables
% - Returns varargout: state-variable arrays in the requested order
arguments
    self WVTransform {mustBeNonempty}
end
arguments (Repeating)
    variableNames char
end

while ~all(isKey(self.variableCache,variableNames))
    missingIndex = find(~isKey(self.variableCache,variableNames),1);
    missingName = variableNames{missingIndex};
    if ~isKey(self.operationVariableNameMap,missingName)
        error("No variable named '%s' is registered with this transform.",missingName)
    end
    annotation = self.operationVariableNameMap(missingName);
    self.performOperation(annotation.modelOp);
end

varargout = cell(size(variableNames));
[varargout{:}] = self.fetchFromVariableCache(variableNames{:});
end
