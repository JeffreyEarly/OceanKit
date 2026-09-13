function varargout = performOperation(self,modelOp)
% computes (runs) the operation
%
% - Topic: Internal
arguments (Input)
    self WVTransform
    modelOp WVOperation
end
arguments (Output,Repeating)
    varargout
end

varNames = cell(1,modelOp.nVarOut);
varargout = cell(1,modelOp.nVarOut);
for iVar=1:length(modelOp.outputVariables)
    varNames{iVar} = modelOp.outputVariables(iVar).name;
end

if all(isKey(self.variableCache,varNames))
    [varargout{:}] = self.fetchFromVariableCache(varNames{:});
else
    scope = self.scopedEvaluation(); %#ok<NASGU>
    if self.isCompiledBuiltinOperation(modelOp) && ~isa(modelOp,'WVNoMotionProfileOperation')
        values = self.compiledVariables(varNames);
        varargout = reshape(values,size(varargout));
    else
        [varargout{:}] = modelOp.compute(self);
    end
    for iOpOut=1:length(varargout)
        self.addToVariableCache(varNames{iOpOut},varargout{iOpOut})
    end
end
end
