function addOperation(self,operation,options)
% Register one or more operations and their output variables.
%
% The complete request is validated before annotations, lookup maps, or
% cached values are changed. Existing operations are replaced only when
% `shouldOverwriteExisting` is true. For an operation array, replacements
% are evaluated in caller order and the last conflicting operation wins.
%
% - Topic: Utility function — Metadata
arguments
    self WVTransform {mustBeNonempty}
    operation (1,:) WVOperation {mustBeNonempty}
    options.shouldOverwriteExisting logical = false
    options.shouldSuppressWarning logical = false
end

registeredOperationNames = self.operationNameMap.keys;
registeredOperations = cell(1,length(registeredOperationNames));
for iOperation = 1:length(registeredOperationNames)
    registeredOperations{iOperation} = self.operationNameMap{registeredOperationNames(iOperation)};
end
stagedOperations = [registeredOperations cell(1,length(operation))];
stagedOperationCount = length(registeredOperations);

for iOperation = 1:length(operation)
    candidate = operation(iOperation);
    outputNames = string({candidate.outputVariables.name});
    if length(unique(outputNames)) ~= length(outputNames)
        error("Operation '%s' declares duplicate output-variable names.",candidate.name)
    end

    for iOutput = 1:candidate.nVarOut
        annotation = candidate.outputVariables(iOutput);
        unknownDimensions = setdiff(string(annotation.dimensions),string(self.annotatedDimensionNames));
        if ~isempty(unknownDimensions)
            error("Operation '%s' output '%s' refers to unknown dimensions: %s.",candidate.name,annotation.name,join(unknownDimensions,", "))
        end
        if ismember(annotation.name,self.annotatedPropertyNames) && ~isKey(self.operationVariableNameMap,annotation.name)
            existingAnnotation = self.propertyAnnotationWithName(annotation.name);
            if ~isa(existingAnnotation,"WVVariableAnnotation")
                error("Operation '%s' cannot replace the non-operation annotation '%s'.",candidate.name,annotation.name)
            end
        end
    end

    conflictingIndices = false(1,stagedOperationCount);
    for iExisting = 1:stagedOperationCount
        existing = stagedOperations{iExisting};
        existingOutputNames = string({existing.outputVariables.name});
        conflictingIndices(iExisting) = strcmp(existing.name,candidate.name) || any(ismember(outputNames,existingOutputNames));
    end
    if any(conflictingIndices) && ~options.shouldOverwriteExisting
        activeOperations = stagedOperations(1:stagedOperationCount);
        conflictingNames = string(cellfun(@(existing) existing.name,activeOperations(conflictingIndices),UniformOutput=false));
        error("Operation '%s' conflicts with registered operation(s): %s. Set shouldOverwriteExisting=true to replace them.",candidate.name,join(conflictingNames,", "))
    end

    remainingOperations = stagedOperations(1:stagedOperationCount);
    remainingOperations(conflictingIndices) = [];
    stagedOperationCount = length(remainingOperations)+1;
    stagedOperations(1:stagedOperationCount-1) = remainingOperations;
    stagedOperations{stagedOperationCount} = candidate;
end
stagedOperations = stagedOperations(1:stagedOperationCount);

removedOperations = operationsAbsentFrom(registeredOperations,stagedOperations);
addedOperations = operationsAbsentFrom(stagedOperations,registeredOperations);
if isempty(removedOperations) && isempty(addedOperations)
    return
end

removedAnnotations = outputAnnotationsForOperations(removedOperations);
if ~isempty(removedAnnotations)
    self.removePropertyAnnotation(removedAnnotations);
end
for iOperation = 1:length(removedOperations)
    removed = removedOperations{iOperation};
    for iOutput = 1:removed.nVarOut
        outputName = removed.outputVariables(iOutput).name;
        if isKey(self.operationVariableNameMap,outputName)
            self.operationVariableNameMap(outputName) = [];
        end
        self.removeFromVariableCache(outputName);
    end
    if isKey(self.operationNameMap,removed.name)
        self.operationNameMap(removed.name) = [];
    end
end

addedAnnotations = outputAnnotationsForOperations(addedOperations);
if ~isempty(addedAnnotations)
    self.addPropertyAnnotation(addedAnnotations);
end
for iOperation = 1:length(addedOperations)
    added = addedOperations{iOperation};
    for iOutput = 1:added.nVarOut
        outputVariable = added.outputVariables(iOutput);
        self.removeFromVariableCache(outputVariable.name);
        self.operationVariableNameMap(outputVariable.name) = outputVariable;
    end
    self.operationNameMap{added.name} = added;
end

if ~options.shouldSuppressWarning && ~isempty(removedOperations)
    removedNames = string(cellfun(@(removed) removed.name,removedOperations,UniformOutput=false));
    addedNames = string(cellfun(@(added) added.name,addedOperations,UniformOutput=false));
    fprintf("Replaced operation(s) {%s} with {%s}.\n",join(removedNames,", "),join(addedNames,", "))
end
end

function absentOperations = operationsAbsentFrom(candidates,reference)
absentOperations = cell(1,length(candidates));
absentOperationCount = 0;
for iCandidate = 1:length(candidates)
    candidate = candidates{iCandidate};
    isPresent = any(cellfun(@(other) other == candidate,reference));
    if ~isPresent
        absentOperationCount = absentOperationCount+1;
        absentOperations{absentOperationCount} = candidate;
    end
end
absentOperations = absentOperations(1:absentOperationCount);
end

function annotations = outputAnnotationsForOperations(operations)
if isempty(operations)
    annotations = WVVariableAnnotation.empty(1,0);
    return
end
annotationCount = sum(cellfun(@(operation) operation.nVarOut,operations));
annotations = repmat(operations{1}.outputVariables(1),1,annotationCount);
nextAnnotation = 1;
for iOperation = 1:length(operations)
    operationAnnotations = operations{iOperation}.outputVariables;
    annotationIndices = nextAnnotation:nextAnnotation+length(operationAnnotations)-1;
    annotations(annotationIndices) = operationAnnotations;
    nextAnnotation = annotationIndices(end)+1;
end
end
