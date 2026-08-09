function val = operationWithName(self,name)
% retrieve a WVOperation by name
%
% - Topic: Utility function — Metadata
arguments (Input)
    self WVTransform {mustBeNonempty}
    name char {mustBeNonempty}
end
arguments (Output)
    val WVOperation 
end
if ~isKey(self.operationNameMap,name)
    error("No operation named '%s' is registered with this transform.",name)
end
val = self.operationNameMap{name};
end
