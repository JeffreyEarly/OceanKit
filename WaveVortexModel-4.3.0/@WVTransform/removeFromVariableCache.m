function removeFromVariableCache(self,name)
% remove one variable from the internal cache
%
% - Topic: Internal
if isKey(self.variableCache,name)
    self.variableCache(name) = [];
end
end
