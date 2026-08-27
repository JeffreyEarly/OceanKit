function addFlowComponent(self,flowComponent)
% add a flow component and its standard variables
%
% - Topic: Flow components
% - Declaration: addFlowComponent(flowComponent)
% - Parameter flowComponent: one or more WVFlowComponent objects
%
% The standard variables supported by the concrete transform are
% registered for each component. These include the three-dimensional state
% variables and the sea-surface variables.
arguments
    self WVTransform {mustBeNonempty}
    flowComponent (1,:) WVFlowComponent {mustBeNonempty}
end

standardVariableNames = {'u','v','w','eta','p','ssh','ssu','ssv'};
supportedVariableNames = self.namesOfTransformVariables();
variableNames = standardVariableNames(ismember(standardVariableNames,supportedVariableNames));

for i=1:length(flowComponent)
    self.flowComponentNameMap{flowComponent(i).shortName} = flowComponent(i);
    self.addOperation(self.operationForKnownVariable(variableNames{:},flowComponent=flowComponent(i)));
end
end
