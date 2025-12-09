function [energy, t] = exactEnergyOverTime(self, options)
% Exact Energy Over Time.
%
% exactEnergyOverTime is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy — Time series, [t 1]
% - Declaration: [energy, t] = exactEnergyOverTime(self, options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns energy: diagnosed energy as a function of time and/or scale
% - Returns t: Summary table of energy flux diagnostics
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter = @(v) v;
else
    filter = @(v) v(options.timeIndices);
end
[ke,ape] =self.diagfile.readVariables('ke','ape');
energy = filter(ke+ape);
t = filter(self.diagfile.readVariables('t'));
end