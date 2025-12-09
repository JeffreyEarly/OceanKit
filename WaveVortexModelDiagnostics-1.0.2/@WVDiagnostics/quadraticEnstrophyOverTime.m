function [enstrophy, t] = quadraticEnstrophyOverTime(self, options)
% Quadratic Enstrophy Over Time.
%
% quadraticEnstrophyOverTime is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Potential Enstrophy — Time series, [t 1]
% - Declaration: [enstrophy, t] = quadraticEnstrophyOverTime(self, options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns enstrophy: output value `enstrophy`
% - Returns t: Summary table of enstrophy flux diagnostics
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
if isinf(options.timeIndices)
    filter = @(v) v;
else
    filter = @(v) v(options.timeIndices);
end
enstrophy = filter(self.diagfile.readVariables('enstrophy_quadratic'));
t = filter(self.diagfile.readVariables('t'));
end