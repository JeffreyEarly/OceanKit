function [M_wwg, omegaAxis] = quadraticEnergyMirrorTriadFluxes1D_omega(self,options)
% Quadratic Energy Mirror Triad Fluxes1 D omega.
%
% quadraticEnergyMirrorTriadFluxes1D_omega is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
% - Declaration: [M_wwg, omegaAxis] = quadraticEnergyMirrorTriadFluxes1D_omega(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns M_wwg: output value `M_wwg`
% - Returns omegaAxis: output value `omegaAxis`
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end
% The mirror fluxes of the the [w{\nabla}w]_g triad
% components are computed on a custom sparse omega grid

val= self.diagfile.readVariables('pi_w_wwg_omega');
M_wwg= diff(cat(1,zeros(1,size(val,2)),val));

omegaAxis =  reshape(self.diagfile.readVariables('omegaAxis'),[],1);

if isinf(options.timeIndices)
    M_wwg = mean(M_wwg,2);
else
    M_wwg = mean(M_wwg(:,options.timeIndices),2);
end

end