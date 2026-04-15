function [M_ggw, kePeAxis] = quadraticEnergyMirrorTriadFluxes1D_kepe(self,options)
% Quadratic Energy Mirror Triad Fluxes1 D kepe.
%
% quadraticEnergyMirrorTriadFluxes1D_kepe is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Diagnostics — Energy Fluxes — Temporal averages, 1D axes
% - Declaration: [M_ggw, kePeAxis] = quadraticEnergyMirrorTriadFluxes1D_kepe(self,options)
% - Parameter self: WVDiagnostics object
% - Parameter timeIndices: (optional) indices specifying which time indices to use (default: Inf)
% - Returns M_ggw: output value `M_ggw`
% - Returns kePeAxis: output value `kePeAxis`
arguments
    self WVDiagnostics
    options.timeIndices = Inf;
end

val= self.diagfile.readVariables('pi_g_ggw_kepe');
M_ggw= diff(cat(1,zeros(1,size(val,2)),val));

kePeAxis =  reshape(self.diagfile.readVariables('kePeAxis'),[],1);

if isinf(options.timeIndices)
    M_ggw = mean(M_ggw,2);
else
    M_ggw = mean(M_ggw(:,options.timeIndices),2);
end

end