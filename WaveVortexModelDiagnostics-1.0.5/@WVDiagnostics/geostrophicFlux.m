function [spatialFlux, spectralFlux] = geostrophicFlux(self)
% Compute the spatial and spectral fluxes for constant stratification
%
% This is from the source-vector form of the fluxes, derived by Jonathan
% Lilly and Jeffrey Early, and is only valid for constant stratification at
% the moment.
%
%
% - Topic: Diagnostics — General — Misc — Fluxes
% - Declaration: [spatialFlux, spectralFlux] = geostrophicFlux(self)
% - Parameter self: WVDiagnostics object
% - Returns spatialFlux: struct containing `ggg`, `ggw`, `ggw_tx`, `wwg_tx`
% - Returns spectralFlux: struct containing `ggg`, `ggw`, `ggw_tx`, `wwg_tx`
arguments
    self WVDiagnostics
end
if self.diagnosticsHasExplicitAntialiasing
    wvt = self.wvt_aa;
else
    wvt = self.wvt;
end

%ggg, ggw, ggw_tx, gww_tx

%% Source-vector g
chi_1 = wvt.v_g;
chi_2 = -wvt.u_g;
chi_3 = -wvt.f*(wvt.eta_g + wvt.eta_mda);

% ggg cascade
u = wvt.u_g;
v = wvt.v_g;
w = wvt.w_g;

t1 = wvt.diffX( ( u .* wvt.diffX(chi_1) + v .* wvt.diffY(chi_1) + w .* wvt.diffZF(chi_1) ) );
t2 = wvt.diffY( ( u .* wvt.diffX(chi_2) + v .* wvt.diffY(chi_2) + w .* wvt.diffZF(chi_2) ) );
t3 = wvt.diffZG(( u .* wvt.diffX(chi_3) + v .* wvt.diffY(chi_3) + w .* wvt.diffZG(chi_3) ) );

spatialFlux.ggg = wvt.psi .* (t1 + t2 + t3);
spectralFlux.ggg = wvt.crossSpectrumWithFgTransform(wvt.psi, t1 + t2 + t3);

% ggw cascade
u = wvt.u_w + wvt.u_io;
v = wvt.v_w + wvt.v_io;
w = wvt.w_w;

t1 = wvt.diffX( ( u .* wvt.diffX(chi_1) + v .* wvt.diffY(chi_1) + w .* wvt.diffZF(chi_1) ) );
t2 = wvt.diffY( ( u .* wvt.diffX(chi_2) + v .* wvt.diffY(chi_2) + w .* wvt.diffZF(chi_2) ) );
t3 = wvt.diffZG(( u .* wvt.diffX(chi_3) + v .* wvt.diffY(chi_3) + w .* wvt.diffZG(chi_3) ) );

spatialFlux.ggw = wvt.psi .* (t1 + t2 + t3);
spectralFlux.ggw = wvt.crossSpectrumWithFgTransform(wvt.psi, t1 + t2 + t3 );

%% Source-vector w
chi_1 = wvt.v_w + wvt.v_io;
chi_2 = -(wvt.u_w + wvt.u_io);
chi_3 = -wvt.f*(wvt.eta_w + wvt.eta_io);

% ggw transfer
u = wvt.u_g;
v = wvt.v_g;
w = wvt.w_g;

t1 = wvt.diffX( ( u .* wvt.diffX(chi_1) + v .* wvt.diffY(chi_1) + w .* wvt.diffZF(chi_1) ) );
t2 = wvt.diffY( ( u .* wvt.diffX(chi_2) + v .* wvt.diffY(chi_2) + w .* wvt.diffZF(chi_2) ) );
t3 = wvt.diffZG(( u .* wvt.diffX(chi_3) + v .* wvt.diffY(chi_3) + w .* wvt.diffZG(chi_3) ) );

spatialFlux.ggw_tx = wvt.psi .* (t1 + t2 + t3);
spectralFlux.ggw_tx = wvt.crossSpectrumWithFgTransform(wvt.psi, t1 + t2 + t3);

% wwg transfer
u = wvt.u_w + wvt.u_io;
v = wvt.v_w + wvt.v_io;
w = wvt.w_w;

t1 = wvt.diffX( ( u .* wvt.diffX(chi_1) + v .* wvt.diffY(chi_1) + w .* wvt.diffZF(chi_1) ) );
t2 = wvt.diffY( ( u .* wvt.diffX(chi_2) + v .* wvt.diffY(chi_2) + w .* wvt.diffZF(chi_2) ) );
t3 = wvt.diffZG(( u .* wvt.diffX(chi_3) + v .* wvt.diffY(chi_3) + w .* wvt.diffZG(chi_3) ) );

spatialFlux.wwg_tx = wvt.psi .* (t1 + t2 + t3);
spectralFlux.wwg_tx = wvt.crossSpectrumWithFgTransform(wvt.psi, t1 + t2 + t3);

end