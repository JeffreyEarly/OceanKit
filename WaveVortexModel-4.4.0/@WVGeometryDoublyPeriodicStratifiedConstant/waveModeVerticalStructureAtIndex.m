function [F,dGdz] = waveModeVerticalStructureAtIndex(self,iZ)
% Return wave vertical-structure factors at one vertical grid index.
%
% For a wave coefficient matrix `Apm`, the horizontal Fourier values of
% the corresponding $$F$$ field and the vertical derivative of the
% corresponding $$G$$ field at `z(iZ)` are
%
% $$
% \widehat F(z_{iZ}) = \sum_j F_{j\boldsymbol{k}}(z_{iZ})A_{j\boldsymbol{k}},
% \qquad
% \partial_z\widehat G(z_{iZ}) =
% \sum_j \partial_zG_{j\boldsymbol{k}}(z_{iZ})A_{j\boldsymbol{k}}.
% $$
%
% Both returned arrays have dimensions `[Nj Nkl]`. Request only `F` when
% the derivative factors are not needed.
%
% - Topic: Operations — Transformations
% - Declaration: [F,dGdz] = waveModeVerticalStructureAtIndex(iZ)
% - Parameter iZ: integer vertical grid index between 1 and `Nz`
% - Returns F: wave $$F$$ factors at `z(iZ)`
% - Returns dGdz: wave $$\partial_zG$$ factors at `z(iZ)`
arguments (Input)
    self WVGeometryDoublyPeriodicStratifiedConstant {mustBeNonempty}
    iZ (1,1) double
end

if ~isfinite(iZ) || iZ ~= round(iZ) || iZ < 1 || iZ > self.Nz
    error("WVStratification:InvalidVerticalIndex", "iZ must be an integer between 1 and Nz (%d).", self.Nz)
end

F = self.iDCT(iZ,:).'.*(self.F_g./self.F_wg);
if nargout < 2
    return
end

verticalWavenumber = pi*self.j/self.Lz;
dGdz = self.iDCT(iZ,:).'.*verticalWavenumber.*(self.G_g./self.G_wg);
end
