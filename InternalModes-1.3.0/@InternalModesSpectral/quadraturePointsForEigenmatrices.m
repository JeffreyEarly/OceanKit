function z_g = quadraturePointsForEigenmatrices(self,nPoints,A,B)
% Return quadrature points inferred from a generalized EVP matrix pair.
%
% - Topic: Developer topics
% - Developer: true
% - Declaration: z_g = quadraturePointsForEigenmatrices(self,nPoints,A,B)
% - Parameter self: InternalModesSpectral instance
% - Parameter nPoints: number of quadrature points requested
% - Parameter A: left generalized-eigenproblem matrix
% - Parameter B: right generalized-eigenproblem matrix
% - Returns z_g: depth locations of the quadrature points

if 2*nPoints >= self.nEVP
    error('InternalModesSpectral:NeedMorePoints', ...
        'You need at least twice as many nEVP points as quadrature points requested.');
end

[G_cheb,h] = self.solveGEP(A, B);
z_g = self.quadraturePointsForModes(nPoints, G_cheb, h);
end
