classdef IMGeostrophicZeroAPVModesBasis
    % Store canonical or rotated geostrophic zero-APV modes.
    %
    % The canonical columns have unit endpoint responses. For every
    % wavenumber page, the selected response matrix is the identity. The
    % diagnostic structure is
    %
    % $$
    % G(z)=-\frac{g}{N^2(z)}\frac{\partial F}{\partial z}(z).
    % $$
    %
    % `F(z)` and `G(z)` have dimensions `nZ x nEndpoints x nK`. The four
    % public quadratic-form properties have dimensions
    % `nEndpoints x nEndpoints x nK`. Rotation methods apply the same
    % boundary-coordinate matrix to both physical variables.
    %
    % If $$\mathbf b_{\mathrm s}$$ and $$\mathbf b_{\mathrm d}$$ are the
    % surface and bottom response rows and $$\mathbf F_{\mathrm s}$$ and
    % $$\mathbf F_{\mathrm d}$$ are the corresponding value rows, then
    % each canonical wavenumber page stores
    %
    % $$
    % \mathsf R_B=\mathbf B^T\mathbf B,
    % \qquad
    % \mathsf H=\frac{f_0^2}{2gk^4}
    % \operatorname{sym}\left(-\mathbf F_{\mathrm s}^T\mathbf b_{\mathrm s}+\mathbf F_{\mathrm d}^T\mathbf b_{\mathrm d}\right),
    % $$
    %
    % $$
    % \mathsf B_0=\frac{f_0^2}{2g^2k^4}\mathbf b_{\mathrm s}^T\mathbf b_{\mathrm s},
    % \qquad
    % \mathsf B_d=\frac{f_0^2}{2g^2k^4}\mathbf b_{\mathrm d}^T\mathbf b_{\mathrm d}.
    % $$
    %
    % An absent endpoint has a zero form matrix. Generalized energy is
    % $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.
    %
    % ```matlab
    % boundaryModes = solver.solveGeostrophicZeroAPVModes(problem);
    % depthModes = boundaryModes.rotateBoundaryDepth(g0=-0.035,gd=0);
    % surfaceModes = boundaryModes.rotateSurfaceBuoyancy(g0=-0.035,gd=0);
    % Hg = boundaryModes.generalizedEnergyMatrix(g0=-0.035,gd=0);
    % customModes = boundaryModes.rotateWithPencil(name="custom",leftMatrix=Hg,rightMatrix=boundaryModes.endpointResponseMetric);
    % ```
    %
    % - Topic: Evaluate geostrophic zero-APV modes
    % - Topic: Inspect geostrophic zero-APV modes
    % - Topic: Form geostrophic zero-APV quadratic forms
    % - Topic: Rotate geostrophic zero-APV modes
    % - Declaration: classdef IMGeostrophicZeroAPVModesBasis

    properties (SetAccess = private)
        % Geostrophic zero-APV problem descriptor.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        problem

        % Numerical solver used to compute the canonical structures.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        solver

        % Horizontal-wavenumber pages.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        k

        % Canonically ordered endpoint-coordinate labels.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        endpoints

        % Endpoint-response Gram matrix $$\mathsf R_B$$.
        %
        % Canonical unit boundary coordinates give $$\mathsf R_B=\mathsf I$$.
        % Rotated bases expose the congruence-transformed matrix.
        %
        % - Topic: Form geostrophic zero-APV quadratic forms
        endpointResponseMetric

        % Physical-energy matrix $$\mathsf H$$.
        %
        % `energyMatrix(:,:,iK)` is the physical-energy bilinear form in
        % the current boundary coordinates for `k(iK)`.
        %
        % - Topic: Form geostrophic zero-APV quadratic forms
        energyMatrix

        % Surface-buoyancy matrix $$\mathsf B_0$$.
        %
        % This matrix is exactly zero when the surface endpoint is absent.
        %
        % - Topic: Form geostrophic zero-APV quadratic forms
        surfaceBuoyancyMatrix

        % Bottom-buoyancy matrix $$\mathsf B_d$$.
        %
        % This matrix is exactly zero when the bottom endpoint is absent.
        %
        % - Topic: Form geostrophic zero-APV quadratic forms
        bottomBuoyancyMatrix

        % Name of the current boundary-coordinate rotation.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        rotationName

        % Canonical-to-current column transformation.
        %
        % `rotationMatrix` has dimensions
        % `nEndpoints x nEndpoints x nK` and acts identically on `F` and
        % `G`.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        rotationMatrix

        % Pagewise generalized-pencil eigenvalues.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        rotationEigenvalues

        % Pagewise quadratic-form signatures.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        signatures

        % Signed boundary depths.
        %
        % `h0` is populated only when `rotationName` is
        % `"boundaryDepth"`, where $$h_0^a=2k^2\gamma_a$$.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        h0

        % Rotation normalization convention.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        normalizationConvention

        % Pagewise pencil and normalization residuals.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        rotationResiduals

        % Coriolis parameter.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        f0

        % Gravitational acceleration.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        g

        % Surface endpoint convention.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        surfaceBoundary

        % Physical vertical domain.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        N2

        % Additional creation and rotation metadata.
        %
        % - Topic: Inspect geostrophic zero-APV modes
        metadata
    end

    properties (Access = private)
        canonicalNativeModes
        currentNativeModes
        canonicalEndpointResponseMetric
        canonicalEnergyMatrix
        canonicalSurfaceBuoyancyMatrix
        canonicalBottomBuoyancyMatrix
    end

    methods
        function self = IMGeostrophicZeroAPVModesBasis(options)
            % Create a canonical boundary-normalized zero-APV basis.
            %
            % - Topic: Evaluate geostrophic zero-APV modes
            % - Declaration: basisSet = IMGeostrophicZeroAPVModesBasis(options)
            % - Parameter options.problem: geostrophic zero-APV problem
            % - Parameter options.solver: configured numerical solver
            % - Parameter options.nativeModes: canonical native structures
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: canonical boundary-normalized basis
            arguments
                options.problem IMGeostrophicZeroAPVModes
                options.solver IMSolver
                options.nativeModes double {mustBeReal, mustBeFinite}
                options.metadata struct = struct()
            end

            nEndpoints = numel(options.problem.endpoints);
            nK = numel(options.problem.k);
            if size(options.nativeModes,2) ~= nEndpoints || size(options.nativeModes,3) ~= nK
                error("IMGeostrophicZeroAPVModesBasis:InvalidModeShape", "nativeModes must have dimensions nNative x nEndpoints x nK.");
            end

            self.problem = options.problem;
            self.solver = options.solver;
            self.k = options.problem.k;
            self.endpoints = options.problem.endpoints;
            self.f0 = options.problem.f0;
            self.g = options.problem.g;
            self.surfaceBoundary = options.problem.surfaceBoundary;
            self.zDomain = options.problem.zDomain;
            self.N2 = options.problem.N2;
            self.canonicalNativeModes = options.nativeModes;
            self.currentNativeModes = options.nativeModes;

            [self.canonicalEndpointResponseMetric, self.canonicalEnergyMatrix, self.canonicalSurfaceBuoyancyMatrix, self.canonicalBottomBuoyancyMatrix] = self.formCanonicalQuadraticMatrices();
            self.endpointResponseMetric = self.canonicalEndpointResponseMetric;
            self.energyMatrix = self.canonicalEnergyMatrix;
            self.surfaceBuoyancyMatrix = self.canonicalSurfaceBuoyancyMatrix;
            self.bottomBuoyancyMatrix = self.canonicalBottomBuoyancyMatrix;
            self.rotationName = "boundaryNormalized";
            self.rotationMatrix = repmat(eye(nEndpoints),1,1,nK);
            self.rotationEigenvalues = [];
            self.signatures = [];
            self.h0 = [];
            self.normalizationConvention = "unitEndpointResponse";
            self.rotationResiduals = struct();
            self.metadata = options.metadata;
            self.metadata.rotationName = self.rotationName;
            self.metadata.normalizationConvention = self.normalizationConvention;
        end

        function values = F(self, z)
            % Evaluate streamfunction structures $$F(z)$$.
            %
            % The result has dimensions `nZ x nEndpoints x nK`.
            %
            % - Topic: Evaluate geostrophic zero-APV modes
            % - Declaration: values = F(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: page-shaped `F` values
            arguments
                self IMGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            values = zeros(length(z),nEndpoints,nK);
            for iK = 1:nK
                values(:,:,iK) = self.solver.evaluateNativeModes(self.currentNativeModes(:,:,iK),z);
            end
        end

        function values = G(self, z)
            % Evaluate diagnostic displacement structures $$G(z)$$.
            %
            % The result has dimensions `nZ x nEndpoints x nK` and uses
            % $$G=-gN^{-2}\partial_zF$$.
            %
            % - Topic: Evaluate geostrophic zero-APV modes
            % - Declaration: values = G(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: page-shaped `G` values
            arguments
                self IMGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            N2Values = self.N2(z(:));
            N2Values = N2Values(:);
            if numel(N2Values) ~= numel(z) || any(~isfinite(N2Values)) || any(N2Values <= 0)
                error("IMGeostrophicZeroAPVModesBasis:InvalidStratification", "N2 must return one finite positive value for each z point.");
            end
            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            values = zeros(length(z),nEndpoints,nK);
            for iK = 1:nK
                dFdz = self.solver.evaluatePhysicalDerivative(self.currentNativeModes(:,:,iK),z,1);
                values(:,:,iK) = -(self.g./N2Values).*dFdz;
            end
        end

        function matrix = generalizedEnergyMatrix(self, options)
            % Return $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.
            %
            % Coefficients for endpoints absent from `endpoints` have no
            % effect because their form matrices are exactly zero.
            %
            % - Topic: Form geostrophic zero-APV quadratic forms
            % - Declaration: matrix = generalizedEnergyMatrix(basisSet,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns matrix: generalized-energy matrix pages
            arguments
                self IMGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            matrix = self.energyMatrix + options.g0*self.surfaceBuoyancyMatrix + options.gd*self.bottomBuoyancyMatrix;
            matrix = self.symmetrizePages(matrix);
        end

        function basisSet = rotateBoundaryDepth(self, options)
            % Diagonalize generalized energy relative to endpoint response.
            %
            % For each wavenumber page, solve
            %
            % $$
            % \mathsf H_g\mathbf c^a=\gamma_a\mathsf R_B\mathbf c^a,
            % \qquad
            % (\mathbf c^a)^T\mathsf R_B\mathbf c^b=\delta^{ab},
            % \qquad
            % h_0^a=2k^2\gamma_a.
            % $$
            %
            % ```matlab
            % depthModes = boundaryModes.rotateBoundaryDepth(g0=-0.035,gd=0.01);
            % ```
            %
            % - Topic: Rotate geostrophic zero-APV modes
            % - Declaration: basisSet = rotateBoundaryDepth(boundaryModes,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns basisSet: boundary-depth rotated basis
            arguments
                self IMGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            leftMatrix = self.canonicalGeneralizedEnergyMatrix(options.g0,options.gd);
            rightMatrix = self.canonicalEndpointResponseMetric;
            [modeRotation,eigenvalues,~,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,"descending");
            modeSignatures = sign(eigenvalues);
            h0Values = 2*(self.k.^2).*eigenvalues;
            basisSet = self.basisWithRotation("boundaryDepth",modeRotation,eigenvalues,modeSignatures,h0Values,"unitEndpointResponseMetric",residuals);
        end

        function basisSet = rotateSurfaceBuoyancy(self, options)
            % Diagonalize surface buoyancy relative to generalized energy.
            %
            % For every wavenumber page, solve
            %
            % $$
            % g\mathsf B_0\mathbf c^a=\chi^a\mathsf H_g\mathbf c^a.
            % $$
            %
            % The surface-buoyancy-carrying direction is ordered first.
            % This rotation leaves `h0` empty.
            %
            % ```matlab
            % surfaceModes = boundaryModes.rotateSurfaceBuoyancy(g0=-0.035,gd=0.01);
            % ```
            %
            % - Topic: Rotate geostrophic zero-APV modes
            % - Declaration: basisSet = rotateSurfaceBuoyancy(boundaryModes,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns basisSet: surface-buoyancy rotated basis
            arguments
                self IMGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            if ~ismember("surface",self.endpoints)
                error("IMGeostrophicZeroAPVModesBasis:SurfaceEndpointRequired", "rotateSurfaceBuoyancy requires a surface endpoint coordinate.");
            end
            leftMatrix = self.g*self.canonicalSurfaceBuoyancyMatrix;
            rightMatrix = self.canonicalGeneralizedEnergyMatrix(options.g0,options.gd);
            [modeRotation,eigenvalues,rightSignatures,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,"nonzeroFirst");
            basisSet = self.basisWithRotation("surfaceBuoyancy",modeRotation,eigenvalues,rightSignatures,[],"unitSignedGeneralizedEnergy",residuals);
        end

        function basisSet = rotateWithPencil(self, options)
            % Apply a custom symmetric matrix-pencil rotation.
            %
            % `leftMatrix` and `rightMatrix` are expressed in canonical
            % boundary coordinates. Two-dimensional matrices are broadcast
            % across wavenumber pages; otherwise their third dimension must
            % equal `nK`.
            % The pagewise pencil is
            %
            % $$
            % \mathsf L\mathbf c^a=\lambda_a\mathsf R\mathbf c^a,
            % \qquad
            % |(\mathbf c^a)^T\mathsf R\mathbf c^a|=1.
            % $$
            %
            % ```matlab
            % Hg = boundaryModes.generalizedEnergyMatrix(g0=-0.035,gd=0.01);
            % customModes = boundaryModes.rotateWithPencil(name="custom",leftMatrix=Hg,rightMatrix=boundaryModes.endpointResponseMetric);
            % ```
            %
            % - Topic: Rotate geostrophic zero-APV modes
            % - Declaration: basisSet = rotateWithPencil(boundaryModes,options)
            % - Parameter options.name: custom rotation name
            % - Parameter options.leftMatrix: symmetric left matrix pages
            % - Parameter options.rightMatrix: symmetric right matrix pages
            % - Returns basisSet: custom-pencil rotated basis
            arguments
                self IMGeostrophicZeroAPVModesBasis
                options.name {mustBeTextScalar}
                options.leftMatrix double {mustBeReal, mustBeFinite}
                options.rightMatrix double {mustBeReal, mustBeFinite}
            end

            name = string(options.name);
            if strlength(name) == 0
                error("IMGeostrophicZeroAPVModesBasis:InvalidRotationName", "Custom rotation name must not be empty.");
            end
            leftMatrix = self.expandMatrixPages(options.leftMatrix,"leftMatrix");
            rightMatrix = self.expandMatrixPages(options.rightMatrix,"rightMatrix");
            ordering = "descending";
            if self.pagesEqual(leftMatrix,self.g*self.canonicalSurfaceBuoyancyMatrix)
                ordering = "nonzeroFirst";
            end
            [modeRotation,eigenvalues,rightSignatures,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,ordering);
            basisSet = self.basisWithRotation(name,modeRotation,eigenvalues,rightSignatures,[],"unitSignedRightMetric",residuals);
        end

        function summarize(self)
            % Print a readable geostrophic zero-APV basis summary.
            %
            % - Topic: Inspect geostrophic zero-APV modes
            % - Declaration: summarize(basisSet)
            arguments
                self IMGeostrophicZeroAPVModesBasis
            end

            fprintf("%s\n",class(self));
            fprintf("  zDomain: [%g, %g]\n",self.zDomain(1),self.zDomain(2));
            fprintf("  endpoints: %s\n",join(self.endpoints,", "));
            fprintf("  surfaceBoundary: %s\n",self.surfaceBoundary);
            fprintf("  rotationName: %s\n",self.rotationName);
            fprintf("  nWavenumbers: %d\n",numel(self.k));
            fprintf("  solver: %s\n",class(self.solver));
        end
    end

    methods (Hidden)
        function tf = sharesCanonicalBasisWith(self, other)
            arguments
                self IMGeostrophicZeroAPVModesBasis
                other IMGeostrophicZeroAPVModesBasis
            end

            tf = isequal(self.canonicalNativeModes,other.canonicalNativeModes) ...
                && isequal(self.k,other.k) ...
                && isequal(self.endpoints,other.endpoints) ...
                && isequal(self.zDomain,other.zDomain) ...
                && isequal(self.f0,other.f0) ...
                && isequal(self.g,other.g) ...
                && isequal(self.surfaceBoundary,other.surfaceBoundary);
        end
    end

    methods (Access = private)
        function [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy] = formCanonicalQuadraticMatrices(self)
            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            responseMetric = zeros(nEndpoints,nEndpoints,nK);
            physicalEnergy = zeros(nEndpoints,nEndpoints,nK);
            surfaceBuoyancy = zeros(nEndpoints,nEndpoints,nK);
            bottomBuoyancy = zeros(nEndpoints,nEndpoints,nK);
            zSurface = self.zDomain(2);
            zBottom = self.zDomain(1);
            endpointRows = zeros(1,nEndpoints);
            endpointRows(self.endpoints == "surface") = 1;
            endpointRows(self.endpoints == "bottom") = 2;

            for iK = 1:nK
                nativeModes = self.canonicalNativeModes(:,:,iK);
                FSurface = self.solver.evaluateNativeModes(nativeModes,zSurface);
                FBottom = self.solver.evaluateNativeModes(nativeModes,zBottom);
                dFSurface = self.solver.evaluatePhysicalDerivative(nativeModes,zSurface,1);
                dFBottom = self.solver.evaluatePhysicalDerivative(nativeModes,zBottom,1);
                N2Surface = self.N2(zSurface);
                N2Bottom = self.N2(zBottom);
                GSurface = -(self.g/N2Surface)*dFSurface;
                GBottom = -(self.g/N2Bottom)*dFBottom;
                surfaceResponse = GSurface;
                if self.surfaceBoundary == "freeSurface"
                    surfaceResponse = surfaceResponse - FSurface;
                end
                bottomResponse = GBottom;
                fullResponse = [surfaceResponse; bottomResponse];
                selectedResponse = fullResponse(endpointRows,:);

                responseMetric(:,:,iK) = self.symmetrize(selectedResponse.'*selectedResponse);
                energyScale = self.f0^2/(2*self.g*self.k(iK)^4);
                energyBoundaryForm = -FSurface.'*surfaceResponse + FBottom.'*bottomResponse;
                physicalEnergy(:,:,iK) = energyScale*self.symmetrize(energyBoundaryForm);
                buoyancyScale = self.f0^2/(2*self.g^2*self.k(iK)^4);
                if ismember("surface",self.endpoints)
                    surfaceBuoyancy(:,:,iK) = buoyancyScale*self.symmetrize(surfaceResponse.'*surfaceResponse);
                end
                if ismember("bottom",self.endpoints)
                    bottomBuoyancy(:,:,iK) = buoyancyScale*self.symmetrize(bottomResponse.'*bottomResponse);
                end
            end
        end

        function matrix = canonicalGeneralizedEnergyMatrix(self,g0,gd)
            matrix = self.canonicalEnergyMatrix + g0*self.canonicalSurfaceBuoyancyMatrix + gd*self.canonicalBottomBuoyancyMatrix;
            matrix = self.symmetrizePages(matrix);
        end

        function basisSet = basisWithRotation(self,name,rotationMatrix,eigenvalues,signatures,h0Values,normalizationConvention,residuals)
            IMGeostrophicZeroAPVFormTools.validateRotationPages(rotationMatrix,numel(self.endpoints),numel(self.k),"IMGeostrophicZeroAPVModesBasis");
            basisSet = IMGeostrophicZeroAPVModesBasis(problem=self.problem,solver=self.solver,nativeModes=self.canonicalNativeModes,metadata=self.metadata);
            nK = numel(self.k);
            for iK = 1:nK
                C = rotationMatrix(:,:,iK);
                basisSet.currentNativeModes(:,:,iK) = self.canonicalNativeModes(:,:,iK)*C;
                basisSet.endpointResponseMetric(:,:,iK) = self.symmetrize(C.'*self.canonicalEndpointResponseMetric(:,:,iK)*C);
                basisSet.energyMatrix(:,:,iK) = self.symmetrize(C.'*self.canonicalEnergyMatrix(:,:,iK)*C);
                basisSet.surfaceBuoyancyMatrix(:,:,iK) = self.symmetrize(C.'*self.canonicalSurfaceBuoyancyMatrix(:,:,iK)*C);
                basisSet.bottomBuoyancyMatrix(:,:,iK) = self.symmetrize(C.'*self.canonicalBottomBuoyancyMatrix(:,:,iK)*C);
            end
            basisSet.rotationName = string(name);
            basisSet.rotationMatrix = rotationMatrix;
            basisSet.rotationEigenvalues = eigenvalues;
            basisSet.signatures = signatures;
            basisSet.h0 = h0Values;
            basisSet.normalizationConvention = string(normalizationConvention);
            basisSet.rotationResiduals = residuals;
            basisSet.metadata.rotationName = basisSet.rotationName;
            basisSet.metadata.normalizationConvention = basisSet.normalizationConvention;
            basisSet.metadata.rotationEigenvalues = eigenvalues;
            basisSet.metadata.signatures = signatures;
        end

        function pages = expandMatrixPages(self,matrix,name)
            pages = IMGeostrophicZeroAPVFormTools.expandMatrixPages(matrix,numel(self.endpoints),numel(self.k),name,"IMGeostrophicZeroAPVModesBasis");
        end

        function [vectors,eigenvalues,rightSignatures,residuals] = solveSymmetricPencil(self,leftMatrix,rightMatrix,ordering)
            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            [vectors,eigenvalues,rightSignatures,residuals] = IMGeostrophicZeroAPVFormTools.solveSymmetricPencil(leftMatrix,rightMatrix,nEndpoints,nK,ordering,"IMGeostrophicZeroAPVModesBasis");
        end

        function tf = pagesEqual(~,left,right)
            tf = IMGeostrophicZeroAPVFormTools.pagesEqual(left,right);
        end
    end

    methods (Static, Access = private)
        function matrix = symmetrize(matrix)
            matrix = IMGeostrophicZeroAPVFormTools.symmetrize(matrix);
        end

        function pages = symmetrizePages(pages)
            pages = IMGeostrophicZeroAPVFormTools.symmetrizePages(pages);
        end
    end
end
