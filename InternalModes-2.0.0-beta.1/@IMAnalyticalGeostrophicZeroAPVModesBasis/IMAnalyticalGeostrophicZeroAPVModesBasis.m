classdef IMAnalyticalGeostrophicZeroAPVModesBasis
    % Store exact canonical or rotated geostrophic zero-APV modes.
    %
    % Concrete analytical stratification families create the canonical
    % columns. Their endpoint responses satisfy
    %
    % $$
    % \mathbf B[F_0^{\mathrm{sur}}]=(1,0)^T,
    % \qquad
    % \mathbf B[F_0^{\mathrm{bot}}]=(0,1)^T,
    % $$
    %
    % for the requested subset of surface and bottom endpoints. The exact
    % diagnostic variable is
    %
    % $$
    % G(z)=-\frac{g}{N^2(z)}\frac{\partial F}{\partial z}(z).
    % $$
    %
    % `F(z)` and `G(z)` have dimensions `nZ x nEndpoints x nK`. Quadratic
    % forms have dimensions `nEndpoints x nEndpoints x nK`. Rotations apply
    % the same pagewise column map to both exact variables.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-4000 0],f0=1e-4);
    % exactModes = solution.geostrophicZeroAPVModesAtWavenumber(1e-4);
    % depthModes = exactModes.rotateBoundaryDepth(g0=-0.03,gd=0.01);
    % ```
    %
    % - Topic: Evaluate exact geostrophic zero-APV modes
    % - Topic: Inspect exact geostrophic zero-APV modes
    % - Topic: Form exact geostrophic zero-APV quadratic forms
    % - Topic: Rotate exact geostrophic zero-APV modes
    % - Topic: Developer topics
    % - Declaration: classdef IMAnalyticalGeostrophicZeroAPVModesBasis

    properties (SetAccess = private)
        % Analytical solution family that created this basis.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        solution

        % Canonical geostrophic zero-APV problem represented by the formulas.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        problem

        % Horizontal-wavenumber pages.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        k

        % Canonically ordered endpoint-coordinate labels.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        endpoints

        % Endpoint-response Gram matrix $$\mathsf R_B$$.
        %
        % Canonical unit endpoint coordinates give $$\mathsf R_B=\mathsf I$$.
        %
        % - Topic: Form exact geostrophic zero-APV quadratic forms
        endpointResponseMetric

        % Physical-energy matrix $$\mathsf H$$.
        %
        % - Topic: Form exact geostrophic zero-APV quadratic forms
        energyMatrix

        % Surface-buoyancy matrix $$\mathsf B_0$$.
        %
        % This matrix is exactly zero when the surface endpoint is absent.
        %
        % - Topic: Form exact geostrophic zero-APV quadratic forms
        surfaceBuoyancyMatrix

        % Bottom-buoyancy matrix $$\mathsf B_d$$.
        %
        % This matrix is exactly zero when the bottom endpoint is absent.
        %
        % - Topic: Form exact geostrophic zero-APV quadratic forms
        bottomBuoyancyMatrix

        % Name of the current boundary-coordinate rotation.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        rotationName

        % Canonical-to-current pagewise column transformation.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        rotationMatrix

        % Pagewise generalized-pencil eigenvalues.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        rotationEigenvalues

        % Pagewise quadratic-form signatures.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        signatures

        % Signed boundary depths from the boundary-depth rotation.
        %
        % Only `rotateBoundaryDepth` populates `h0`, using
        % $$h_0^a=2k^2\gamma_a$$.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        h0

        % Rotation normalization convention.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        normalizationConvention

        % Pagewise pencil and normalization residuals.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        rotationResiduals

        % Coriolis parameter.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        f0

        % Gravitational acceleration.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        g

        % Surface endpoint convention.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        surfaceBoundary

        % Physical vertical domain.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        N2

        % Additional creation and rotation metadata.
        %
        % - Topic: Inspect exact geostrophic zero-APV modes
        metadata
    end

    properties (Access = private)
        canonicalFFunction
        canonicalGFunction
        canonicalEndpointResponseMetric
        canonicalEnergyMatrix
        canonicalSurfaceBuoyancyMatrix
        canonicalBottomBuoyancyMatrix
    end

    methods
        function self = IMAnalyticalGeostrophicZeroAPVModesBasis(options)
            % Create an exact canonical boundary-normalized basis.
            %
            % Concrete analytical solution families supply exact `F` and
            % `G` evaluators. Users normally construct this basis through
            % `geostrophicZeroAPVModesAtWavenumber` on one of those families.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(options)
            % - Parameter options.solution: analytical solution family
            % - Parameter options.problem: canonical zero-APV problem
            % - Parameter options.FFunction: exact canonical `F` evaluator
            % - Parameter options.GFunction: exact canonical `G` evaluator
            % - Parameter options.metadata: creation metadata
            % - Returns basisSet: exact canonical basis
            % - Developer: true
            arguments
                options.solution IMAnalyticalSolution
                options.problem IMGeostrophicZeroAPVModes
                options.FFunction (1,1) function_handle
                options.GFunction (1,1) function_handle
                options.metadata struct = struct()
            end

            self.solution = options.solution;
            self.problem = options.problem;
            self.k = options.problem.k;
            self.endpoints = options.problem.endpoints;
            self.f0 = options.problem.f0;
            self.g = options.problem.g;
            self.surfaceBoundary = options.problem.surfaceBoundary;
            self.zDomain = options.problem.zDomain;
            self.N2 = options.problem.N2;
            self.canonicalFFunction = options.FFunction;
            self.canonicalGFunction = options.GFunction;

            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            self.rotationName = "boundaryNormalized";
            self.rotationMatrix = repmat(eye(nEndpoints),1,1,nK);
            self.rotationEigenvalues = [];
            self.signatures = [];
            self.h0 = [];
            self.normalizationConvention = "unitEndpointResponse";
            self.rotationResiduals = struct();

            [self.canonicalEndpointResponseMetric,self.canonicalEnergyMatrix,self.canonicalSurfaceBuoyancyMatrix,self.canonicalBottomBuoyancyMatrix,responseResidual] = self.formCanonicalQuadraticMatrices();
            self.endpointResponseMetric = self.canonicalEndpointResponseMetric;
            self.energyMatrix = self.canonicalEnergyMatrix;
            self.surfaceBuoyancyMatrix = self.canonicalSurfaceBuoyancyMatrix;
            self.bottomBuoyancyMatrix = self.canonicalBottomBuoyancyMatrix;
            self.metadata = options.metadata;
            self.metadata.rotationName = self.rotationName;
            self.metadata.normalizationConvention = self.normalizationConvention;
            self.metadata.canonicalEndpointResponseResidual = responseResidual;
        end

        function values = F(self,z)
            % Evaluate exact streamfunction structures $$F(z)$$.
            %
            % The result has dimensions `nZ x nEndpoints x nK`.
            %
            % - Topic: Evaluate exact geostrophic zero-APV modes
            % - Declaration: values = F(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: page-shaped exact `F` values
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.evaluateCanonical(self.canonicalFFunction,z);
            values = self.applyRotation(values);
        end

        function values = G(self,z)
            % Evaluate exact diagnostic displacement structures $$G(z)$$.
            %
            % The evaluator uses the closed-form derivative relation
            % $$G=-gN^{-2}\partial_zF$$ without numerical differentiation.
            %
            % - Topic: Evaluate exact geostrophic zero-APV modes
            % - Declaration: values = G(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: page-shaped exact `G` values
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.evaluateCanonical(self.canonicalGFunction,z);
            values = self.applyRotation(values);
        end

        function matrix = generalizedEnergyMatrix(self,options)
            % Return $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.
            %
            % A coefficient for an absent endpoint has no effect because its
            % form matrix is exactly zero.
            %
            % - Topic: Form exact geostrophic zero-APV quadratic forms
            % - Declaration: matrix = generalizedEnergyMatrix(basisSet,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns matrix: generalized-energy matrix pages
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            matrix = self.energyMatrix + options.g0*self.surfaceBuoyancyMatrix + options.gd*self.bottomBuoyancyMatrix;
            matrix = IMGeostrophicZeroAPVFormTools.symmetrizePages(matrix);
        end

        function basisSet = rotateBoundaryDepth(self,options)
            % Diagonalize generalized energy relative to endpoint response.
            %
            % For each wavenumber page,
            %
            % $$
            % \mathsf H_g\mathbf c^a=\gamma_a\mathsf R_B\mathbf c^a,
            % \qquad h_0^a=2k^2\gamma_a.
            % $$
            %
            % - Topic: Rotate exact geostrophic zero-APV modes
            % - Declaration: basisSet = rotateBoundaryDepth(exactModes,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns basisSet: boundary-depth rotated exact basis
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            leftMatrix = self.canonicalGeneralizedEnergyMatrix(options.g0,options.gd);
            rightMatrix = self.canonicalEndpointResponseMetric;
            [rotation,eigenvalues,~,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,"descending");
            h0Values = 2*(self.k.^2).*eigenvalues;
            basisSet = self.basisWithRotation("boundaryDepth",rotation,eigenvalues,sign(eigenvalues),h0Values,"unitEndpointResponseMetric",residuals);
        end

        function basisSet = rotateSurfaceBuoyancy(self,options)
            % Diagonalize surface buoyancy relative to generalized energy.
            %
            % The surface-carrying direction is ordered first. This rotation
            % reports its eigenvalues through `rotationEigenvalues` and leaves
            % `h0` empty.
            %
            % - Topic: Rotate exact geostrophic zero-APV modes
            % - Declaration: basisSet = rotateSurfaceBuoyancy(exactModes,options)
            % - Parameter options.g0: finite signed surface coefficient
            % - Parameter options.gd: finite signed bottom coefficient
            % - Returns basisSet: surface-buoyancy rotated exact basis
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.gd (1,1) double {mustBeReal, mustBeFinite} = 0
            end

            if ~ismember("surface",self.endpoints)
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:SurfaceEndpointRequired", "rotateSurfaceBuoyancy requires a surface endpoint coordinate.");
            end
            leftMatrix = self.g*self.canonicalSurfaceBuoyancyMatrix;
            rightMatrix = self.canonicalGeneralizedEnergyMatrix(options.g0,options.gd);
            [rotation,eigenvalues,rightSignatures,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,"nonzeroFirst");
            basisSet = self.basisWithRotation("surfaceBuoyancy",rotation,eigenvalues,rightSignatures,[],"unitSignedGeneralizedEnergy",residuals);
        end

        function basisSet = rotateWithPencil(self,options)
            % Apply a custom symmetric matrix-pencil rotation.
            %
            % Input matrices are expressed in canonical endpoint coordinates.
            % Two-dimensional matrices broadcast across wavenumber pages.
            %
            % - Topic: Rotate exact geostrophic zero-APV modes
            % - Declaration: basisSet = rotateWithPencil(exactModes,options)
            % - Parameter options.name: custom rotation name
            % - Parameter options.leftMatrix: symmetric left matrix pages
            % - Parameter options.rightMatrix: symmetric right matrix pages
            % - Returns basisSet: custom-pencil rotated exact basis
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.name {mustBeTextScalar}
                options.leftMatrix double {mustBeReal, mustBeFinite}
                options.rightMatrix double {mustBeReal, mustBeFinite}
            end

            name = string(options.name);
            if strlength(name) == 0
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidRotationName", "Custom rotation name must not be empty.");
            end
            leftMatrix = self.expandMatrixPages(options.leftMatrix,"leftMatrix");
            rightMatrix = self.expandMatrixPages(options.rightMatrix,"rightMatrix");
            ordering = "descending";
            if IMGeostrophicZeroAPVFormTools.pagesEqual(leftMatrix,self.g*self.canonicalSurfaceBuoyancyMatrix)
                ordering = "nonzeroFirst";
            end
            [rotation,eigenvalues,rightSignatures,residuals] = self.solveSymmetricPencil(leftMatrix,rightMatrix,ordering);
            basisSet = self.basisWithRotation(name,rotation,eigenvalues,rightSignatures,[],"unitSignedRightMetric",residuals);
        end

        function summarize(self)
            % Print a readable exact zero-APV basis summary.
            %
            % - Topic: Inspect exact geostrophic zero-APV modes
            % - Declaration: summarize(exactModes)
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
            end

            fprintf("%s\n",class(self));
            fprintf("  solution: %s\n",class(self.solution));
            fprintf("  zDomain: [%g, %g]\n",self.zDomain(1),self.zDomain(2));
            fprintf("  endpoints: %s\n",join(self.endpoints,", "));
            fprintf("  surfaceBoundary: %s\n",self.surfaceBoundary);
            fprintf("  rotationName: %s\n",self.rotationName);
            fprintf("  nWavenumbers: %d\n",numel(self.k));
        end
    end

    methods (Access = private)
        function values = evaluateCanonical(self,evaluator,z)
            values = evaluator(z(:));
            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            if size(values,1) ~= numel(z) || size(values,2) ~= nEndpoints || size(values,3) ~= nK
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidEvaluatorShape", "Exact evaluators must return nZ x nEndpoints x nK arrays.");
            end
            if ~isreal(values) || any(~isfinite(values(:)))
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidEvaluatorValues", "Exact evaluators must return finite real values.");
            end
        end

        function values = applyRotation(self,values)
            for iK = 1:numel(self.k)
                values(:,:,iK) = values(:,:,iK)*self.rotationMatrix(:,:,iK);
            end
        end

        function [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy,responseResidual] = formCanonicalQuadraticMatrices(self)
            zSurface = self.zDomain(2);
            zBottom = self.zDomain(1);
            FSurface = self.evaluateCanonical(self.canonicalFFunction,zSurface);
            FBottom = self.evaluateCanonical(self.canonicalFFunction,zBottom);
            GSurface = self.evaluateCanonical(self.canonicalGFunction,zSurface);
            GBottom = self.evaluateCanonical(self.canonicalGFunction,zBottom);
            surfaceResponse = GSurface;
            if self.surfaceBoundary == "freeSurface"
                surfaceResponse = surfaceResponse-FSurface;
            end
            bottomResponse = GBottom;
            fullResponse = cat(1,surfaceResponse,bottomResponse);
            endpointRows = zeros(1,numel(self.endpoints));
            endpointRows(self.endpoints == "surface") = 1;
            endpointRows(self.endpoints == "bottom") = 2;
            selectedResponse = fullResponse(endpointRows,:,:);
            responseResidual = zeros(1,numel(self.k));
            for iK = 1:numel(self.k)
                responseResidual(iK) = norm(selectedResponse(:,:,iK)-eye(numel(self.endpoints)),"fro")/max(1,sqrt(numel(self.endpoints)));
            end
            if any(responseResidual > 1e-10)
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:EndpointNormalizationFailure", "Exact canonical modes did not satisfy unit endpoint responses.");
            end
            [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy] = IMGeostrophicZeroAPVFormTools.formMatrices(FSurface,FBottom,surfaceResponse,bottomResponse,selectedResponse,self.endpoints,self.k,self.f0,self.g);
        end

        function matrix = canonicalGeneralizedEnergyMatrix(self,g0,gd)
            matrix = self.canonicalEnergyMatrix + g0*self.canonicalSurfaceBuoyancyMatrix + gd*self.canonicalBottomBuoyancyMatrix;
            matrix = IMGeostrophicZeroAPVFormTools.symmetrizePages(matrix);
        end

        function basisSet = basisWithRotation(self,name,rotation,eigenvalues,signatures,h0Values,normalizationConvention,residuals)
            IMGeostrophicZeroAPVFormTools.validateRotationPages(rotation,numel(self.endpoints),numel(self.k),"IMAnalyticalGeostrophicZeroAPVModesBasis");
            basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(solution=self.solution,problem=self.problem,FFunction=self.canonicalFFunction,GFunction=self.canonicalGFunction,metadata=self.metadata);
            for iK = 1:numel(self.k)
                C = rotation(:,:,iK);
                basisSet.endpointResponseMetric(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(C.'*self.canonicalEndpointResponseMetric(:,:,iK)*C);
                basisSet.energyMatrix(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(C.'*self.canonicalEnergyMatrix(:,:,iK)*C);
                basisSet.surfaceBuoyancyMatrix(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(C.'*self.canonicalSurfaceBuoyancyMatrix(:,:,iK)*C);
                basisSet.bottomBuoyancyMatrix(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(C.'*self.canonicalBottomBuoyancyMatrix(:,:,iK)*C);
            end
            basisSet.rotationName = string(name);
            basisSet.rotationMatrix = rotation;
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
            pages = IMGeostrophicZeroAPVFormTools.expandMatrixPages(matrix,numel(self.endpoints),numel(self.k),name,"IMAnalyticalGeostrophicZeroAPVModesBasis");
        end

        function [vectors,eigenvalues,rightSignatures,residuals] = solveSymmetricPencil(self,leftMatrix,rightMatrix,ordering)
            [vectors,eigenvalues,rightSignatures,residuals] = IMGeostrophicZeroAPVFormTools.solveSymmetricPencil(leftMatrix,rightMatrix,numel(self.endpoints),numel(self.k),ordering,"IMAnalyticalGeostrophicZeroAPVModesBasis");
        end
    end
end
