classdef IMDiscreteTransform
    % Store forward and inverse matrices for a scalar modal transform.
    %
    % Let $$n_z$$ be the number of sample points, $$n_m$$ the number of
    % retained modes, and $$n_p$$ the number of sampled profiles. The
    % sampled modal basis is
    %
    % $$
    % A_{\mathrm i}=\Phi\in\mathbb{R}^{n_z\times n_m},
    % \qquad
    % \Phi_{ij}=u_j(z_i),
    % $$
    %
    % where the columns contain the retained, normalized modes. For the
    % sample-space metric $$W$$, the sampled Gram matrix and Galerkin
    % forward matrix are
    %
    % $$
    % \Gamma=\Phi^\mathsf{T}W\Phi,
    % \qquad
    % A_{\mathrm f}=\Gamma^{-1}\Phi^\mathsf{T}W
    % =\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
    % A_{\mathrm i}^\mathsf{T}W.
    % $$
    %
    % For sampled profiles $$X\in\mathbb{R}^{n_z\times n_p}$$ and modal
    % coefficients $$A\in\mathbb{R}^{n_m\times n_p}$$, the forward and
    % back transforms are
    %
    % $$
    % A=A_{\mathrm f}X,
    % \qquad
    % \widehat{X}=A_{\mathrm i}A.
    % $$
    %
    % The matrices are generally rectangular. For a full-rank retained
    % basis they obey
    %
    % $$
    % A_{\mathrm f}A_{\mathrm i}=I,
    % \qquad
    % A_{\mathrm i}A_{\mathrm f}=P_W,
    % $$
    %
    % where $$P_W$$ is the sampled-space Galerkin projector onto the
    % retained modal subspace. Thus transforming retained coefficients back
    % and then forward recovers those coefficients exactly, while a general
    % sampled profile is projected rather than necessarily reproduced. When
    % $$W$$ is positive definite, $$P_W$$ is the $$W$$-orthogonal projector.
    %
    % Transform construction and quadrature accuracy are related but
    % distinct. If $$\Gamma_0$$ is the continuous target Gram matrix, define
    %
    % $$
    % S=\operatorname{diag}\!\left(
    % \left|\operatorname{diag}\Gamma_0\right|^{-1/2}
    % \right),
    % \qquad
    % E=S(\Gamma-\Gamma_0)S.
    % $$
    %
    % `relativeGramOperatorError` is $$\|E\|_2$$: the largest Gram
    % distortion over any normalized combination of retained modes. By
    % contrast, `roundTripError` measures only how accurately the two
    % transform matrices recover retained coefficients. It can be near
    % roundoff even when the quadrature reproduces $$\Gamma_0$$ poorly.
    %
    % Construct transforms from solved scalar modes with
    % `IMBasisSet.discreteTransform`.
    %
    % ```matlab
    % transform = basisSet.discreteTransform(z=z,nModes=8);
    % coefficients = transform.transformForward(values);
    % valuesFit = transform.transformBack(coefficients);
    % ```
    %
    % - Topic: Create discrete transforms
    % - Topic: Use transform matrices
    % - Topic: Apply discrete transforms
    % - Topic: Inspect samples and modes
    % - Topic: Assess transform quality
    % - Topic: Developer topics
    % - Declaration: classdef IMDiscreteTransform

    properties (SetAccess = private)
        % Map sampled profiles to retained modal coefficients.
        %
        % The `forwardMatrix` is the $$n_m\times n_z$$ Galerkin matrix
        %
        % $$
        % A_{\mathrm f}
        % =\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
        % A_{\mathrm i}^\mathsf{T}W
        % =\Gamma^{-1}\Phi^\mathsf{T}W,
        % $$
        %
        % where $$A_{\mathrm i}=\Phi$$ is `inverseMatrix`, $$W$$ is
        % `metricMatrix`, and $$\Gamma=\Phi^\mathsf{T}W\Phi$$ is
        % `gramMatrix`. For sampled profiles $$X$$, the coefficients
        % $$A=A_{\mathrm f}X$$ satisfy the Galerkin normal equations
        %
        % $$
        % \Phi^\mathsf{T}W\left(X-\Phi A\right)=0.
        % $$
        %
        % Direct multiplication and `transformForward` are equivalent:
        %
        % ```matlab
        % coefficientsByMatrix = transform.forwardMatrix*values;
        % coefficients = transform.transformForward(values);
        % ```
        %
        % - Topic: Use transform matrices
        % - nav_order: 10
        forwardMatrix

        % Map retained modal coefficients back to sampled profiles.
        %
        % The `inverseMatrix` is the $$n_z\times n_m$$ sampled modal basis
        %
        % $$
        % A_{\mathrm i}=\Phi,
        % \qquad
        % (A_{\mathrm i})_{ij}=\Phi_{ij}=u_j(z_i).
        % $$
        %
        % Row $$i$$ corresponds to `z(i)`, and column $$j$$ corresponds to
        % `modeNumber(j)`. The sampled modes use the normalization recorded
        % by `normalization`. Direct multiplication and `transformBack` are
        % equivalent:
        %
        % ```matlab
        % valuesByMatrix = transform.inverseMatrix*coefficients;
        % values = transform.transformBack(coefficients);
        % ```
        %
        % - Topic: Use transform matrices
        % - nav_order: 20
        inverseMatrix

        % Physical points at which profiles and modes are sampled.
        %
        % `z` is an $$n_z\times1$$ column vector. Its entries define the row
        % ordering of `inverseMatrix`, the column ordering of
        % `forwardMatrix`, and the expected row ordering of values passed to
        % `transformForward`. Coordinates use the same units as the source
        % basis set, normally meters for vertical-mode problems.
        %
        % ```matlab
        % z = transform.z;
        % values = profile(z);
        % coefficients = transform.transformForward(values);
        % ```
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 10
        z

        % Quadrature weights associated with the sample points.
        %
        % `weights(i)` is the quadrature weight $$w_i$$ associated with
        % `z(i)`. For transforms built by `IMBasisSet`, the interior part of
        % the sampled metric begins with
        %
        % $$
        % W_{\mathrm{int}}
        % =\operatorname{diag}\!\left(r(z_i)w_i\right),
        % $$
        %
        % before supported endpoint terms are added. The weights may be
        % geometric control-volume widths or fitted algebraic weights.
        % `hasNegativeWeights` reports whether any are negative.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 20
        weights

        % Physical labels for the retained modal rows and columns.
        %
        % `modeNumber` is a $$1\times n_m$$ row vector. Entry $$j$$ labels
        % column $$j$$ of `inverseMatrix`, row $$j$$ of `forwardMatrix`, and
        % row $$j$$ of arrays returned by `transformForward`. These are mode
        % labels, not MATLAB array indices, and may include negative or zero
        % values when the source basis contains such modes.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 30
        modeNumber

        % Name of the normalization captured by this transform.
        %
        % The columns of `inverseMatrix` were sampled using this basis-set
        % normalization, so modal coefficients are defined relative to the
        % same scaling. The value is a snapshot taken when the transform was
        % built; subsequently changing the source basis normalization does
        % not modify an existing transform.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 40
        normalization

        % Measure the worst normalized Gram distortion.
        %
        % Let $$\Gamma$$ be `gramMatrix`, let $$\Gamma_0$$ be
        % `targetGramMatrix`, and define
        %
        % $$
        % S=\operatorname{diag}\!\left(
        % \left|\operatorname{diag}\Gamma_0\right|^{-1/2}
        % \right).
        % $$
        %
        % With
        %
        % $$
        % E=S(\Gamma-\Gamma_0)S,
        % $$
        %
        % the reported error is
        %
        % $$
        % \left\|S(\Gamma-\Gamma_0)S\right\|_2.
        % $$
        %
        % Zero means the sampled metric reproduces every retained modal
        % inner product exactly. The operator norm asks for the single
        % normalized combination of retained modes whose Gram value is most
        % distorted. This differs from the default quadrature-fit
        % `residualNorm`, which uses the Frobenius norm $$\|E\|_{\mathrm F}$$
        % and aggregates errors over all individual mode pairs. When
        % `targetGramIsPositiveDefinite` is true, this is the worst-case
        % relative quadratic-form, or Parseval, error. For a signed target
        % it remains a useful magnitude-scaled discrepancy, but not a
        % positive-norm error.
        %
        % ```matlab
        % operatorError = transform.relativeGramOperatorError;
        % ```
        %
        % - Topic: Assess transform quality
        % - nav_order: 10
        relativeGramOperatorError

        % Measure recovery of retained modal coefficients.
        %
        % The reported value is
        %
        % $$
        % \left\|A_{\mathrm f}A_{\mathrm i}-I_{n_m}\right\|_2.
        % $$
        %
        % A small value means coefficients transformed back to sample space
        % and then forward are recovered accurately. This is an algebraic
        % consistency check on the transform matrices, not a quadrature
        % accuracy metric. It can be near roundoff even when
        % `relativeGramOperatorError` is appreciable. It also does not
        % measure the reconstruction error of an arbitrary sampled profile;
        % that profile is generally projected by
        % $$A_{\mathrm i}A_{\mathrm f}$$.
        %
        % - Topic: Assess transform quality
        % - nav_order: 20
        roundTripError

        % Two-norm condition number of the sampled modal basis.
        %
        % This is
        %
        % $$
        % \kappa_2(A_{\mathrm i})
        % =\frac{\sigma_{\max}(A_{\mathrm i})}
        % {\sigma_{\min}(A_{\mathrm i})}.
        % $$
        %
        % Large values indicate that retained mode columns are nearly
        % linearly dependent on the selected sample points. The definition
        % applies to rectangular `inverseMatrix` matrices through their
        % singular values.
        %
        % - Topic: Assess transform quality
        % - nav_order: 30
        inverseMatrixConditionNumber

        % Numerical rank of the sampled modal Gram matrix.
        %
        % `sampledGramRank` is the number of singular values of
        % $$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$ above the
        % scale-aware numerical rank tolerance. A value smaller than the
        % number of retained modes means the sampled rule cannot distinguish
        % the complete modal family. In that case `forwardMatrix` is formed
        % with a pseudoinverse and `relativeGramOperatorError` is `Inf` so an
        % automatic retained-band policy can reject the deficient prefix.
        %
        % - Topic: Assess transform quality
        % - nav_order: 35
        sampledGramRank

        % Two-norm condition number of the sampled Gram matrix.
        %
        % This is $$\kappa_2(\Gamma)$$ for
        % $$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$. Large values
        % indicate that the Galerkin normal equations used to construct
        % `forwardMatrix` are sensitive to perturbations. This diagnostic
        % depends on both the sampled modes and the metric, unlike
        % `inverseMatrixConditionNumber`.
        %
        % - Topic: Assess transform quality
        % - nav_order: 40
        gramConditionNumber

        % Whether the target modal Gram matrix defines a positive norm.
        %
        % `targetGramMatrix` is required to be diagonal, so this property is
        % true exactly when every target modal norm
        % $$C_j=(\Gamma_0)_{jj}$$ is positive. A false value does not make
        % the transform invalid: canonical endpoint weights can produce a
        % signed metric. It changes the interpretation of
        % `relativeGramOperatorError` from a relative norm error to a signed,
        % magnitude-scaled Gram discrepancy.
        %
        % - Topic: Assess transform quality
        % - nav_order: 50
        targetGramIsPositiveDefinite

        % Whether at least one quadrature weight is negative.
        %
        % This is `true` when any $$w_i<0$$. A negative weight means the
        % quadrature rule is algebraic rather than a positive weighted sum,
        % but it does not by itself invalidate the transform. Inspect
        % `targetGramIsPositiveDefinite` and the Gram diagnostics to assess
        % the resulting metric.
        %
        % - Topic: Assess transform quality
        % - nav_order: 60
        hasNegativeWeights

        % Sample-space bilinear-form matrix $$W$$.
        %
        % The metric defines the sampled bilinear form
        %
        % $$
        % \langle x,y\rangle_W=x^\mathsf{T}Wy.
        % $$
        %
        % For transforms built by `IMBasisSet`, its structure is
        %
        % $$
        % W=\operatorname{diag}\!\left(r(z_i)w_i\right)
        % +W_{\mathrm{endpoint}},
        % $$
        %
        % where supported value-only endpoint terms are represented in
        % $$W_{\mathrm{endpoint}}$$. The matrix is $$n_z\times n_z$$,
        % symmetric, and may be indefinite.
        %
        % ```matlab
        % metricSymmetryError = norm(transform.metricMatrix-transform.metricMatrix.',2);
        % ```
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 10
        % - Developer: true
        metricMatrix

        % Sampled modal Gram matrix $$\Gamma$$.
        %
        % The transform computes
        %
        % $$
        % \Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}
        % =\Phi^\mathsf{T}W\Phi.
        % $$
        %
        % It is an $$n_m\times n_m$$ matrix of sampled inner products among
        % the retained modes. It enters the full definition
        % $$A_{\mathrm f}=\Gamma^{-1}\Phi^\mathsf{T}W$$ and is compared with
        % `targetGramMatrix` by `relativeGramOperatorError`.
        %
        % ```matlab
        % gramDifference = transform.gramMatrix-transform.targetGramMatrix;
        % ```
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 20
        % - Developer: true
        gramMatrix

        % Continuous diagonal Gram matrix targeted by the quadrature rule.
        %
        % The matrix
        %
        % $$
        % \Gamma_0=\operatorname{diag}(C_1,\ldots,C_{n_m})
        % $$
        %
        % contains the continuous full-domain inner products of the retained
        % normalized modes. A fitted quadrature rule seeks to make
        % `gramMatrix` reproduce this target. The entries $$C_j$$ are finite
        % and nonzero, but may be negative for a signed canonical metric.
        %
        % ```matlab
        % targetNorms = diag(transform.targetGramMatrix);
        % ```
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 30
        % - Developer: true
        targetGramMatrix
    end

    methods
        function self = IMDiscreteTransform(options)
            % Create a scalar discrete transform from canonical matrices.
            %
            % Let $$n_z$$ be the number of samples and $$n_m$$ the number of
            % retained modes. `inverseMatrix` must be $$n_z\times n_m$$;
            % `z` and `weights` must each contain $$n_z$$ entries;
            % `metricMatrix` must be a symmetric $$n_z\times n_z$$ matrix;
            % and `targetGramMatrix` must be a diagonal $$n_m\times n_m$$
            % matrix with finite, nonzero diagonal entries. `modeNumber`
            % supplies one label for each retained mode column.
            %
            % The constructor derives `gramMatrix`, `forwardMatrix`, and all
            % quality diagnostics from these inputs. Most users should build
            % this object from a solved basis with
            % `IMBasisSet.discreteTransform`; direct construction is useful
            % for alternative transform builders.
            %
            % ```matlab
            % z = [-1; -0.5; 0];
            % weights = [0.25; 0.5; 0.25];
            % inverseMatrix = [1 0; 1 1; 0 1];
            % metricMatrix = diag(weights);
            % transform = IMDiscreteTransform(z=z,weights=weights,modeNumber=[1 2],normalization="unity",inverseMatrix=inverseMatrix,metricMatrix=metricMatrix,targetGramMatrix=eye(2));
            % ```
            %
            % - Topic: Create discrete transforms
            % - Declaration: transform = IMDiscreteTransform(options)
            % - Parameter options.z: physical sample points
            % - Parameter options.weights: quadrature weights
            % - Parameter options.modeNumber: retained mode labels
            % - Parameter options.normalization: basis normalization name
            % - Parameter options.inverseMatrix: inverse transform matrix containing the sampled modes
            % - Parameter options.metricMatrix: sampled metric matrix
            % - Parameter options.targetGramMatrix: continuous diagonal Gram target
            % - Returns transform: initialized scalar discrete transform
            arguments
                options.z (:,1) double {mustBeReal, mustBeFinite}
                options.weights (:,1) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.normalization {mustBeTextScalar}
                options.inverseMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.metricMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.targetGramMatrix (:,:) double {mustBeReal, mustBeFinite}
            end

            z = options.z(:);
            weights = options.weights(:);
            inverseMatrix = options.inverseMatrix;
            metricMatrix = options.metricMatrix;
            targetGramMatrix = options.targetGramMatrix;
            nSamples = length(z);
            nModes = size(inverseMatrix,2);

            if length(weights) ~= nSamples || size(inverseMatrix,1) ~= nSamples
                error("IMDiscreteTransform:InvalidShape", "z, weights, and inverseMatrix rows must describe the same sample count.");
            end
            if length(options.modeNumber) ~= nModes
                error("IMDiscreteTransform:InvalidShape", "modeNumber must contain one label for each inverseMatrix column.");
            end
            if ~isequal(size(metricMatrix), [nSamples nSamples])
                error("IMDiscreteTransform:InvalidShape", "metricMatrix must have one row and column for each sample point.");
            end
            if ~isequal(size(targetGramMatrix), [nModes nModes])
                error("IMDiscreteTransform:InvalidShape", "targetGramMatrix must have one row and column for each retained mode.");
            end

            metricTolerance = 100*eps(max(1,norm(metricMatrix,2)));
            if norm(metricMatrix - metricMatrix.',2) > metricTolerance
                error("IMDiscreteTransform:NonSymmetricMetric", "metricMatrix must be symmetric.");
            end
            targetTolerance = 100*eps(max(1,norm(targetGramMatrix,2)));
            if norm(targetGramMatrix - diag(diag(targetGramMatrix)),2) > targetTolerance
                error("IMDiscreteTransform:InvalidTargetGramMatrix", "targetGramMatrix must be diagonal.");
            end
            targetNorms = diag(targetGramMatrix);
            normTolerance = 100*eps(max(1,max(abs(targetNorms))));
            if any(~isfinite(targetNorms)) || any(abs(targetNorms) <= normTolerance)
                error("IMDiscreteTransform:InvalidTargetGramMatrix", "targetGramMatrix must have finite, numerically nonzero diagonal entries.");
            end

            metricMatrix = 0.5*(metricMatrix + metricMatrix.');
            gramMatrix = inverseMatrix.'*metricMatrix*inverseMatrix;
            gramMatrix = 0.5*(gramMatrix + gramMatrix.');
            singularValues = svd(gramMatrix);
            rankTolerance = max(size(gramMatrix))*eps(max(1,norm(gramMatrix,2)));
            sampledGramRank = sum(singularValues > rankTolerance);
            pairingMatrix = inverseMatrix.'*metricMatrix;
            if sampledGramRank < nModes
                forwardMatrix = pinv(gramMatrix,rankTolerance)*pairingMatrix;
            else
                forwardMatrix = gramMatrix \ pairingMatrix;
            end

            scale = 1./sqrt(abs(targetNorms));
            scaledDifference = scale.*(gramMatrix - targetGramMatrix).*scale.';
            relativeGramOperatorError = norm(scaledDifference,2);
            if sampledGramRank < nModes
                relativeGramOperatorError = Inf;
            end

            self.z = z;
            self.weights = weights;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.normalization = string(options.normalization);
            self.inverseMatrix = inverseMatrix;
            self.metricMatrix = metricMatrix;
            self.forwardMatrix = forwardMatrix;
            self.gramMatrix = gramMatrix;
            self.targetGramMatrix = targetGramMatrix;
            self.relativeGramOperatorError = relativeGramOperatorError;
            self.roundTripError = norm(forwardMatrix*inverseMatrix - eye(nModes),2);
            self.inverseMatrixConditionNumber = cond(inverseMatrix);
            self.sampledGramRank = sampledGramRank;
            self.gramConditionNumber = cond(gramMatrix);
            self.targetGramIsPositiveDefinite = all(targetNorms > 0);
            self.hasNegativeWeights = any(weights < 0);
        end

        function coefficients = transformForward(self, values)
            % Transform sampled profiles forward to modal coefficients.
            %
            % For $$n_p$$ profiles arranged as
            % $$X\in\mathbb{R}^{n_z\times n_p}$$, this method returns
            % $$A\in\mathbb{R}^{n_m\times n_p}$$ using
            %
            % $$
            % A=A_{\mathrm f}X
            % =\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
            % A_{\mathrm i}^\mathsf{T}WX.
            % $$
            %
            % Each input column is transformed independently. The resulting
            % coefficients minimize the sampled quadratic residual when
            % $$W$$ is positive definite and, more generally, satisfy
            %
            % $$
            % A_{\mathrm i}^\mathsf{T}W
            % \left(X-A_{\mathrm i}A\right)=0.
            % $$
            %
            % ```matlab
            % coefficients = transform.transformForward(values);
            % coefficientsByMatrix = transform.forwardMatrix*values;
            % ```
            %
            % - Topic: Apply discrete transforms
            % - nav_order: 10
            % - Declaration: coefficients = transformForward(transform,values)
            % - Parameter values: `nSamples`-by-`nProfiles` sampled profile array with rows aligned to `z`
            % - Returns coefficients: `nModes`-by-`nProfiles` retained modal coefficient array
            arguments
                self IMDiscreteTransform
                values double
            end

            if size(values,1) ~= length(self.z)
                error("IMDiscreteTransform:InvalidSampleCount", "values must have one row for each sample point in z.");
            end
            coefficients = self.forwardMatrix*values;
        end

        function values = transformBack(self, coefficients)
            % Transform modal coefficients back to sampled profiles.
            %
            % For $$n_p$$ coefficient sets arranged as
            % $$A\in\mathbb{R}^{n_m\times n_p}$$, this method returns
            % $$\widehat{X}\in\mathbb{R}^{n_z\times n_p}$$ using
            %
            % $$
            % \widehat{X}=A_{\mathrm i}A=\Phi A.
            % $$
            %
            % Each coefficient column is transformed independently. For
            % coefficients obtained from `transformForward`, the result is
            % the sampled-space projection
            % $$\widehat{X}=A_{\mathrm i}A_{\mathrm f}X$$, not generally the
            % original profile unless it lies in the retained modal subspace.
            %
            % ```matlab
            % values = transform.transformBack(coefficients);
            % valuesByMatrix = transform.inverseMatrix*coefficients;
            % ```
            %
            % - Topic: Apply discrete transforms
            % - nav_order: 20
            % - Declaration: values = transformBack(transform,coefficients)
            % - Parameter coefficients: `nModes`-by-`nProfiles` coefficient array with rows aligned to `modeNumber`
            % - Returns values: `nSamples`-by-`nProfiles` reconstructed profile array sampled on `z`
            arguments
                self IMDiscreteTransform
                coefficients double
            end

            if size(coefficients,1) ~= length(self.modeNumber)
                error("IMDiscreteTransform:InvalidCoefficientCount", "coefficients must have one row for each retained mode.");
            end
            values = self.inverseMatrix*coefficients;
        end
    end
end
