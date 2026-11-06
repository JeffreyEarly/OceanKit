classdef IMGeostrophicModes < IMInternalModes
    % Describe regular paired geostrophic eigenvalue problems.
    %
    % Numerical solves retain both fields in
    % $$F_\xi=-N^2G/g,\quad G_\xi=rF,\quad r=1/h.$$
    % Nonnegative stratification, including neutral intervals, is supported.
    % Spectral solves share one coordinate for both fields. Finite differences
    % use F on the supplied grid and G on a conservative staggered grid.
    % The inherited scalar coefficients describe the equivalent problem
    % only where N2 is positive; numerical assembly never evaluates them.
    % Independent zero-energy neutral-label states are not represented.
    %
    % Construct through `IMInternalModes.geostrophicAPVModes` or
    % `IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes`.
    %
    % ```matlab
    % evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=0,gd=-0.01,surfaceBoundary="freeSurface");
    % basis = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=20);
    % ```
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Inspect internal-mode inner products
    % - Topic: Inspect mode selection
    % - Topic: Developer topics
    % - Declaration: classdef IMGeostrophicModes < IMInternalModes

    methods
        function self = IMGeostrophicModes(options)
            % Create a paired descriptor from a geostrophic factory configuration.
            %
            % - Topic: Developer topics
            % - Declaration: evp = IMGeostrophicModes(options)
            % - Developer: true
            arguments
                options.name
                options.zDomain (1,2) double
                options.N2 function_handle
                options.formulation = "F"
                options.modeFamily = "hydrostatic"
                options.p
                options.q
                options.r
                options.surfaceBoundary IMBoundaryCondition
                options.bottomBoundary IMBoundaryCondition
                options.rotationRate (1,1) double = 7.2921e-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                options.hFromEigenvalue function_handle = @(lambda) 1./lambda
                options.parameters struct
            end
            args = namedargs2cell(options);
            self@IMInternalModes(args{:});
            self.stratificationValues(linspace(self.zDomain(1),self.zDomain(2),2001).');
        end

        function basis = makeBasisSet(self,solver,V,lambda,modeNumber,diagnostics)
            % Construct aligned F/G fields with stable equivalent depths.
            %
            % - Topic: Developer topics
            % - Declaration: basis = makeBasisSet(evp,solver,V,lambda,modeNumber,diagnostics)
            % - Developer: true
            if isa(solver,"IMSolverFiniteDifference")
                basis = IMGeostrophicFiniteDifferenceTools.makeBasis(self,solver,V,lambda,modeNumber,diagnostics);
                return
            end
            n = numel(solver.zNative);
            if size(V,1) ~= 2*n
                error("IMGeostrophicModes:InvalidNativeModes","Paired native modes must contain n F coefficients followed by n G coefficients.");
            end
            if isfield(diagnostics,"geostrophicInverseDepth")
                inverseDepth = diagnostics.geostrophicInverseDepth;
                diagnostics = rmfield(diagnostics,"geostrophicInverseDepth");
            elseif self.name == "geostrophicAPVModes"
                inverseDepth = lambda;
            else
                inverseDepth = self.g*(lambda-self.parameters.k^2)/self.f0^2;
            end
            basis = IMGeostrophicModesBasis(solver=solver,evp=self,nativeModes=V(1:n,:), ...
                nativeGModes=V(n+1:end,:),eigenvalues=lambda,inverseDepth=reshape(inverseDepth,1,[]), ...
                modeNumber=modeNumber,modeSelectionDiagnostics=diagnostics,metadata=self.parameters);
        end

        function spec = innerProduct(self,variable)
            % Return the manuscript pairing, including companion endpoint fields.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            arguments
                self IMGeostrophicModes
                variable {mustBeTextScalar,mustBeMember(variable,["F","G"])} = "F"
            end
            spec = innerProduct@IMInternalModes(self,variable);
            if self.name == "geostrophicAPVModes" && string(variable) == "G"
                % Physical endpoint coefficients come from the declared
                % accelerations. Relative tests on dimensional canonical
                % boundary coefficients can erase finite value observations.
                spec.surfaceWeights = struct([]);
                spec.bottomWeights = struct([]);
                spec.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();
                g0 = self.parameters.g0;
                gd = self.parameters.gd;
                if self.parameters.surfaceBoundary == "freeSurface"
                    if isinf(g0)
                        surfaceCoefficient = 1;
                    elseif g0 == -self.g
                        surfaceCoefficient = 0; % Gs is constrained to zero.
                    else
                        surfaceCoefficient = g0/(self.g+g0);
                        if isinf(self.g+g0)
                            surfaceCoefficient = 1/(1+self.g/g0);
                        end
                    end
                elseif isfinite(g0)
                    surfaceCoefficient = g0/self.g;
                else
                    surfaceCoefficient = 0;
                end
                if isfinite(gd), bottomCoefficient = gd/self.g; else, bottomCoefficient = 0; end
                surfaceActive = isfinite(g0) && g0 ~= 0 && (self.parameters.surfaceBoundary ~= "freeSurface" || g0 ~= -self.g);
                surfaceActive = surfaceActive || (self.parameters.surfaceBoundary == "freeSurface" && isinf(g0));
                coefficients = [surfaceCoefficient,bottomCoefficient];
                if any(~isfinite(coefficients)) || (surfaceActive && surfaceCoefficient == 0) ...
                        || (isfinite(gd) && gd ~= 0 && bottomCoefficient == 0)
                    error("IMGeostrophicModes:UnrepresentablePairing","A nonzero physical APV endpoint coefficient exceeds floating-point range. Choose representable physical units.");
                end
                locations = ["surface","bottom"];
                for i = find(coefficients ~= 0)
                    term = struct(location=locations(i),coefficient=coefficients(i),variable="G", ...
                        description="G("+locations(i)+")",catalogCase="declaredAPVAcceleration");
                    spec.endpointInnerProductTerms(end+1,1) = term;
                end
                spec.hasInnerProduct = true;
                spec.reason = "The declared APV accelerations define the signed G value pairing, including literal inactive and constrained-zero limits.";
            elseif self.name == "geostrophicGeneralizedPotentialEnstrophyModes" && string(variable) == "G"
                spec.hasInnerProduct = true;
                spec.reason = "The paired generalized-enstrophy problem supplies companion F endpoint traces.";
                k = self.parameters.k;
                spec.endpointInnerProductTerms = struct("location",{"surface","bottom"}, ...
                    "variable",{"F","F"},"coefficient",{1-self.g*k^2/self.parameters.alpha0,-self.g*k^2/self.parameters.alphaD});
            end
        end

        function diagnostics = modeSelectionDiagnostics(self,solver,~)
            % Classify observable signed modes from the regular endpoint problem.
            %
            % - Topic: Inspect mode selection
            % - Declaration: diagnostics = modeSelectionDiagnostics(evp,solver,A)
            solver = solver.configuredForEVP(self);
            [nNegative,zeroPresent] = self.endpointClassification(solver);
            diagnostics = struct("assessmentLevel","endpointInertia", ...
                "negativeEndpointWeightCount",nNegative,"minNegativeEigenvalueCount",nNegative, ...
                "maxNegativeEigenvalueCount",nNegative,"zeroModeStatus","absent", ...
                "zeroModeCount",0,"zeroModeSingularValue",NaN,"zeroModeTolerance",NaN, ...
                "reason","Regular geostrophic endpoint classification; no inverse N2 samples.");
            diagnostics.zeroInverseDepthCount = double(zeroPresent);
            if self.name == "geostrophicAPVModes" && zeroPresent
                diagnostics.zeroModeStatus = "present";
                diagnostics.zeroModeCount = 1;
            end
        end

        function bounds = negativeEigenvalueBounds(self,solver,A)
            % Bound signed physical modes using endpoint inertia.
            %
            % - Topic: Inspect mode selection
            % - Declaration: bounds = negativeEigenvalueBounds(evp,solver,A)
            if nargin < 3, A = []; end
            bounds = self.modeSelectionDiagnostics(solver,A);
        end

        function diagnostics = definitenessDiagnostics(self,solver)
            % Assess the regular positive F metric and endpoint energy.
            %
            % - Topic: Inspect mode selection
            % - Declaration: diagnostics = definitenessDiagnostics(evp,solver)
            solver = solver.configuredForEVP(self);
            N2 = self.stratificationValues(solver.zNative);
            [negativeCount,~] = self.endpointClassification(solver);
            diagnostics = struct("assessmentLevel","regularPair","pMin",NaN,"qMin",0,"rMin",1, ...
                "pPositive",all(N2>0),"qNonnegative",true,"rPositive",true, ...
                "endpointWeights",self.endpointWeights(),"negativeEndpointWeightCount",0, ...
                "metricPositive",true,"hasDegenerateEndpointMetric",false, ...
                "endpointNumeratorNegativeDirections",negativeCount,"endpointNumeratorNonnegative",negativeCount==0, ...
                "interiorNonnegative",true,"quadraticFormNonnegative",negativeCount==0, ...
                "reason","The regular paired problem has a positive F metric, including neutral intervals.");
        end

        function selection = selectModes(self,lambda,nModes,solver,A,options)
            % Select using physical units without an absolute eigenvalue floor.
            %
            % - Topic: Developer topics
            % - Declaration: selection = selectModes(evp,lambda,nModes,solver,A,options)
            % - Developer: true
            arguments
                self IMGeostrophicModes
                lambda (:,1) double {mustBeReal,mustBeFinite}
                nModes (1,1) double {mustBeInteger,mustBePositive}
                solver IMSolver
                A (:,:) double
                options.diagnostics struct = struct()
            end
            diagnostics = options.diagnostics;
            if isempty(fieldnames(diagnostics)), diagnostics = self.modeSelectionDiagnostics(solver,A); end
            negative = find(lambda<0);
            if numel(negative)<min(nModes,diagnostics.minNegativeEigenvalueCount) ...
                    || (nModes>diagnostics.minNegativeEigenvalueCount && diagnostics.zeroModeStatus=="present" && ~any(lambda==0))
                error("IMGeostrophicModes:UnresolvedPhysicalSector", ...
                    "The requested band is missing a negative or zero mode required by the endpoint problem. Increase nEVP.");
            end
            [~,order] = sort(lambda(negative));
            negative = negative(order(1:min(numel(order),diagnostics.maxNegativeEigenvalueCount)));
            zero = find(lambda==0,1);
            if diagnostics.zeroModeStatus ~= "present", zero = []; end
            positive = find(lambda>0);
            [~,order] = sort(lambda(positive));
            positive = positive(order);
            index = [negative(:);zero(:);positive(:)];
            labels = [-1:-1:-numel(negative),zeros(1,numel(zero)),1:numel(positive)];
            count = min(nModes,numel(index));
            selection = struct("sortIndex",index(1:count),"modeNumber",labels(1:count),"modeSelectionDiagnostics",diagnostics);
        end

        function summarize(self,solver)
            % Summarize the regular pair without evaluating scalar inverse N2.
            %
            % - Topic: Inspect mode selection
            % - Declaration: summarize(evp,solver)
            if nargin < 2, solver = []; end
            fprintf('%s: F_xi = -N2 G/g, G_xi = F/h\n',self.name);
            fprintf('  domain: [%g %g] m; f0: %g s^-1; g: %g m s^-2\n',self.zDomain,self.f0,self.g);
            fprintf('  formulation: paired F/G; nonnegative N2; h from solved inverse depth\n');
            if ~isempty(solver), disp(self.modeSelectionDiagnostics(solver,[])); end
        end
    end

    methods (Hidden)
        [A,B,samples] = assembleConfigured(self,solver)
        [A,B,map] = prepareSolvePencil(self,solver,A,B)
        [V,r] = prepareGeneralizedEigenpairs(self,solver,A,B,V,r)
        function diagnostics = preparedModeSelectionDiagnostics(self,~,A) %#ok<INUSD>
            diagnostics = struct();
        end
    end

    methods (Access = private)
        [A,B] = inverseDepthPencil(self,solver)
        [negativeCount,zeroPresent,IN] = endpointClassification(self,solver)
        lambda = publicEigenvalues(self,solver,V,r)
        valid = qualifiedInverseDepthPairs(self,solver,V,r)
        [V,r] = refineSmallInverseDepthModes(self,solver,V,r,IN)
        [r,V] = fullyNeutralMode(self,solver)
        function values = stratificationValues(self,z)
            values = self.N2(z);
            if ~isreal(values) || numel(values)~=numel(z) || any(~isfinite(values(:))) || any(values(:)<0)
                error("IMGeostrophicModes:InvalidStratification","N2 must return finite real nonnegative values at every requested depth.");
            end
            values = values(:);
        end
    end
end
