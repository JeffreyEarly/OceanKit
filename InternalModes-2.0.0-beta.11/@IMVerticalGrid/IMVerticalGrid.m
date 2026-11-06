classdef IMVerticalGrid
    % Fixed physical sample points and quadrature weights for one column.
    %
    % A vertical grid is the sampling rule shared by every mode family and
    % horizontal wavenumber of a model. It is independent of the numerical
    % coordinate used to solve any eigenvalue problem: modes may be solved
    % at much higher resolution and then evaluated on these points.
    %
    % The default rule is the smoothed-bandwidth grid. It places mapped
    % Clenshaw-Curtis points in the coordinate $$x=\int\sqrt{Q}\,dz$$, where
    % $$Q-L_s^2Q_{zz}=B$$ with homogeneous Neumann conditions, and
    % $$B=\sqrt{(N^2+\tilde f^2-f_0^2)^2+4f_0^2\tilde f^2}$$ is the local
    % squared-frequency bandwidth. The default smoothing length is
    % $$L_s=2.5\,D/(N_z-1)$$.
    %
    % ```matlab
    % background = IMBackground(N2=@(z) 1e-5+5e-5*exp(z/300),zDomain=[-2000 0],surfaceCondition="freeSurface");
    % grid = IMVerticalGrid.smoothedBandwidth(background,65);
    % grid.z; grid.weights;
    % ```
    %
    % - Topic: Create a vertical grid
    % - Topic: Inspect a vertical grid
    % - Declaration: classdef IMVerticalGrid

    properties (SetAccess = private)
        % Strictly increasing physical sample points, including both endpoints.
        %
        % `z(1)` and `z(end)` equal `zDomain` exactly. The constructor
        % accepts boundary samples within 1e-10*D of the domain endpoints
        % and stores them snapped to those endpoints.
        %
        % - Topic: Inspect a vertical grid
        z

        % Nonnegative physical quadrature weights, summing to the column depth.
        %
        % Smoothed-bandwidth and physical-depth weights are strictly positive.
        % Mode-root weights can vanish where every fitting mode vanishes.
        %
        % - Topic: Inspect a vertical grid
        weights

        % Increasing physical domain `[bottom surface]` in meters.
        %
        % - Topic: Inspect a vertical grid
        zDomain

        % Grid rule: `"smoothedBandwidth"`, `"wkb"`, `"physicalDepth"`, `"modeRoots"` or `"custom"`.
        %
        % - Topic: Inspect a vertical grid
        kind

        % Construction parameters and diagnostics for this rule.
        %
        % - Topic: Inspect a vertical grid
        design

        % Explicit reconstruction and calculus for this grid, when supplied.
        % A geometry-only grid has empty calculus regardless of its kind.
        % - Topic: Inspect a vertical grid
        calculus
    end

    properties (Dependent, SetAccess = private)
        % Number of sample points.
        %
        % - Topic: Inspect a vertical grid
        nPoints

        % Positive column depth.
        %
        % - Topic: Inspect a vertical grid
        D
    end

    methods
        function self = IMVerticalGrid(options)
            % Create a grid from explicit points and weights.
            %
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid(options)
            % - Parameter options.z: strictly increasing points whose first and last values are the domain endpoints, within 1e-10*D; they are stored exactly as `zDomain`
            % - Parameter options.weights: nonnegative weights summing to the depth
            % - Parameter options.zDomain: increasing `[bottom surface]`; defaults to the point range
            % - Parameter options.kind: rule name, default `"custom"`
            % - Parameter options.design: provenance struct
            % - Parameter options.calculus: optional IMLobattoCalculus on exactly these nodes
            % - Returns self: immutable vertical grid
            arguments
                options.z (:,1) double {mustBeReal, mustBeFinite}
                options.weights (:,1) double {mustBeReal, mustBeFinite}
                options.zDomain double {mustBeReal} = []
                options.kind {mustBeTextScalar} = "custom"
                options.design (1,1) struct = struct()
                options.calculus = []
            end
            z = options.z;
            if isempty(options.zDomain), zDomain = [z(1) z(end)]; else, zDomain = options.zDomain; end
            if ~isequal(size(zDomain),[1 2]) || any(~isfinite(zDomain))
                error("IMVerticalGrid:InvalidDomain","zDomain must be a finite [bottom surface] pair.");
            end
            if numel(z) < 2 || any(diff(z) <= 0)
                error("IMVerticalGrid:InvalidPoints","Grid points must be strictly increasing, with at least two points.");
            end
            if zDomain(2) <= zDomain(1)
                error("IMVerticalGrid:InvalidDomain","zDomain must be strictly increasing.");
            end
            endpointTolerance = 1e-10*diff(zDomain);
            if abs(z(1)-zDomain(1)) > endpointTolerance || abs(z(end)-zDomain(2)) > endpointTolerance
                error("IMVerticalGrid:MissingEndpoint","The first and last points must be the domain endpoints zDomain(1) and zDomain(2).");
            end
            % Store the boundary samples exactly: endpoint-trace consumers
            % match zDomain to machine precision, not to this tolerance.
            z([1 end]) = zDomain;
            if any(diff(z) <= 0)
                error("IMVerticalGrid:InvalidPoints","Grid points must remain strictly increasing after snapping the endpoints to zDomain.");
            end
            if numel(options.weights) ~= numel(z) || any(options.weights < 0)
                error("IMVerticalGrid:InvalidWeights","Supply one nonnegative weight per point.");
            end
            depth = diff(zDomain);
            if abs(sum(options.weights)-depth) > 1e-10*depth
                error("IMVerticalGrid:InvalidWeights","Weights must integrate a constant exactly: sum(weights) must equal the depth.");
            end
            self.z = z;
            self.weights = options.weights;
            self.zDomain = zDomain;
            self.kind = string(options.kind);
            self.design = options.design;
            if ~isempty(options.calculus)
                if ~isa(options.calculus,"IMLobattoCalculus") || ~isscalar(options.calculus) || ...
                        ~isequal(options.calculus.depthAtCoordinate(options.calculus.nativeCoordinate),z)
                    error("IMVerticalGrid:IncompatibleCalculus","The supplied calculus must reconstruct these exact depths in order.");
                end
            end
            self.calculus = options.calculus;
        end

        function value = get.nPoints(self)
            value = numel(self.z);
        end

        function value = get.D(self)
            value = diff(self.zDomain);
        end

        function value = integrate(self,values)
            % Apply the quadrature rule along the first dimension.
            %
            % - Topic: Inspect a vertical grid
            % - Declaration: value = integrate(grid,values)
            % - Parameter values: samples with one row per grid point
            % - Returns value: weighted sums over the first dimension
            if size(values,1) ~= self.nPoints
                error("IMVerticalGrid:SampleMismatch","values must have one row per grid point.");
            end
            value = reshape(self.weights.'*reshape(values,self.nPoints,[]),[1 size(values,2:ndims(values))]);
        end

        function state = snapshot(self)
            % Export geometry, provenance and any explicit calculus.
            % - Topic: Inspect a vertical grid
            % - Declaration: state = grid.snapshot()
            % - Returns state: versioned scalar struct
            calculusState = [];
            if ~isempty(self.calculus), calculusState = self.calculus.snapshot(); end
            state = struct(version=1,z=self.z,zDomain=self.zDomain,weights=self.weights, ...
                kind=self.kind,design=self.design,calculus=calculusState);
        end
    end

    methods (Static)
        function grid = fromSnapshot(state)
            % Restore saved geometry and calculus without rebuilding the map from a profile.
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid.fromSnapshot(state)
            % - Parameter state: scalar struct from snapshot
            % - Returns grid: restored vertical grid
            arguments
                state (1,1) struct
            end
            fields = ["version","z","zDomain","weights","kind","design","calculus"];
            if ~all(isfield(state,fields)) || ~isequal(state.version,1)
                error("IMVerticalGrid:InvalidSnapshot","A supported vertical grid snapshot is required.");
            end
            calculus = [];
            if ~isempty(state.calculus), calculus = IMLobattoCalculus.fromSnapshot(state.calculus); end
            grid = IMVerticalGrid(z=state.z(:),weights=state.weights(:),zDomain=state.zDomain(:).', ...
                kind=state.kind,design=state.design,calculus=calculus);
        end

        function Ls = defaultSmoothingLength(D,nPoints)
            % Return the provisional default smoothing length 2.5*D/(nPoints-1).
            %
            % The rule scales with the mean grid spacing. Joint-recovery
            % comparisons (study joint correction, retained-mode data, one
            % 1000 m column with a neutral surface layer, 10 km wavelength)
            % found a recovery plateau of 2-4.7 spacings at 65 points and the
            % smallest tested boundary errors near 2.4-3.4 spacings at 33
            % points, with mandatory boundary modes failing at 65 points
            % below about 1.2 spacings. This is provisional.
            %
            % - Topic: Create a vertical grid
            % - Declaration: Ls = IMVerticalGrid.defaultSmoothingLength(D,nPoints)
            arguments
                D (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                nPoints (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPoints,3)}
            end
            Ls = 2.5*D/(nPoints-1);
        end

        function grid = physicalDepth(zDomain,nPoints)
            % Create Clenshaw-Curtis points and weights in physical depth.
            %
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid.physicalDepth(zDomain,nPoints)
            % - Parameter zDomain: increasing `[bottom surface]`, or an IMBackground
            % - Parameter nPoints: number of points, at least three
            % - Returns grid: physical-depth Chebyshev-Lobatto grid
            arguments
                zDomain
                nPoints (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPoints,3)}
            end
            if isa(zDomain,"IMBackground"), zDomain = zDomain.zDomain; end
            [z,w] = IMVerticalGrid.clenshawCurtis(nPoints,zDomain);
            metric = chebfun(1,zDomain);
            calculus = IMLobattoCalculus.fromMetric(metric,z,z-z(1));
            grid = IMVerticalGrid(z=z,weights=w,zDomain=zDomain,kind="physicalDepth", ...
                design=struct(rule="Clenshaw-Curtis in physical depth"),calculus=calculus);
        end

        function grid = wkb(background,nPoints)
            % Create a Lobatto grid in the classical WKB coordinate.
            %
            % The coordinate metric is $$dx/dz=\sqrt{N^2(z)}$$ and must be
            % positive throughout the reference column.
            %
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid.wkb(background,nPoints)
            % - Parameter background: column background with positive stratification
            % - Parameter nPoints: number of points, at least three
            % - Returns grid: WKB grid with mapped quadrature and calculus
            arguments
                background (1,1) IMBackground
                nPoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nPoints,3)}
            end
            metric = chebfun(@(z) sqrt(background.N2(z)),background.zDomain);
            [calculus,z,w] = IMLobattoCalculus.forMetric(metric,nPoints);
            design = struct(rule="wkb",weightRule="Mapped Clenshaw-Curtis, normalized to integrate the column depth");
            grid = IMVerticalGrid(z=z,weights=w,zDomain=background.zDomain,kind="wkb",design=design,calculus=calculus);
        end

        function grid = smoothedBandwidth(background,nPoints,options)
            % Create the default smoothed-bandwidth grid for a background.
            %
            % The smoothing changes only where points are placed; it does not
            % alter the physical stratification. `tildeF` defaults to the
            % full horizontal Coriolis component at the background latitude.
            % Supply `breakpoints` at known kinks of N2 (for example the ends
            % of a neutral layer) to represent the profile accurately.
            %
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid.smoothedBandwidth(background,nPoints,options)
            % - Parameter background: IMBackground
            % - Parameter nPoints: number of points, at least three
            % - Parameter options.smoothingLength: positive length in meters; default 2.5*D/(nPoints-1)
            % - Parameter options.breakpoints: interior depths where N2 is not smooth
            % - Returns grid: smoothed-bandwidth grid with filter and map diagnostics in `design`
            arguments
                background (1,1) IMBackground
                nPoints (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPoints,3)}
                options.smoothingLength double {mustBeReal} = []
                options.breakpoints (1,:) double {mustBeReal, mustBeFinite} = zeros(1,0)
            end
            domain = background.zDomain;
            D = background.D;
            if isempty(options.smoothingLength)
                Ls = IMVerticalGrid.defaultSmoothingLength(D,nPoints);
            else
                Ls = options.smoothingLength;
                if ~isscalar(Ls) || ~isfinite(Ls) || Ls <= 0
                    error("IMVerticalGrid:InvalidSmoothingLength","smoothingLength must be a positive finite scalar.");
                end
            end
            interior = unique(options.breakpoints(options.breakpoints > domain(1) & options.breakpoints < domain(2)));
            breaks = [domain(1) interior domain(2)];
            f = background.f0;
            tildeF = 2*background.rotationRate*cosd(background.latitude);
            N2 = background.N2;
            BFunction = @(z) hypot(N2(z)+tildeF^2-f^2,2*f*tildeF);
            zCheck = unique([linspace(domain(1),domain(2),2001),breaks]);
            BCheck = BFunction(zCheck);
            if any(~isfinite(BCheck)) || ~isreal(BCheck) || max(BCheck) <= 0
                error("IMVerticalGrid:InvalidBandwidth","The bandwidth must be finite, real and positive somewhere; supply N2 or nonzero rotation.");
            end
            B = chebfun(BFunction,breaks,'splitting','on');
            meanB = sum(B)/D;
            % Residuals are evaluated one-sided at interior breakpoints: at a
            % jump in N2, Q_zz jumps and a two-sided evaluation would mix the
            % averaged Chebfun value with one side of the source.
            zResidual = IMVerticalGrid.oneSidedCheckPoints(domain,interior,2001);
            BResidual = BFunction(zResidual);
            [valueQ,derivativeQ,secondQ,conditioning] = IMVerticalGrid.neumannFilter(B,breaks,Ls,D,"direct");
            filterResidual = max(abs(valueQ(zResidual)-Ls^2*secondQ(zResidual)-BResidual))/max(BCheck);
            neumannResidual = D*max(abs(derivativeQ(domain)))/max(abs(valueQ(zCheck)));
            if filterResidual > 1e-9 || neumannResidual > 1e-9
                % A nearly constant Q can lose small curvature when simplified;
                % solve for its variation about the preserved mean instead.
                [valueQ,derivativeQ,secondQ,conditioning] = IMVerticalGrid.neumannFilter(B-meanB,breaks,Ls,D,"centered",meanB);
                filterResidual = max(abs(valueQ(zResidual)-Ls^2*secondQ(zResidual)-BResidual))/max(BCheck);
                neumannResidual = D*max(abs(derivativeQ(domain)))/max(abs(valueQ(zCheck)));
            end
            QCheck = valueQ(zCheck);
            if min(QCheck) <= 0 || filterResidual > 1e-9 || neumannResidual > 1e-9
                error("IMVerticalGrid:UnqualifiedFilter","The smoothing filter failed its positivity or 1e-9 residual check.");
            end
            qFunction = @(z) sqrt(valueQ(z));
            q = chebfun(qFunction,breaks,'splitting','on');
            qMin = min(q); qMax = max(q);
            map = cumsum(q);
            map = map-feval(map,domain(1));
            lengthX = feval(map,domain(2));
            if qMax-qMin <= 50*eps(qMax)
                [z,rawWeights] = IMVerticalGrid.clenshawCurtis(nPoints,domain);
                x = (z-domain(1))*lengthX/D;
            else
                [x,wx] = IMVerticalGrid.clenshawCurtis(nPoints,[0 lengthX]);
                z = zeros(nPoints,1);
                z([1 end]) = domain(:);
                for j = 2:nPoints-1
                    z(j) = fzero(@(depth) feval(map,depth)-x(j),domain);
                end
                rawWeights = wx./feval(q,z);
            end
            rawSum = sum(rawWeights);
            w = (D/rawSum)*rawWeights;
            inversionResidual = max(abs(feval(map,z)-x))/lengthX;
            if any(diff(z) <= 0) || any(~isfinite(w) | w <= 0) || inversionResidual > 1e-10
                error("IMVerticalGrid:UnqualifiedMap","The coordinate map must be invertible with positive weights.");
            end
            design = struct(rule="smoothedBandwidth",smoothingLength=Ls, ...
                smoothingLengthInGridSpacings=Ls*(nPoints-1)/D,tildeF=tildeF,breakpoints=breaks, ...
                filterEquation="Q-Ls^2*Q_zz=B, Q_z=0 at both endpoints",filterConditioning=conditioning, ...
                filterRelativeResidual=filterResidual,filterNeumannResidual=neumannResidual, ...
                qMin=qMin,qMax=qMax,nativeCoordinate=x,coordinateLength=lengthX,mapResidual=inversionResidual, ...
                weightRule="Mapped Clenshaw-Curtis with physical Jacobian, normalized to integrate a constant exactly", ...
                unadjustedRelativeDepthError=abs(rawSum-D)/D);
            calculus = IMLobattoCalculus.fromMetric(q,z,x);
            grid = IMVerticalGrid(z=z,weights=w,zDomain=domain,kind="smoothedBandwidth",design=design,calculus=calculus);
        end

        function grid = modeRoots(basisSet,nPoints,options)
            % Create the classical mode-root grid with fitted weights.
            %
            % Points are both endpoints plus the roots of a generating mode
            % from `basisSet`, as in Early, Lelong and Smith (2020). Weights
            % are the nonnegative, full-depth fit of `quadratureWeightsForPoints`
            % to the leading `nFitModes` modes.
            %
            % - Topic: Create a vertical grid
            % - Declaration: grid = IMVerticalGrid.modeRoots(basisSet,nPoints,options)
            % - Parameter basisSet: numerical IMBasisSet supplying the generating modes
            % - Parameter nPoints: exact number of points
            % - Parameter options.nFitModes: modes in the weight fit; default nPoints-2, capped by the basis
            % - Returns grid: mode-root grid; `design` holds root provenance and the weight fit
            arguments
                basisSet (1,1) IMBasisSet
                nPoints (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPoints,3)}
                options.nFitModes double = []
            end
            [z,gridDesign] = basisSet.modeRootGrid(nPoints=nPoints);
            nFit = options.nFitModes;
            if isempty(nFit), nFit = min(nPoints-2,size(basisSet.nativeModes,2)); end
            [w,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nFit);
            zDomain = basisSet.zDomain;
            w = w*(diff(zDomain)/sum(w));
            design = struct(rule="modeRoots",gridDesign=gridDesign,nFitModes=nFit,weightFit=weightFit);
            grid = IMVerticalGrid(z=z(:),weights=w(:),zDomain=zDomain,kind="modeRoots",design=design);
        end
    end

    methods (Static, Access = private)
        function [z,w] = clenshawCurtis(nPoints,domain)
            n = nPoints-1;
            theta = (0:n).'*pi/n;
            weights = zeros(nPoints,1);
            interior = 2:n;
            v = ones(n-1,1);
            if mod(n,2) == 0
                weights([1 end]) = 1/(n^2-1);
                for k = 1:n/2-1
                    v = v-2*cos(2*k*theta(interior))/(4*k^2-1);
                end
                v = v-cos(n*theta(interior))/(n^2-1);
            else
                weights([1 end]) = 1/n^2;
                for k = 1:(n-1)/2
                    v = v-2*cos(2*k*theta(interior))/(4*k^2-1);
                end
            end
            weights(interior) = 2*v/n;
            halfLength = diff(domain)/2;
            z = flipud(halfLength*cos(theta)+mean(domain));
            z([1 end]) = domain(:);
            w = halfLength*flipud(weights);
        end

        function z = oneSidedCheckPoints(domain,interiorBreaks,nPoints)
            % Uniform check points with each interior breakpoint replaced by
            % points just inside its two adjacent pieces.
            offset = 1e-8*diff(domain);
            z = linspace(domain(1),domain(2),nPoints);
            for breakpoint = interiorBreaks
                z = z(abs(z-breakpoint) > offset);
                z = [z breakpoint-offset breakpoint+offset]; %#ok<AGROW>
            end
            z = sort(z);
        end

        function [valueQ,derivativeQ,secondQ,conditioning] = neumannFilter(source,breaks,Ls,D,form,meanB)
            domain = breaks([1 end]);
            if form == "direct"
                filter = chebop(domain);
                filter.op = @(z,Q) Q-Ls^2*diff(Q,2);
                filter.lbc = @(Q) diff(Q);
                filter.rbc = @(Q) diff(Q);
                Q = filter\source;
                Qz = diff(Q); Qzz = diff(Q,2);
                valueQ = @(z) feval(Q,z);
                derivativeQ = @(z) feval(Qz,z);
                secondQ = @(z) feval(Qzz,z);
                conditioning = "direct";
            else
                centered = chebop(breaks);
                centered.op = @(z,U) U/Ls^2-diff(U,2);
                centered.lbc = @(U) diff(U);
                centered.rbc = @(U) diff(U);
                preferences = cheboppref();
                preferences.bvpTol = eps;
                preferences.discretization = 'coeffs';
                U = solvebvp(centered,source,preferences);
                Uz = diff(U); Uzz = diff(U,2);
                meanU = sum(U)/D;
                valueQ = @(z) meanB+(feval(U,z)-meanU)/Ls^2;
                derivativeQ = @(z) feval(Uz,z)/Ls^2;
                secondQ = @(z) feval(Uzz,z)/Ls^2;
                conditioning = "centeredConstantPlusVariation";
            end
        end
    end
end
