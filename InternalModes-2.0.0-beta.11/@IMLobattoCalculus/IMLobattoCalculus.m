classdef IMLobattoCalculus
    % Reconstruct and differentiate mapped Chebyshev-Lobatto samples.
    %
    % One global polynomial in the increasing native coordinate reconstructs
    % each column of samples. A positive metric $$q=dx/dz$$ defines the map
    % $$x(z)=\int_{z_b}^z q(s)\,ds$$. The metric may have several pieces.
    %
    % ```matlab
    % grid = IMVerticalGrid.physicalDepth([-1000 0],33);
    % Dz = grid.calculus.firstDerivativeMatrix();
    % ```
    %
    % - Topic: Apply vertical calculus
    % - Topic: Restore vertical calculus
    % - Declaration: classdef IMLobattoCalculus

    properties (SetAccess = private)
        % Native coordinates of the Lobatto nodes in increasing order.
        % - Topic: Apply vertical calculus
        nativeCoordinate
    end

    properties (Access = private)
        z
        metric
        map
        barycentricWeights
    end

    methods
        function Dz = firstDerivativeMatrix(self)
            % Map nodal values to first derivatives in reference depth.
            %
            % The nodal matrix is $$D_z=\operatorname{diag}(dx/dz)D_x$$.
            % Applying this matrix repeatedly reconstructs each intermediate
            % derivative from its nodal samples.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: Dz = calculus.firstDerivativeMatrix()
            % - Returns Dz: square nodal derivative matrix in inverse meters
            Dz = self.derivativeAtNodes(self.z,self.nativeCoordinate);
        end

        function P = interpolationMatrix(self,zQuery)
            % Evaluate the native interpolant at depths in the domain.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: P = calculus.interpolationMatrix(zQuery)
            % - Parameter zQuery: query depths in meters
            % - Returns P: matrix with one row per depth and one column per sample
            x = self.coordinateAtDepth(zQuery);
            P = self.interpolationAtNativePoints(x(:));
        end

        function x = coordinateAtDepth(self,zQuery)
            % Evaluate the native coordinate at in-domain depths.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: x = calculus.coordinateAtDepth(zQuery)
            % - Parameter zQuery: query depths in meters
            % - Returns x: native coordinates with the shape of zQuery
            arguments
                self (1,1) IMLobattoCalculus
                zQuery double {mustBeReal,mustBeFinite}
            end
            if any(zQuery < self.z(1) | zQuery > self.z(end),"all")
                error("IMLobattoCalculus:OutsideDomain","Query depths must lie within the grid domain.");
            end
            x = feval(self.map,zQuery(:));
            [atNode,index] = ismember(zQuery(:),self.z);
            x(atNode) = self.nativeCoordinate(index(atNode));
            x = reshape(x,size(zQuery));
        end

        function zQuery = depthAtCoordinate(self,x)
            % Invert the native coordinate at in-domain points.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: zQuery = calculus.depthAtCoordinate(x)
            % - Parameter x: native coordinates
            % - Returns zQuery: reference depths in meters with the shape of x
            arguments
                self (1,1) IMLobattoCalculus
                x double {mustBeReal,mustBeFinite}
            end
            if any(x < self.nativeCoordinate(1) | x > self.nativeCoordinate(end),"all")
                error("IMLobattoCalculus:OutsideDomain","Query coordinates must lie within the native domain.");
            end
            [atNode,index] = ismember(x(:),self.nativeCoordinate);
            zQuery = zeros(numel(x),1);
            zQuery(~atNode) = IMLobattoCalculus.invertMap(self.map,x(~atNode));
            zQuery(atNode) = self.z(index(atNode));
            zQuery = reshape(zQuery,size(x));
        end

        function [zQuery,weights] = quadratureRule(self,nPoints)
            % Map a Gauss-Legendre rule into reference depth.
            %
            % The positive weights include the physical Jacobian. Their sum
            % approximates the column depth without a constant renormalization.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: [zQuery,weights] = calculus.quadratureRule(nPoints)
            % - Parameter nPoints: number of quadrature points
            % - Returns zQuery: increasing interior depths in meters
            % - Returns weights: corresponding quadrature weights in meters
            arguments
                self (1,1) IMLobattoCalculus
                nPoints (1,1) double {mustBeInteger,mustBePositive}
            end
            [x,w] = legpts(nPoints,self.nativeCoordinate([1 end]).');
            zQuery = self.depthAtCoordinate(x(:));
            weights = w(:)./feval(self.metric,zQuery);
        end

        function [zQuery,weights,Dz] = lobattoRule(self,nPoints)
            % Resample the same map at another Lobatto resolution.
            %
            % The map and its smoothing length are fixed. The returned
            % Clenshaw-Curtis weights integrate constants exactly.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: [zQuery,weights,Dz] = calculus.lobattoRule(nPoints)
            % - Parameter nPoints: number of Lobatto nodes, including endpoints
            % - Returns zQuery: increasing reference depths in meters
            % - Returns weights: mapped quadrature weights in meters
            % - Returns Dz: physical first-derivative matrix on zQuery
            arguments
                self (1,1) IMLobattoCalculus
                nPoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nPoints,3)}
            end
            [x,w] = chebpts(nPoints,self.nativeCoordinate([1 end]).');
            zQuery = self.depthAtCoordinate(x);
            weights = w(:)./feval(self.metric,zQuery);
            weights = diff(self.z([1 end]))*weights/sum(weights);
            if nargout >= 3
                Dz = self.derivativeAtNodes(zQuery,x);
            end
        end

        function values = integrateFromSurface(self,samples,zQuery)
            % Integrate the reconstructed samples from each depth to the surface.
            %
            % Gauss-Legendre quadrature in the native coordinate accounts for
            % $$dz/dx$$ and splits at the metric's internal breakpoints.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: values = calculus.integrateFromSurface(samples,zQuery)
            % - Parameter samples: columns sampled on the Lobatto grid
            % - Parameter zQuery: lower integration limits in meters
            % - Returns values: one integral per depth and sample column
            arguments
                self (1,1) IMLobattoCalculus
                samples (:,:) double
                zQuery double {mustBeReal,mustBeFinite}
            end
            self.validateSamples(samples);
            x = self.coordinateAtDepth(zQuery(:));
            xBreaks = self.coordinateAtDepth(self.metric.domain);
            [nodes,w] = legpts(max(32,2*numel(self.z)));
            values = zeros(numel(x),size(samples,2),"like",samples);
            for j = 1:numel(x)
                if x(j) == xBreaks(end), continue; end
                bounds = [x(j),xBreaks(xBreaks > x(j))];
                halfWidths = diff(bounds)/2;
                centers = (bounds(1:end-1)+bounds(2:end))/2;
                xq = nodes(:)*halfWidths+centers;
                zq = self.depthAtCoordinate(xq(:));
                weights = reshape(w(:)*halfWidths,[],1)./feval(self.metric,zq);
                values(j,:) = (weights.'*self.interpolationAtNativePoints(xq(:)))*samples;
            end
        end

        function values = endpointTraces(self,samples,options)
            % Return bottom and surface values or reference-depth derivatives.
            %
            % - Topic: Apply vertical calculus
            % - Declaration: values = calculus.endpointTraces(samples,options)
            % - Parameter samples: columns sampled on the Lobatto grid
            % - Parameter options.derivativeOrder: zero for values or one for first derivatives
            % - Returns values: two rows in bottom-to-surface order
            arguments
                self (1,1) IMLobattoCalculus
                samples (:,:) double
                options.derivativeOrder (1,1) double {mustBeMember(options.derivativeOrder,[0 1])} = 0
            end
            self.validateSamples(samples);
            if options.derivativeOrder == 1
                Dz = self.firstDerivativeMatrix();
                values = Dz([1 end],:)*samples;
            else
                values = samples([1 end],:);
            end
        end

        function state = snapshot(self)
            % Export a numeric representation of the metric and sample grid.
            %
            % - Topic: Restore vertical calculus
            % - Declaration: state = calculus.snapshot()
            % - Returns state: versioned scalar struct of numeric fields
            pieces = cellfun(@chebcoeffs,self.metric.funs,UniformOutput=false);
            state = struct(version=1,z=self.z,nativeCoordinate=self.nativeCoordinate, ...
                metricBreaks=self.metric.domain,metricCoefficientCounts=cellfun(@numel,pieces), ...
                metricCoefficients=vertcat(pieces{:}),metricPointValues=feval(self.metric,self.metric.domain));
        end
    end

    methods (Static)
        function self = fromSnapshot(state)
            % Rebuild the saved reconstruction without a profile or solve.
            %
            % - Topic: Restore vertical calculus
            % - Declaration: calculus = IMLobattoCalculus.fromSnapshot(state)
            % - Parameter state: numeric data from snapshot
            % - Returns self: calculus on the saved grid
            arguments
                state (1,1) struct
            end
            fields = ["version","z","nativeCoordinate","metricBreaks","metricCoefficientCounts", ...
                "metricCoefficients","metricPointValues"];
            if ~all(isfield(state,fields)) || ~isequal(state.version,1)
                error("IMLobattoCalculus:InvalidSnapshot","A supported mapped Lobatto calculus snapshot is required.");
            end
            breaks = state.metricBreaks(:).';
            counts = state.metricCoefficientCounts(:);
            coefficients = state.metricCoefficients(:);
            pointValues = state.metricPointValues(:);
            if numel(breaks) ~= numel(counts)+1 || isempty(counts) || ...
                    ~isreal(breaks) || any(~isfinite(breaks)) || any(diff(breaks) <= 0) || ...
                    ~isreal(counts) || any(~isfinite(counts) | counts < 1 | counts ~= fix(counts)) || ...
                    sum(counts) ~= numel(coefficients) || ~isreal(coefficients) || any(~isfinite(coefficients)) || ...
                    numel(pointValues) ~= numel(breaks) || ~isreal(pointValues) || any(~isfinite(pointValues) | pointValues <= 0)
                error("IMLobattoCalculus:InvalidSnapshot","The metric piece boundaries and coefficients are inconsistent.");
            end
            slices = mat2cell(coefficients,counts,1);
            pieces = cell(numel(counts),1);
            for j = 1:numel(counts)
                piece = chebfun(slices{j},breaks(j:j+1),'coeffs');
                pieces{j} = piece.funs{1};
            end
            metric = chebfun(pieces);
            % Chebfun can preserve independently sampled breakpoint values.
            % Keep the saved values so nodal derivatives round-trip exactly.
            if max(abs(metric.pointValues-pointValues)) > 1e-10*max(pointValues)
                error("IMLobattoCalculus:InvalidSnapshot","The metric point values disagree with the represented pieces.");
            end
            metric.pointValues = pointValues;
            self = IMLobattoCalculus.fromMetric(metric,state.z(:),state.nativeCoordinate(:));
        end
    end

    methods (Static, Hidden)
        function self = fromMetric(metric,z,x)
            % Bind an existing set of nodes to the metric that made it.
            map = IMLobattoCalculus.integratedMetric(metric);
            self = IMLobattoCalculus(metric,map,z,x);
        end

        function [self,z,weights] = forMetric(metric,nPoints)
            % Prepare a Lobatto rule directly from a positive metric.
            map = IMLobattoCalculus.integratedMetric(metric);
            [x,w] = chebpts(nPoints,[0 feval(map,map.domain(end))]);
            z = IMLobattoCalculus.invertMap(map,x);
            weights = w(:)./feval(metric,z);
            weights = diff(z([1 end]))*weights/sum(weights);
            self = IMLobattoCalculus(metric,map,z,x);
        end
    end

    methods (Access = private)
        function self = IMLobattoCalculus(metric,map,z,x)
            domain = metric.domain([1 end]);
            nativeDomain = [0 feval(map,domain(2))];
            nPoints = numel(z);
            if nPoints < 3 || numel(x) ~= nPoints || ~isreal(z) || ~isreal(x) || ...
                    any(~isfinite(z)) || any(~isfinite(x)) || any(diff(z) <= 0) || any(diff(x) <= 0) || ...
                    ~isequal(z([1 end]).',domain) || any(~isfinite(nativeDomain)) || nativeDomain(2) <= 0
                error("IMLobattoCalculus:InvalidSnapshot","Lobatto samples and map domain are inconsistent.");
            end
            [expected,~,barycentricWeights] = chebpts(nPoints,nativeDomain);
            mapError = max(abs(feval(map,z)-x));
            if max(abs(x-expected)) > 64*eps(nativeDomain(2)) || mapError > 1e-10*nativeDomain(2)
                error("IMLobattoCalculus:InvalidSnapshot","Samples must be the mapped Chebyshev-Lobatto nodes.");
            end
            self.metric = metric;
            self.map = map;
            self.z = z;
            self.nativeCoordinate = x;
            self.barycentricWeights = barycentricWeights(:);
        end

        function Dz = derivativeAtNodes(self,z,x)
            Dz = feval(self.metric,z).*diffmat(numel(z),1,x([1 end]).');
        end

        function P = interpolationAtNativePoints(self,x)
            [atNode,index] = ismember(x,self.nativeCoordinate);
            P = zeros(numel(x),numel(self.z));
            P(sub2ind(size(P),find(atNode),index(atNode))) = 1;
            terms = self.barycentricWeights.'./(x(~atNode)-self.nativeCoordinate.');
            P(~atNode,:) = terms./sum(terms,2);
        end

        function validateSamples(self,samples)
            if size(samples,1) ~= numel(self.z)
                error("IMLobattoCalculus:SampleMismatch","Supply one sample row per Lobatto point.");
            end
        end
    end

    methods (Static, Access = private)
        function map = integratedMetric(metric)
            if ~isreal(metric)
                error("IMLobattoCalculus:InvalidMetric","The represented coordinate metric must be real.");
            end
            minimum = min(metric);
            if ~isfinite(minimum) || minimum <= 0
                error("IMLobattoCalculus:InvalidMetric","The represented coordinate metric must be finite and strictly positive.");
            end
            map = cumsum(metric);
            map = map-feval(map,metric.domain(1));
        end

        function z = invertMap(map,x)
            if isempty(x), z = x; return; end
            domain = map.domain([1 end]);
            lower = repmat(domain(1),size(x));
            upper = repmat(domain(2),size(x));
            for iteration = 1:55
                middle = (lower+upper)/2;
                below = feval(map,middle) < x;
                lower(below) = middle(below);
                upper(~below) = middle(~below);
            end
            z = (lower+upper)/2;
            z(x == 0) = domain(1);
            z(x == feval(map,domain(2))) = domain(2);
        end
    end
end
