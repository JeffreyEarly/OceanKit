classdef BSpline < handle
    % Create, evaluate, and manipulate one-dimensional terminated B-splines.
    %
    % `BSpline` is the low-level one-dimensional spline object used by the
    % higher-level interpolation and fitting classes. It stores a spline
    % degree `S`, a terminated knot sequence `knotPoints`, and a coefficient
    % vector `xi`, then caches an equivalent piecewise-polynomial
    % representation for fast evaluation.
    %
    % Mathematically, the stored spline is
    %
    % $$
    % f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau),
    % $$
    %
    % where $$\tau$$ is the terminated knot sequence, $$B_{j,S}$$ are the
    % one-dimensional B-spline basis functions of degree $$S$$, and
    % `xMean` is added only for zero-order evaluation.
    %
    % ## Basic usage
    %
    % In most workflows you first build a knot sequence from sample
    % locations, assemble the spline basis matrix, solve for coefficients,
    % and then evaluate the resulting spline object.
    %
    % ```matlab
    % t = linspace(0,1,20)';
    % x = sin(2*pi*t);
    % knotPoints = BSpline.knotPointsForDataPoints(t, S=3);
    % X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3);
    % spline = BSpline(S=3, knotPoints=knotPoints, xi=X\x);
    %
    % xq = spline(linspace(0,1,100)');
    % ```
    %
    % - Topic: Create a spline
    % - Topic: Inspect spline properties
    % - Topic: Evaluate the spline
    % - Topic: Transform the spline
    % - Topic: Build spline bases
    % - Topic: Represent piecewise polynomials
    % - Topic: Maintain cached state
    %
    % - Declaration: classdef BSpline < handle
    properties (SetAccess = private)
        % Spline order K, where polynomial degree is S = K - 1.
        %
        % The order is the number of coefficients in each local polynomial
        % piece. On any one interval, a spline of order `K` is a polynomial
        % of the form
        %
        % $$
        % p_i(t) = a_{i,0} + a_{i,1}(t-t_i) + \cdots + a_{i,K-1}(t-t_i)^{K-1}.
        % $$
        %
        % So `K=1` gives piecewise constants, `K=2` gives piecewise
        % linear splines, and `K=4` gives cubic splines. The matching degree
        % property is [`S`](/spline-core/classes/bspline/s.html), with
        % `S = K - 1`.
        %
        % ```matlab
        % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
        % spline.K
        % % returns 4, meaning each piece is cubic
        % ```
        %
        % - Topic: Inspect spline properties
        K
    end

    properties (Access = private)
        % Internal spline coefficients, stored as an Mx1 vector.
        %
        % - Topic: Inspect spline properties
        xi_

        % Internal knot sequence for the spline basis.
        %
        % - Topic: Inspect spline properties
        tKnot_
    end

    properties (GetAccess=public, SetAccess=protected)
        % Basis values and derivatives sampled at piecewise breakpoints.
        %
        % This cached array stores the basis and its derivatives evaluated at
        % the piecewise-polynomial breakpoints [`t_pp`](/spline-core/classes/bspline/t_pp.html):
        %
        % $$
        % \mathrm{Xtpp}(i,j,d+1) = B_{j,S}^{(d)}(t_{\mathrm{pp},i};\tau).
        % $$
        %
        % [`ppCoefficientsFromSplineCoefficients`](/spline-core/classes/bspline/ppcoefficientsfromsplinecoefficients.html)
        % uses `Xtpp` to convert the spline coefficients `xi` into the cached
        % interval coefficients [`C`](/spline-core/classes/bspline/c.html),
        % and [`evaluateFromPPCoefficients`](/spline-core/classes/bspline/evaluatefromppcoefficients.html)
        % consumes the resulting PP representation for fast evaluation.
        %
        % ```matlab
        % [C, tpp, Xtpp] = BSpline.ppCoefficientsFromSplineCoefficients( ...
        %     xi=spline.xi, knotPoints=spline.knotPoints, S=spline.S);
        % size(Xtpp)
        % % numel(tpp) x numel(xi) x (S+1)
        % ```
        %
        % - Topic: Represent piecewise polynomials
        % - Developer: true
        Xtpp = [];

        % Piecewise-polynomial breakpoint locations.
        %
        % These are the breakpoints of the cached interval representation.
        % If `Nk = numel(knotPoints)`, then
        %
        % $$
        % t_{\mathrm{pp}} = \tau_{K:(N_k-K+1)}.
        % $$
        %
        % On interval `i`, the PP cache works in the local coordinate
        % `u = t - t_pp(i)` for `t \in [t_pp(i), t_pp(i+1)]`.
        % The functions
        % [`ppCoefficientsFromSplineCoefficients`](/spline-core/classes/bspline/ppcoefficientsfromsplinecoefficients.html)
        % and [`evaluateFromPPCoefficients`](/spline-core/classes/bspline/evaluatefromppcoefficients.html)
        % are the main consumers of this breakpoint vector.
        %
        % ```matlab
        % [C, tpp] = BSpline.ppCoefficientsFromSplineCoefficients( ...
        %     xi=spline.xi, knotPoints=spline.knotPoints, S=spline.S);
        % values = BSpline.evaluateFromPPCoefficients( ...
        %     queryPoints=tQuery, C=C, tpp=tpp);
        % ```
        %
        % - Topic: Represent piecewise polynomials
        % - Developer: true
        % size(t_pp) = length(knotPoints) - 2*K + 1
        t_pp

        % Piecewise-polynomial coefficients for interval evaluation.
        %
        % For interval `i`, let `u = t - t_pp(i)`. The cached PP form stores
        %
        % $$
        % f_i(u) = \sum_{m=0}^{S} \frac{c_{i,m}}{m!} u^m,
        % $$
        %
        % where row `i` of `C` contains the coefficients in descending power
        % order so they can be passed to `polyval` after factorial scaling.
        % The helper
        % [`ppCoefficientsFromSplineCoefficients`](/spline-core/classes/bspline/ppcoefficientsfromsplinecoefficients.html)
        % builds `C`, and
        % [`evaluateFromPPCoefficients`](/spline-core/classes/bspline/evaluatefromppcoefficients.html)
        % evaluates it.
        %
        % ```matlab
        % [C, tpp] = BSpline.ppCoefficientsFromSplineCoefficients( ...
        %     xi=spline.xi, knotPoints=spline.knotPoints, S=spline.S);
        % xq = BSpline.evaluateFromPPCoefficients(queryPoints=tQuery, C=C, tpp=tpp);
        % ```
        %
        % - Topic: Represent piecewise polynomials
        % - Developer: true
        % size(C) = [length(t_pp)-1, K]
        C       
    end
    
    properties (SetAccess = private)
        % Mean added back to zero-order spline evaluations.
        %
        % `xMean` is an output offset used in the stored spline model
        %
        % $$
        % f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau).
        % $$
        %
        % This affine output normalization is useful for numerical work:
        % large means can be removed before solving for the coefficient
        % vector, then added back only at evaluation time. That keeps the
        % fitted coefficients closer to order one and reduces the need to
        % represent a large constant level directly in `xi`.
        %
        % `xMean` contributes only to zero-order evaluation. Derivatives are
        % unaffected because constants differentiate to zero.
        %
        % ```matlab
        % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi, xMean=12.4, xStd=0.8);
        % values = spline(tQuery);
        % ```
        %
        % - Topic: Inspect spline properties
        xMean (1,1) double {mustBeReal,mustBeFinite} = 0
        % Multiplicative scale applied to spline evaluations.
        %
        % `xStd` is the multiplicative scaling in
        %
        % $$
        % f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau).
        % $$
        %
        % It exists for the same numerical reason as
        % [`xMean`](/spline-core/classes/bspline/xmean.html): when the output
        % amplitude is large or very small, solving for a normalized
        % coefficient vector `xi` is often better conditioned than solving
        % directly in physical units. Unlike `xMean`, `xStd` multiplies both
        % the spline values and all derivatives.
        %
        % ```matlab
        % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi, xMean=12.4, xStd=0.8);
        % d2values = spline.valueAtPoints(tQuery, D=2);
        % ```
        %
        % - Topic: Inspect spline properties
        xStd (1,1) double {mustBeReal,mustBeFinite} = 1
    end

    properties (Dependent)
        % Polynomial degree S = K - 1.
        %
        % The degree is the highest power that appears in each local
        % polynomial piece:
        %
        % $$
        % p_i(t) = a_{i,0} + a_{i,1}(t-t_i) + \cdots + a_{i,S}(t-t_i)^S.
        % $$
        %
        % Degree `S=0` is piecewise constant, `S=1` is piecewise linear,
        % `S=2` is quadratic, and `S=3` is cubic. The matching order is
        % [`K`](/spline-core/classes/bspline/k.html), with `K = S + 1`.
        %
        % ```matlab
        % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
        % spline.S
        % % returns 3 for a cubic spline
        % ```
        %
        % - Topic: Inspect spline properties
        % A cubic spline is K=4, S=3
        S

        % Minimum and maximum values of the spline domain.
        %
        % For a terminated spline basis, the domain is the interval covered
        % by the repeated end knots:
        %
        % $$
        % \mathrm{domain} = [\tau_1,\ \tau_{N_k}].
        % $$
        %
        % Outside this interval the terminated basis is zero, so this
        % property is the natural plotting and evaluation range for the
        % spline.
        %
        % ```matlab
        % spline.domain
        % % returns [tMin tMax]
        % ```
        %
        % - Topic: Inspect spline properties
        domain

        % Spline coefficients as an Mx1 vector.
        %
        % The coefficient vector weights the terminated B-spline basis in
        %
        % $$
        % f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau).
        % $$
        %
        % So `xi(j)` is the weight on the `j`th basis function. The basis
        % itself comes from
        % [`matrixForDataPoints`](/spline-core/classes/bspline/matrixfordatapoints.html),
        % and evaluation is handled by
        % [`valueAtPoints`](/spline-core/classes/bspline/valueatpoints.html).
        % For a knot sequence `tau` and order `K`, the coefficient count is
        % `M = numel(knotPoints) - K`.
        %
        % ```matlab
        % X = BSpline.matrixForDataPoints(t, knotPoints=knotPoints, S=3);
        % xi = X \ x;
        % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
        % ```
        %
        % - Topic: Inspect spline properties
        xi
    end

    properties (Dependent, SetAccess = private)

        % Knot sequence used to define the spline basis.
        %
        % This is the terminated knot vector
        %
        % $$
        % \tau = [\tau_1,\ldots,\tau_{N_k}]^\mathsf{T}
        % $$
        %
        % that defines the basis functions $$B_{j,S}(t;\tau)$$. Repeating the
        % first and last knot values `K` times terminates the basis at the
        % endpoints; interior multiplicity controls continuity across knot
        % locations.
        %
        % For example, the cubic knot vector
        %
        % ```matlab
        % knotPoints = [0; 0; 0; 0; 0.3; 0.7; 1; 1; 1; 1];
        % ```
        %
        % defines four cubic basis functions on the domain `[0, 1]`. Use
        % [`knotPointsForDataPoints`](/spline-core/classes/bspline/knotpointsfordatapoints.html)
        % to build a terminated knot sequence from sample locations.
        %
        % - Topic: Inspect spline properties
        knotPoints
    end
    
    methods
        varargout = subsref(self, index)
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = BSpline(options)
            % Create a one-dimensional spline from degree, knots, and coefficients.
            %
            % Use this constructor when you already know the terminated knot
            % sequence `knotPoints` and the coefficient vector `xi`.
            %
            % The constructed spline is
            %
            % $$
            % f(t) = x_{\mathrm{Mean}} + x_{\mathrm{Std}} \sum_{j=1}^{M} \xi_j B_{j,S}(t;\tau).
            % $$
            %
            % ```matlab
            % knotPoints = [0; 0; 0; 0; 1; 1; 1; 1];
            % xi = [1; -0.5; 0.25; 0];
            % spline = BSpline(S=3, knotPoints=knotPoints, xi=xi);
            % x = spline(linspace(0,1,50)');
            % ```
            %
            % - Topic: Create a spline
            % - Declaration: spline = BSpline(options)
            % - Parameter options.S: spline degree
            % - Parameter options.knotPoints: knot points
            % - Parameter options.xi: (optional) spline coefficients
            % - Parameter options.Xtpp: optional cached basis values at piecewise breakpoints
            % - Parameter options.xMean: optional additive output offset
            % - Parameter options.xStd: optional multiplicative output scale
            % - Returns spline: BSpline instance
            arguments
                options.S (1,1) double {mustBeInteger,mustBeNonnegative}
                options.knotPoints (:,1) double {mustBeNumeric,mustBeReal}
                options.xi {mustBeNumeric,mustBeReal,mustBeFinite} = []
                options.Xtpp (:,:,:) double = []
                options.xMean (1,1) double {mustBeReal,mustBeFinite} = 0
                options.xStd (1,1) double {mustBeReal,mustBeFinite} = 1
            end

            K = options.S + 1;
            knotPoints = options.knotPoints;
            if isempty(options.xi)
                xi = zeros(length(knotPoints) - K, 1);
            else
                xi = reshape(options.xi, [], 1);
            end

            self.K = K;
            self.tKnot_ = knotPoints;
            self.Xtpp = options.Xtpp;
            self.xMean = options.xMean;
            self.xStd = options.xStd;
            self.xi = xi;
        end
        
        function S = get.S(self)
            % Return the spline polynomial degree.
            %
            % - Topic: Inspect spline properties
            % - Declaration: S = get.S(self)
            % - Parameter self: BSpline instance
            % - Returns S: double scalar equal to K - 1
            S = self.K-1;
        end
        
        function domain = get.domain(self)
            % Return the spline domain endpoints.
            %
            % - Topic: Inspect spline properties
            % - Declaration: domain = get.domain(self)
            % - Parameter self: BSpline instance
            % - Returns domain: 1x2 vector [tMin tMax]
            domain = [self.knotPoints(1) self.knotPoints(end)];
        end

        function xi = get.xi(self)
            % Return the current spline coefficients.
            %
            % - Topic: Inspect spline properties
            % - Declaration: xi = get.xi(self)
            % - Parameter self: BSpline instance
            % - Returns xi: spline coefficient column vector
            xi = self.xi_;
        end

        function set.xi(self, xi)
            % Update spline coefficients and refresh cached polynomial forms.
            %
            % - Topic: Inspect spline properties
            % - Declaration: set.xi(self,xi)
            % - Parameter self: BSpline instance
            % - Parameter xi: spline coefficient column vector
            arguments
                self (1,1) BSpline
                xi {mustBeNumeric,mustBeReal,mustBeFinite}
            end
            xi = reshape(xi, [], 1);
            if ~isempty(xi) && numel(xi) ~= numel(self.tKnot_) - self.K
                error('BSpline:InvalidCoefficientCount', 'xi must contain exactly numel(knotPoints) - K coefficients.');
            end
            self.xi_ = xi;
            self.splineCoefficientsDidChange();
        end

        function knotPoints = get.knotPoints(self)
            % Return the current knot sequence.
            %
            % - Topic: Inspect spline properties
            % - Declaration: knotPoints = get.knotPoints(self)
            % - Parameter self: BSpline instance
            % - Returns knotPoints: knot vector
            knotPoints = self.tKnot_;
        end
        
        x_out = valueAtPoints(self, t, options)
        tKnotDidChange(self)
        splineCoefficientsDidChange(self)
    end

    methods (Access = private)
        transformedSpline = affineOutputTransform(self, scale, offset)
    end
    
    methods (Static)
        knotPoints = knotPointsForDataPoints(dataPoints, options)
        t = pointsOfSupportFromKnotPoints(knotPoints, options)
        [C,tpp,Xtpp] = ppCoefficientsFromSplineCoefficients(options)
        f = evaluateFromPPCoefficients(options)
        B = matrixForDataPoints(dataPoints, options)
    end
end
