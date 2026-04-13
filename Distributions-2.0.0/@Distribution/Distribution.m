classdef (Abstract) Distribution < handle & matlab.mixin.Heterogeneous & CAAnnotatedClass
    % Represent a univariate probability distribution and related noise model.
    %
    % `Distribution` is the abstract base class for a univariate
    % [probability distribution](https://en.wikipedia.org/wiki/Probability_distribution)
    % with density $$p(z)$$ and cumulative distribution
    % $$F(z) = \int_{-\infty}^{z} p(\zeta)\,\mathrm{d}\zeta.$$ Subclasses
    % provide those distribution functions, and this base class adds
    % percentile inversion, range-limited second-moment integrals, random
    % sampling, correlated-noise generation, and goodness-of-fit
    % diagnostics.
    %
    % The method `locationOfCDFPercentile` returns $$z_{\alpha}$$ satisfying
    % $$F(z_{\alpha}) = \alpha.$$ The method `varianceInRange` evaluates
    % $$\int_{z_{\min}}^{z_{\max}} z^{2} p(z)\,\mathrm{d}z,$$ which
    % coincides with the variance contribution on that interval for
    % centered distributions. The Anderson-Darling and Kolmogorov-Smirnov
    % methods compare the empirical distribution of a sample against the
    % theoretical CDF $$F(z)$$ through standard AD and KS statistics.
    %
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Declaration: classdef (Abstract) Distribution
    
    properties
        % Autocorrelation function $$\rho(\tau)$$ used by `noise`.
        %
        % When nonempty, `rho` is evaluated on the lag grid passed to
        % `noise`, and the resulting Toeplitz correlation model is used to
        % correlate otherwise independent draws.
        %
        % - Topic: Model correlated noise
        rho

        % Weight function $$w(z)$$ for residual-based reweighting.
        %
        % Robust-fitting code uses `w` as a residual-dependent variance
        % proxy. When defined from the density, a common form is
        % $$w(z) = -\frac{p(z)}{(\partial_z p(z))/z}.$$
        %
        % - Topic: Inspect distribution properties
        w

        % Logarithm of the probability density $$\log p(z)$$.
        %
        % `logPDF` can be evaluated directly in the tails without first
        % forming `pdf(z)`, which is often more stable than computing
        % `log(pdf(z))`.
        %
        % - Topic: Inspect distribution properties
        logPDF
    end

    properties (SetAccess = protected)
        % Probability density function $$p(z)$$.
        %
        % `pdf` stores a function handle that evaluates the density at one
        % or more query locations. A valid density integrates to one across
        % the support of the distribution.
        %
        % - Topic: Inspect distribution properties
        pdf

        % Cumulative distribution function $$F(z)$$.
        %
        % `cdf` stores a function handle that evaluates
        % $$F(z) = \int_{-\infty}^{z} p(\zeta)\,\mathrm{d}\zeta$$ at one or
        % more query locations.
        %
        % - Topic: Inspect distribution properties
        cdf

        % Support interval $$[z_{\min}, z_{\max}]$$ for the distribution.
        %
        % `zrange` is used to clamp integration bounds for methods such as
        % `varianceInRange`. The default support is
        % $$(-\infty,\infty).$$
        %
        % - Topic: Inspect distribution properties
        zrange = [-Inf Inf] % range of support

        % Total variance.
        %
        % Each concrete distribution sets `variance` according to its
        % mathematical family. For centered distributions this also matches
        % the full-range integral from `varianceInRange`.
        %
        % - Topic: Inspect distribution properties
        variance % total variance
    end

    properties (Dependent)
        % Initialization scale $$\sigma_0$$ satisfying $$w(\sigma_0) = \sigma_0^2$$.
        %
        % `sigma0` is computed from the current weight function and is used
        % to seed iterative reweighting methods in downstream packages.
        %
        % - Topic: Inspect distribution properties
        sigma0 % initial seed for weight function (set by default to w(sigma0)=sigma0^2)
    end

    properties (Access = protected)
        dPDFOverZ % derivative of the pdf wrt z, divided by z
        logPDFNorm
    end
    
    
    methods (Access = protected)
        function self = Distribution()
            % Initialize base storage for a distribution subclass.
            %
            % Concrete subclasses call this constructor implicitly before
            % they set the distribution function handles and cached
            % statistics stored by the base class.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: self = Distribution()
            % - Returns self: Distribution base-class instance
            self@CAAnnotatedClass();
        end

        function restoreOptionalPersistedPropertiesFromGroup(self,group)
            arguments
                self (1,1) Distribution
                group (1,1) NetCDFGroup
            end
            vars = CAAnnotatedClass.propertyValuesFromGroup(group, {'rho'}, shouldIgnoreMissingProperties=true);
            if ~isempty(vars) && isfield(vars, 'rho')
                self.rho = vars.rho;
            end
        end
    end

    methods
        function tf = isequal(self,other)
            arguments
                self (1,1) Distribution
                other (1,1) Distribution
            end
            tf = isequal@CAAnnotatedClass(self, other);
            if tf
                if isempty(self.rho) || isempty(other.rho)
                    tf = isempty(self.rho) && isempty(other.rho);
                else
                    tf = isequal(func2str(self.rho), func2str(other.rho));
                end
            end
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            % Write this distribution into an existing NetCDF group.
            %
            % Nested persistence paths use `writeToGroup(...)` rather than
            % `writeToFile(...)`, so this override ensures an optional
            % nonempty `rho` is preserved there as well.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: writeToGroup(self,group,propertyAnnotations,attributes)
            % - Parameter group: target NetCDFGroup
            % - Parameter propertyAnnotations: property annotations selected for persistence
            % - Parameter attributes: optional group attributes
            arguments
                self (1,1) Distribution
                group (1,1) NetCDFGroup
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end

            selectedNames = string({propertyAnnotations.name});
            if ~isempty(self.rho) && ~any(selectedNames == "rho")
                propertyAnnotations = cat(2, propertyAnnotations, self.propertyAnnotationWithName("rho"));
            end
            writeToGroup@CAAnnotatedClass(self, group, propertyAnnotations, attributes);
        end

        function ncfile = writeToFile(self,path,properties,options)
            % Write this distribution to a NetCDF file.
            %
            % `writeToFile` persists the canonical constructor state needed
            % to reconstruct the concrete distribution subtype. When `rho`
            % is nonempty, it is also serialized through the annotated
            % function-handle path.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: ncfile = writeToFile(self,path,properties,options)
            % - Parameter path: path to the NetCDF file to create
            % - Parameter properties: optional additional property names to persist
            % - Parameter options.shouldOverwriteExisting: optional logical scalar that overwrites an existing file when true
            % - Parameter options.attributes: optional global NetCDF attributes to add to the file
            % - Returns ncfile: NetCDFFile handle for the written file
            arguments (Input)
                self (1,1) Distribution
                path char {mustBeNonempty}
            end
            arguments (Input,Repeating)
                properties char
            end
            arguments (Input)
                options.shouldOverwriteExisting logical = false
                options.attributes = configureDictionary("string","string")
            end
            arguments (Output)
                ncfile NetCDFFile
            end

            propertyNames = union(string(properties), string(self.requiredProperties()));
            if ~isempty(self.rho)
                propertyNames = union(propertyNames, "rho");
            end
            propertyAnnotations = self.propertyAnnotationWithName(propertyNames);
            selectedPropertyNames = string({propertyAnnotations.name});
            dimensionNames = string.empty(0, 1);
            for iProperty = 1:length(propertyAnnotations)
                if isa(propertyAnnotations(iProperty), 'CANumericProperty')
                    dimensionNames = union(dimensionNames, string(propertyAnnotations(iProperty).dimensions));
                end
            end
            dimensionNames = setdiff(dimensionNames, selectedPropertyNames);
            if ~isempty(dimensionNames)
                propertyAnnotations = cat(2, self.propertyAnnotationWithName(dimensionNames), propertyAnnotations);
            end
            ncfile = CAAnnotatedClass.writeToPath(self, path, shouldOverwriteExisting=options.shouldOverwriteExisting, propertyAnnotations=propertyAnnotations, attributes=options.attributes);
            ncfile.sync();
        end

        function z = locationOfCDFPercentile(self, alpha)
            % Find the location of a CDF percentile.
            %
            % This method returns $$z_{\alpha}$$ satisfying
            % $$F(z_{\alpha}) = \alpha$$ for the current cumulative
            % distribution $$F$$.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: z = locationOfCDFPercentile(self,alpha)
            % - Parameter alpha: target percentile strictly between 0 and 1
            % - Returns z: location where `cdf(z)` equals `alpha`
            arguments
                self (1,1) Distribution
                alpha (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeGreaterThan(alpha,0),mustBeLessThan(alpha,1)}
            end
            z = fminsearch(@(z) abs(self.cdf(z)-alpha),0);
        end
        
        function var = varianceInRange(self,zmin,zmax)
            % Compute the range-limited second moment of the distribution.
            %
            % This method evaluates
            % $$\int_{z_{\min}}^{z_{\max}} z^{2} p(z)\,\mathrm{d}z.$$ For
            % centered distributions, that integral equals the variance
            % contribution on the interval.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: var = varianceInRange(self,zmin,zmax)
            % - Parameter zmin: lower integration bound
            % - Parameter zmax: upper integration bound
            % - Returns var: range-limited second moment
            arguments
                self (1,1) Distribution
                zmin (1,1) {mustBeNumeric,mustBeReal}
                zmax (1,1) {mustBeNumeric,mustBeReal}
            end
            zmin = Distribution.truncateToSupport(zmin,self.zrange);
            zmax = Distribution.truncateToSupport(zmax,self.zrange);
            var = integral( @(z) z.*z.*self.pdf(z),zmin,zmax);
        end
        
        function var = varianceInPercentileRange(self,pctmin,pctmax)
            % Compute the second moment between two CDF percentiles.
            %
            % This method first finds $$z_{\alpha_{\min}}$$ and
            % $$z_{\alpha_{\max}}$$ from the inverse CDF and then evaluates
            % `varianceInRange` between those locations.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: var = varianceInPercentileRange(self,pctmin,pctmax)
            % - Parameter pctmin: lower percentile strictly between 0 and 1
            % - Parameter pctmax: upper percentile strictly between 0 and 1
            % - Returns var: second moment between the requested percentiles
            arguments
                self (1,1) Distribution
                pctmin (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeGreaterThan(pctmin,0),mustBeLessThan(pctmin,1)}
                pctmax (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeGreaterThan(pctmax,0),mustBeLessThan(pctmax,1)}
            end
            zmin = self.locationOfCDFPercentile(pctmin);
            zmax = self.locationOfCDFPercentile(pctmax);
            var = self.varianceInRange(zmin,zmax);
        end
        
        % initial weighting should start at the inflection point of the
        % weight function
        function sigma0 = get.sigma0(self)
            sigma0 = fminsearch( @(x) abs(self.w(x)-x^2),sqrt(self.variance));
        end
        
        function totalError = andersonDarlingError(self,epsilon)
            % Compute the Anderson-Darling fit statistic for a sample.
            %
            % Given the sorted sample $$Y_i$$, this method evaluates the
            % Anderson-Darling statistic formed from the theoretical CDF
            % $$F$$ and the empirical order statistics.
            %
            % - Topic: Evaluate distribution fit
            % - Declaration: totalError = andersonDarlingError(self,epsilon)
            % - Parameter epsilon: sample values to compare against the distribution
            % - Returns totalError: Anderson-Darling discrepancy statistic
            arguments
                self (1,1) Distribution
                epsilon {mustBeNumeric,mustBeReal}
            end
            Y = sort(epsilon);
            n = length(Y);
            
            s = ((2*(1:n)'-1)/n) .* (log(self.cdf(Y)) + log(1-self.cdf(flip(Y))));
            
            totalError = -n-sum(s);
        end
        
        function totalError = andersonDarlingInterquartileError(self,epsilon)
            % Compute an interquartile Anderson-Darling fit statistic.
            %
            % This method forms the Anderson-Darling contribution sequence
            % and retains only the middle half of the ordered sample before
            % summing the error.
            %
            % - Topic: Evaluate distribution fit
            % - Declaration: totalError = andersonDarlingInterquartileError(self,epsilon)
            % - Parameter epsilon: sample values to compare against the distribution
            % - Returns totalError: interquartile Anderson-Darling discrepancy statistic
            arguments
                self (1,1) Distribution
                epsilon {mustBeNumeric,mustBeReal}
            end
            if issorted(epsilon)
                Y = epsilon;
            else
                Y = sort(epsilon);
            end
            n = length(Y);
            
            s = ((2*(1:n)'-1)/n) .* (log(self.cdf(Y)) + log(1-self.cdf(flip(Y))));
            sIQ = s(floor(n/4):ceil(3*n/4));
            
            totalError = length(sIQ)-sum(sIQ);
        end
        
        function totalError = kolmogorovSmirnovError(self,epsilon,zmin,zmax)
            % Compute the Kolmogorov-Smirnov fit statistic for a sample.
            %
            % This method evaluates the scaled KS discrepancy
            % $$\left(\sqrt{n} + 0.12 + \frac{0.11}{\sqrt{n}}\right) D_n,$$
            % where $$D_n$$ is the maximum separation between the empirical
            % CDF and the theoretical CDF. When `zmin` and `zmax` are
            % supplied, the comparison is restricted to that interval and
            % the theoretical CDF is renormalized there.
            %
            % - Topic: Evaluate distribution fit
            % - Declaration: totalError = kolmogorovSmirnovError(self,epsilon,zmin,zmax)
            % - Parameter epsilon: sample values to compare against the distribution
            % - Parameter zmin: optional lower bound for a truncated comparison interval
            % - Parameter zmax: optional upper bound for a truncated comparison interval
            % - Returns totalError: Kolmogorov-Smirnov discrepancy statistic
            arguments
                self (1,1) Distribution
                epsilon {mustBeNumeric,mustBeReal}
                zmin {mustBeNumeric,mustBeReal,mustBeScalarOrEmpty} = []
                zmax {mustBeNumeric,mustBeReal,mustBeScalarOrEmpty} = []
            end

            if nargin == 4
                x = sort(epsilon( epsilon > zmin & epsilon < zmax ));
                A = 1/(self.cdf(zmax)-self.cdf(zmin));
                b = self.cdf(zmin);
            else
                if issorted(epsilon)
                    x = epsilon;
                else
                    x = sort(epsilon);
                end
                A = 1;
                b = 0;
            end
            
            n = length(x);
            y_data = (1:n)'/n;
            y = A*(self.cdf(x)-b);
            
            D = max(abs(y-y_data));
            
            totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
            if isempty(totalError)
                totalError = Inf;
            end
        end
        
        function totalError = kolmogorovSmirnovInterquartileError(self,epsilon)
            % Compute an interquartile Kolmogorov-Smirnov fit statistic.
            %
            % This method evaluates the empirical-versus-theoretical CDF
            % mismatch on the middle half of the sorted sample and returns
            % the scaled interquartile KS discrepancy.
            %
            % - Topic: Evaluate distribution fit
            % - Declaration: totalError = kolmogorovSmirnovInterquartileError(self,epsilon)
            % - Parameter epsilon: sample values to compare against the distribution
            % - Returns totalError: interquartile Kolmogorov-Smirnov discrepancy statistic
            arguments
                self (1,1) Distribution
                epsilon {mustBeNumeric,mustBeReal}
            end
            x = sort(epsilon);
            n = length(x);
            y_data = (1:n)'/n;
            y = self.cdf(x);
            
            ks = y-y_data;
            DIQ = max(abs(ks(floor(n/4):ceil(3*n/4))));
            n = length(DIQ);
            
            totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*DIQ;
        end
        
        function y = rand(self,varargin)
            % Draw random samples from the distribution.
            %
            % `rand` draws samples by mapping uniform random numbers
            % through a discretized approximation to the inverse CDF. Pass
            % either one size vector or one nonnegative integer per output
            % dimension.
            %
            % - Topic: Sample from distributions
            % - Declaration: y = rand(self,varargin)
            % - Parameter varargin: size vector or one nonnegative integer per output dimension
            % - Returns y: random samples with the requested output shape
            arguments
                self (1,1) Distribution
            end
            arguments (Repeating)
                varargin {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
            end
            sz = Distribution.sizeVectorFromInputs(varargin{:});
            n = prod(sz);
            pct = 1/1e6;
            zmin = self.locationOfCDFPercentile(pct/2);
            zmax = self.locationOfCDFPercentile(1-pct/2);
            binEdges = linspace(zmin,zmax,1e6+1)';
            binWidths = diff(binEdges);
            binEdges_cdf = self.cdf(binEdges); % maps the bin edges [0, 1]
            [~, ~, bin] = histcounts(rand(n+10,1),binEdges_cdf);
            bin(bin==0 | bin == length(binWidths)+1) = [];
            y = binEdges(bin) + rand(length(bin),1).*binWidths(bin);
            y = y(1:n);
            y = reshape(y,sz);
        end
        
        function y = noise(self,t)
            % Draw a correlated noise process on a lag grid.
            %
            % When `rho` is empty, `noise` returns independent samples with
            % the same shape as `t`. When `rho` is defined, it evaluates
            % $$\rho(t)$$, forms the corresponding Toeplitz correlation
            % model, and applies its Cholesky factor to the independent
            % draws.
            %
            % - Topic: Model correlated noise
            % - Declaration: y = noise(self,t)
            % - Parameter t: lag grid used to evaluate the autocorrelation model
            % - Returns y: random process samples with the same shape as `t`
            arguments
                self (1,1) Distribution
                t {mustBeNumeric,mustBeReal}
            end
            y = self.rand(size(t));
            
            % Some notes:
            % For a given autocorrelation function, the structure of C away
            % from the end points is uniform for each row. Thus, it's
            % possible to find this structure and convolve it as a kernal
            % across the random increments. That'd be much faster/less
            % memory.
            if ~isempty(self.rho)
                r = self.rho(t);
                C = toeplitz(r);
                
                % this is from maternchol.m and helps keep the matrix
                % positive definite.
                eps=1e-12;
                C=(C+eye(size(C))*eps)./(1+eps);
                
                
                T = chol(C,'lower');
                y = T*y;
            end
        end
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAFunctionProperty('rho', 'Autocorrelation function $$\\rho(\\tau)$$ used to generate correlated noise.');
        end

        function names = classRequiredPropertyNames()
            names = {};
        end

        function distribution = distributionFromFile(path)
            % Initialize a persisted distribution from a NetCDF file.
            %
            % Use this abstract-role factory when the concrete distribution
            % subtype is only known from the persisted `AnnotatedClass`
            % metadata in the file.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: distribution = distributionFromFile(path)
            % - Parameter path: path to an existing NetCDF file
            % - Returns distribution: reconstructed concrete Distribution instance
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                distribution (1,1) Distribution
            end
            atc = CAAnnotatedClass.annotatedClassFromFile(path);
            if ~isa(atc, 'Distribution')
                error('Distribution:UnexpectedAnnotatedClass', 'The persisted AnnotatedClass %s does not inherit from Distribution.', class(atc));
            end
            distribution = atc;
        end

        function distribution = distributionFromGroup(group)
            % Initialize a persisted distribution from a NetCDF group.
            %
            % `distributionFromGroup` reconstructs the concrete subtype
            % stored in the group's `AnnotatedClass` metadata through the
            % shared annotated-class read path.
            %
            % - Topic: Inspect distribution properties
            % - Declaration: distribution = distributionFromGroup(group)
            % - Parameter group: NetCDF group containing a persisted distribution
            % - Returns distribution: reconstructed concrete Distribution instance
            arguments (Input)
                group (1,1) NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                distribution (1,1) Distribution
            end
            atc = CAAnnotatedClass.annotatedClassFromGroup(group);
            if ~isa(atc, 'Distribution')
                error('Distribution:UnexpectedAnnotatedClass', 'The persisted AnnotatedClass %s does not inherit from Distribution.', class(atc));
            end
            distribution = atc;
        end
    end

    methods (Static, Access = protected)
        function sz = sizeVectorFromInputs(varargin)
            if isempty(varargin)
                error('You must specify the size');
            end

            if numel(varargin) == 1
                sz = varargin{1};
                mustBeVector(sz);
                sz = reshape(sz,1,[]);
            else
                sz = zeros(1,numel(varargin));
                for i=1:numel(varargin)
                    validateattributes(varargin{i}, {'numeric'}, {'scalar'});
                    sz(i) = varargin{i};
                end
            end
        end
    end
    
    methods (Static, Access = private)
        function z = truncateToSupport(z,zrange)
            arguments
                z (1,1) {mustBeNumeric,mustBeReal}
                zrange (1,2) {mustBeNumeric,mustBeReal}
            end
            if z < zrange(1)
                z=zrange(1);
            end
            if z > zrange(2)
                z = zrange(2);
            end
        end
    end
end
