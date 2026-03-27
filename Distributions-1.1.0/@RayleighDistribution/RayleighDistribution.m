classdef RayleighDistribution < Distribution
    % Model a Rayleigh distribution with scale $$\sigma$$.
    %
    % `RayleighDistribution` implements the
    % [Rayleigh distribution](https://en.wikipedia.org/wiki/Rayleigh_distribution)
    % with density
    % $$p(z) = \frac{z}{\sigma^{2}}
    % \exp\!\left(-\frac{z^{2}}{2\sigma^{2}}\right), \quad z \ge 0,$$
    % cumulative distribution
    % $$F(z) = 1 - \exp\!\left(-\frac{z^{2}}{2\sigma^{2}}\right), \quad z \ge 0,$$
    % and variance
    % $$\mathrm{variance} = \frac{4-\pi}{2}\sigma^{2}.$$ This family is
    % the radial-distance law obtained from two independent zero-mean
    % Gaussian coordinates with common standard deviation $$\sigma.$$
    %
    % ```matlab
    % distribution = RayleighDistribution(1.0);
    % samples = distribution.rand([1000 1]);
    % ```
    %
    % - Topic: Create distributions
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Declaration: classdef RayleighDistribution < Distribution
    
    properties (SetAccess = private)
        % Scale parameter $$\sigma$$.
        %
        % `sigma` sets the radial scale of the distribution and fixes the
        % total variance at $$\frac{4-\pi}{2}\sigma^{2}.$$
        %
        % - Topic: Inspect distribution properties
        sigma
    end
    
    methods
        function self = RayleighDistribution(sigma)
            % Create a Rayleigh distribution from its scale parameter.
            %
            % Use this constructor for nonnegative radial data whose
            % density follows the Rayleigh family.
            %
            % - Topic: Create distributions
            % - Declaration: self = RayleighDistribution(sigma)
            % - Parameter sigma: positive scale parameter
            % - Returns self: RayleighDistribution instance
            arguments
                sigma (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive}
            end
            self.sigma = sigma;
            self.pdf = @(z) z.*exp(-(z.*z)/(2*sigma*sigma))/(sigma*sigma);
            self.cdf = @(z) 1 - exp(-(z.*z)/(2*sigma*sigma));
%             self.w = @(z) sigma*sigma*ones(size(z));
            self.variance = (4-pi)*sigma*sigma/2;
            self.logPDF = @(z) log(z)-(z.*z)/(2*sigma*sigma);
            
%             self.dPDFoverZ = @(z) -exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi))/(sigma*sigma);
        end
        
    end
end
