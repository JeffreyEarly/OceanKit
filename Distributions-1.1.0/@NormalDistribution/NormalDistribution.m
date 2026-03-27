classdef NormalDistribution < Distribution
    % Model a zero-mean normal distribution with standard deviation $$\sigma$$.
    %
    % `NormalDistribution` implements the
    % [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution)
    % with density
    % $$p(z) = \frac{1}{\sigma \sqrt{2\pi}}
    % \exp\!\left(-\frac{z^{2}}{2\sigma^{2}}\right),$$
    % cumulative distribution
    % $$F(z) = \frac{1}{2}\left(1 + \operatorname{erf}\!\left(\frac{z}{\sigma \sqrt{2}}\right)\right),$$
    % and variance $$\sigma^{2}.$$
    %
    % ```matlab
    % distribution = NormalDistribution(0.5);
    % samples = distribution.rand([1000 1]);
    % ```
    %
    % - Topic: Create distributions
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Declaration: classdef NormalDistribution < Distribution
    
    properties (SetAccess = private)
        % Standard deviation $$\sigma$$.
        %
        % `sigma` sets the scale of the Gaussian density and fixes the
        % total variance at $$\sigma^{2}.$$
        %
        % - Topic: Inspect distribution properties
        sigma
    end
    
    methods
        function self = NormalDistribution(sigma)
            % Create a normal distribution from its standard deviation.
            %
            % Use this constructor for a zero-mean Gaussian model with
            % variance $$\sigma^{2}.$$
            %
            % - Topic: Create distributions
            % - Declaration: self = NormalDistribution(sigma)
            % - Parameter sigma: positive standard deviation
            % - Returns self: NormalDistribution instance
            arguments
                sigma (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive}
            end
            self.sigma = sigma;
            self.pdf = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
            self.cdf = @(z) 0.5*(1 + erf(z/(sigma*sqrt(2))));
            self.w = @(z) sigma*sigma*ones(size(z));
            self.variance = sigma*sigma;
            
            self.dPDFOverZ = @(z) -exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi))/(sigma*sigma);
            self.logPDF = @(z) -(z.*z)/(2*sigma*sigma)-(log(sigma)+log(2*pi)/2);
            self.logPDFNorm = -(log(sigma)+log(2*pi)/2);
        end
        
        function y = rand(self,sz)
            % Draw normally distributed random samples.
            %
            % This override uses MATLAB's `randn` directly and expects one
            % size vector `sz` describing the requested output shape.
            %
            % - Topic: Sample from distributions
            % - Declaration: y = rand(self,sz)
            % - Parameter sz: size vector of nonnegative integers
            % - Returns y: array of Gaussian random samples
            arguments
                self (1,1) NormalDistribution
                sz {mustBeNumeric,mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative,mustBeVector}
            end
            sz = reshape(sz,1,[]);
            y = self.sigma*randn(sz);
        end

    end
end
