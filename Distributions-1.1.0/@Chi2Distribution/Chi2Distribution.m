classdef Chi2Distribution < Distribution
    % Model a chi-squared distribution with $$k$$ degrees of freedom.
    %
    % `Chi2Distribution` implements the
    % [chi-squared distribution](https://en.wikipedia.org/wiki/Chi-squared_distribution)
    % with density
    % $$p(z) = \frac{z^{k/2-1} e^{-z/2}}
    % {2^{k/2}\Gamma(k/2)}, \quad z \ge 0,$$
    % cumulative distribution
    % $$F(z) = P\!\left(\frac{k}{2},\frac{z}{2}\right),$$
    % where $$P$$ is the regularized lower incomplete gamma function, and
    % variance $$\mathrm{variance} = 2k.$$
    %
    % ```matlab
    % distribution = Chi2Distribution(4);
    % samples = distribution.rand([1000 1]);
    % ```
    %
    % - Topic: Create distributions
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Declaration: classdef Chi2Distribution < Distribution
    
    properties (SetAccess = private)
        % Degrees of freedom $$k$$.
        %
        % `k` determines the shape of the chi-squared density and fixes
        % the total variance at $$2k.$$
        %
        % - Topic: Inspect distribution properties
        k
    end
    
    methods
        function self = Chi2Distribution(k)
            % Create a chi-squared distribution from its degrees of freedom.
            %
            % Use this constructor for nonnegative data modeled by a
            % chi-squared law with `k` degrees of freedom.
            %
            % - Topic: Create distributions
            % - Declaration: self = Chi2Distribution(k)
            % - Parameter k: positive degrees of freedom
            % - Returns self: Chi2Distribution instance
            arguments
                k (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive}
            end
            self.k = k;
            self.pdf = @(z) (z.^(k/2-1)).*exp(-z/2)/( (2^(k/2))*gamma(k/2) );
            self.cdf = @(z) gammainc(z/2,k/2);

            % self.w = @(z) sigma*sigma*ones(size(z));
            self.variance = 2*k;
            
            % self.dPDFoverZ = @(z) -exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi))/(sigma*sigma);
            % self.logPDF = @(z) -(z.*z)/(2*sigma*sigma)-(log(sigma)+log(2*pi)/2);
            % self.logPDFNorm = -(log(sigma)+log(2*pi)/2);
        end
        
        % function y = rand(self,sz)
        %    y = self.sigma*randn(sz); 
        % end

    end
end
