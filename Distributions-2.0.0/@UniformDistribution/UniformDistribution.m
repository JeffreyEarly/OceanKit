classdef UniformDistribution < Distribution
    % Model a continuous uniform distribution on the interval $$[a,b]$$.
    %
    % `UniformDistribution` implements the
    % [continuous uniform distribution](https://en.wikipedia.org/wiki/Continuous_uniform_distribution)
    % with density
    % $$p(z) = \frac{1}{b-a}, \quad a \le z \le b,$$
    % cumulative distribution
    % $$F(z) = \frac{z-a}{b-a}, \quad a \le z \le b,$$
    % together with the usual piecewise extensions to $$0$$ below $$a$$
    % and $$1$$ above $$b$$. The total variance is
    % $$\mathrm{variance} = \frac{(b-a)^{2}}{12}.$$
    %
    % ```matlab
    % distribution = UniformDistribution(a=-0.5, b=0.5);
    % samples = distribution.rand([1000 1]);
    % ```
    %
    % - Topic: Create distributions
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Declaration: classdef UniformDistribution < Distribution
    
    properties (SetAccess = private)
        % Lower endpoint $$a$$ of the support interval.
        %
        % `a` stores the inclusive lower bound of the uniform support.
        %
        % - Topic: Inspect distribution properties
        a

        % Upper endpoint $$b$$ of the support interval.
        %
        % `b` stores the inclusive upper bound of the uniform support.
        %
        % - Topic: Inspect distribution properties
        b
    end
    
    methods
        function self = UniformDistribution(options)
            % Create a continuous uniform distribution on $$[a,b]$$.
            %
            % When called with no inputs, the constructor uses the default
            % interval $$[-0.5,0.5].$$ If either endpoint is omitted, the
            % remaining endpoint uses its default value.
            %
            % - Topic: Create distributions
            % - Declaration: self = UniformDistribution(a=<value>,b=<value>)
            % - Parameter options.a: lower endpoint of the interval, default `-0.5`
            % - Parameter options.b: upper endpoint of the interval, default `0.5`
            % - Returns self: UniformDistribution instance
            arguments
                options.a (1,1) {mustBeNumeric,mustBeReal,mustBeFinite} = -0.5
                options.b (1,1) {mustBeNumeric,mustBeReal,mustBeFinite} = 0.5
            end
            a = options.a;
            b = options.b;
            self.a = a;
            self.b = b;
            self.pdf = @(z) (z < self.a | z > self.b)*0 + (z >= self.a & z <= self.b)*(1/(self.b-self.a)) ;
            self.cdf = @(z) (z < self.a)*0 + (z >= self.a & z <= self.b)*((z-self.a)/(self.b-self.a)) +  (z > self.b)*1;
%             self.w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
            self.variance = (1/12)*(self.b-self.a)^2;
            
%             self.dPDFoverZ = @(z) -2*(a*m/c)*(1+z.*z/c).^(-m-1);
%             self.logPDF = @(z) -m*log(1+z.*z/c);
        end
        
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = Distribution.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CANumericProperty('a', {}, '', 'Lower endpoint $$a$$ of the support interval.');
            propertyAnnotations(end+1) = CANumericProperty('b', {}, '', 'Upper endpoint $$b$$ of the support interval.');
        end

        function names = classRequiredPropertyNames()
            names = {'a', 'b'};
        end
    end

end
