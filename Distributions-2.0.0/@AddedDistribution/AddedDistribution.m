classdef AddedDistribution < Distribution
    % Model a weighted mixture of component distributions.
    %
    % `AddedDistribution` implements a
    % [mixture distribution](https://en.wikipedia.org/wiki/Mixture_distribution)
    % by combining component densities $$p_i(z)$$ with nonnegative weights
    % $$\alpha_i$$ that sum to one. The resulting density and cumulative
    % distribution are
    % $$p(z) = \sum_{i=1}^{N} \alpha_i p_i(z),$$
    % $$F(z) = \sum_{i=1}^{N} \alpha_i F_i(z),$$
    % and the stored total variance is the weighted sum
    % $$\mathrm{variance} = \sum_{i=1}^{N} \alpha_i \,
    % \mathrm{variance}_i.$$ When `scalings` is supplied as a scalar, the
    % constructor interprets it as the first weight in a two-component
    % mixture and expands it to $$[\alpha, 1-\alpha].$$
    %
    % ```matlab
    % distribution = AddedDistribution(scalings=[0.2 0.8], distributions=[NormalDistribution(sigma=2.0); StudentTDistribution(nu=4.5, sigma=1.0)]);
    % ```
    %
    % - Topic: Create distributions
    % - Topic: Inspect distribution properties
    % - Topic: Sample from distributions
    % - Topic: Evaluate distribution fit
    % - Topic: Model correlated noise
    % - Topic: Compose distributions
    % - Declaration: classdef AddedDistribution < Distribution
    
    properties (SetAccess = private)
        % Mixture weights $$\alpha_i$$ for the component distributions.
        %
        % `scalings` stores the weights applied to each component in the
        % order they are passed to the constructor.
        %
        % - Topic: Compose distributions
        scalings

        % Component distributions used in the mixture model.
        %
        % `distributions` stores the component distributions whose
        % densities and CDFs are combined by the mixture weights.
        %
        % - Topic: Compose distributions
        distributions
    end

    properties (Dependent, Hidden, SetAccess = private)
        distributionIndex
    end
    
    methods
        function self = AddedDistribution(options)
            % Create a weighted mixture from two or more distributions.
            %
            % Supply either one scalar weight for a two-component mixture
            % or one weight per component. The `distributions` input lists
            % the component distributions in mixture order.
            %
            % - Topic: Compose distributions
            % - Declaration: self = AddedDistribution(scalings=<value>,distributions=<value>)
            % - Parameter options.scalings: mixture weights that sum to one
            % - Parameter options.distributions: component distributions in mixture order
            % - Returns self: AddedDistribution instance
            arguments
                options.scalings {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
                options.distributions {mustBeA(options.distributions,"Distribution"),mustBeVector}
            end
            scalings = reshape(options.scalings, 1, []);
            distributionArray = reshape(options.distributions, [], 1);
            if length(distributionArray) < 2
                error('AddedDistribution:NeedAtLeastTwoComponents', 'You must supply at least two component distributions.');
            end
            if length(scalings) == 1
                if length(distributionArray) ~= 2
                    error('AddedDistribution:ScalarScalingRequiresTwoComponents', 'A scalar scaling is only valid for a two-component mixture.');
                end
                scalings(2) = 1 - scalings;
            end
            if any(scalings < 0) || sum(scalings) ~= 1.0
                error('AddedDistribution:InvalidScalings', 'The scalings must sum to 1.0, and all values must be between 0 and 1.');
            end
            if length(distributionArray) ~= length(scalings)
                error('AddedDistribution:ComponentCountMismatch', 'There must be a scaling value for all distributions.');
            end
            self.scalings = scalings;
            self.distributions = distributionArray;
            
            self.pdf = @(z) self.summedPDF(z);
            self.cdf = @(z) self.summedCDF(z);
            self.dPDFOverZ = @(z) self.summedDPDFOverZ(z);
            self.w = @(z) self.summedWeight(z);
            self.logPDF = @(z) log(self.pdf(z));
            self.logPDFNorm = self.summedPDFNorm();

            self.variance = 0;
            for i = 1:length(self.distributions)
                self.variance = self.variance + self.scalings(i) * self.distributions(i).variance;
            end
        end

        function distributionIndex = get.distributionIndex(self)
            distributionIndex = reshape(1:length(self.distributions), [], 1);
        end
    end

    methods (Static, Hidden)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = Distribution.classDefinedPropertyAnnotations();
            propertyAnnotations(end+1) = CADimensionProperty('distributionIndex', '', 'Mixture component index.');
            propertyAnnotations(end+1) = CANumericProperty('scalings', {'distributionIndex'}, '', 'Mixture weights $$\\alpha_i$$.');
            propertyAnnotations(end+1) = CAObjectProperty('distributions', 'Array of component Distribution objects.');
        end

        function names = classRequiredPropertyNames()
            names = {'scalings', 'distributions'};
        end
    end

    methods (Access = private)
        function pdf = summedPDF(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            pdf = zeros(size(z));
            for i = 1:length(self.distributions)
               pdf = pdf + self.scalings(i) * self.distributions(i).pdf(z);
            end
        end

        function pdfNorm = summedPDFNorm(self)
            arguments
                self (1,1) AddedDistribution
            end
            pdfNorm = 0;
            for i = 1:length(self.distributions)
               pdfNorm = pdfNorm + self.scalings(i) * self.distributions(i).logPDFNorm;
            end
        end
                
        function cdf = summedCDF(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            cdf = zeros(size(z));
            for i = 1:length(self.distributions)
                cdf = cdf + self.scalings(i) * self.distributions(i).cdf(z);
            end
        end
        
        function dPDFOverZ = summedDPDFOverZ(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            dPDFOverZ = zeros(size(z));
            for i = 1:length(self.distributions)
                dPDFOverZ = dPDFOverZ + self.scalings(i) * self.distributions(i).dPDFOverZ(z);
            end
        end
        
        function w = summedWeight(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            w = -self.pdf(z) ./ self.dPDFOverZ(z);
        end
    end
end
