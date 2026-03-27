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
    % distribution = AddedDistribution([0.2 0.8], ...
    %     NormalDistribution(2.0), ...
    %     StudentTDistribution(nu=4.5, sigma=1.0));
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
    
    methods
        function self = AddedDistribution(scalings,distribution1,distribution2,varargin)
            % Create a weighted mixture from two or more distributions.
            %
            % Supply either one scalar weight for a two-component mixture
            % or one weight per component. The remaining positional inputs
            % are the component distributions in mixture order.
            %
            % - Topic: Compose distributions
            % - Declaration: self = AddedDistribution(scalings,distribution1,distribution2,varargin)
            % - Parameter scalings: mixture weights that sum to one
            % - Parameter distribution1: first component distribution
            % - Parameter distribution2: second component distribution
            % - Parameter varargin: additional component distributions
            % - Returns self: AddedDistribution instance
            arguments
                scalings {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
                distribution1 (1,1) Distribution
                distribution2 (1,1) Distribution
            end
            arguments (Repeating)
                varargin (1,1) Distribution
            end
            if length(scalings) == 1
                scalings(2) = 1-scalings;
            end
            if any(scalings<0) || sum(scalings) ~= 1.0
                error('The scalings must sum to 1.0, and all values must be between 0 and 1.');
            end
            self.scalings = scalings;
            
            nDistributions = 2 + length(varargin);
            if nDistributions ~= length(scalings)
                error('There must be a scaling value for all distrubtions.');
            end
            
            self.distributions = cell(nDistributions,1);
            self.distributions{1} = distribution1;
            self.distributions{2} = distribution2;
            for i = 1:length(varargin)
                self.distributions{i+2} = varargin{i};
            end
            
            self.pdf = @(z) self.summedPDF(z);
            self.cdf = @(z) self.summedCDF(z);
            self.dPDFOverZ = @(z) self.summedDPDFOverZ(z);
            self.w = @(z) self.summedWeight(z);
            self.logPDF = @(z) log(self.pdf(z));
            self.logPDFNorm = self.summedPDFNorm();

            self.variance = 0;
            for i=1:length(self.distributions)
                self.variance = self.variance + self.scalings(i)*self.distributions{i}.variance;
            end
            
        end
        
    end

    methods (Access = private)
        function pdf = summedPDF(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            pdf = zeros(size(z));
            for i=1:length(self.distributions)
               pdf = pdf + self.scalings(i)*self.distributions{i}.pdf(z); 
            end
        end

        function pdfNorm = summedPDFNorm(self)
            arguments
                self (1,1) AddedDistribution
            end
            pdfNorm = 0;
            for i=1:length(self.distributions)
               pdfNorm = pdfNorm + self.scalings(i)*self.distributions{i}.logPDFNorm;
            end
        end
                
        function cdf = summedCDF(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            cdf = zeros(size(z));
            for i=1:length(self.distributions)
                cdf = cdf + self.scalings(i)*self.distributions{i}.cdf(z);
            end
        end
        
        function dPDFOverZ = summedDPDFOverZ(self,z)
            arguments
                self (1,1) AddedDistribution
                z {mustBeNumeric,mustBeReal}
            end
            dPDFOverZ = zeros(size(z));
            for i=1:length(self.distributions)
                dPDFOverZ = dPDFOverZ + self.scalings(i)*self.distributions{i}.dPDFOverZ(z);
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
