classdef Chi2Distribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k
    end
    
    methods
        function self = Chi2Distribution(k)
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

