function [energyFrequency,omegaVector] = convertFromWavenumberToFrequency(self)
% Bin wave energy by vertical mode and intrinsic frequency
%
% Redistributes the energy in the stored Ap and Am coefficients from the
% horizontal-wavenumber grid onto uniformly spaced intrinsic-frequency bins.
%
% - Topic: Analyze flow
% - Declaration: [energyFrequency,omegaVector] = wvt.convertFromWavenumberToFrequency()
% - Returns energyFrequency: wave energy for each vertical mode and frequency bin, with dimensions `Nj`-by-`numel(omegaVector)`
% - Returns omegaVector: intrinsic-frequency bin coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%​
  % Defining omega vector based on biggest dOmega %​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=self.Omega;
omegaj1=omega(self.J==1);
dOmega=max(diff(sort(omegaj1(:))));
omegaVector=min(omega(:)):dOmega:max(omega(:));

energyFrequency=zeros(length(self.j),length(omegaVector));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ​
         % Redistributing the energy %         ​
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=(self.Am.*conj(self.Am) + self.Ap.*conj(self.Ap)).*self.Apm_TE_factor;
for iMode=1:self.Nj
    for iOmega=(1:length(omegaVector)-1)   
        % find all the kl point btw the two values of omega
        indForOmega = self.Omega >= omegaVector(iOmega) & self.Omega < omegaVector(iOmega+1) & self.J == iMode;
        energyFrequency(iMode,iOmega) = sum(E(indForOmega));
    end
end
end
