function [Mxx, Myy, Mxy] = momentTensorEvolution(self, Mxx0, Myy0, Mxy0, t, kappa)
% Evolve the second-moment tensor for the linear model.
%
% The returned moments satisfy the matrix equation
% $$\dot{M} = AM + MA^{\top} + 2\kappa I$$ for the affine velocity gradient
% `A` defined by the model parameters.
%
% - Topic: Analyze particle and moment evolution
% - Declaration: [Mxx, Myy, Mxy] = momentTensorEvolution(self,Mxx0,Myy0,Mxy0,t,kappa)
% - Parameter Mxx0: initial $$M_{xx}$$ entry in m^2
% - Parameter Myy0: initial $$M_{yy}$$ entry in m^2
% - Parameter Mxy0: initial $$M_{xy}$$ entry in m^2
% - Parameter t: time vector in seconds
% - Parameter kappa: scalar diffusivity in $$m^2 s^{-1}$$
% - Returns Mxx: time-dependent $$M_{xx}$$ entry in m^2 with the same shape as `t`
% - Returns Myy: time-dependent $$M_{yy}$$ entry in m^2 with the same shape as `t`
% - Returns Mxy: time-dependent $$M_{xy}$$ entry in m^2 with the same shape as `t`
arguments
    self (1,1) LinearVelocityField
    Mxx0 (1,1) double
    Myy0 (1,1) double
    Mxy0 (1,1) double
    t {mustBeNumeric}
    kappa (1,1) double
end

if self.sigma == 0 && self.zeta == 0
    Mxx = 2 * kappa * t + Mxx0;
    Myy = 2 * kappa * t + Myy0;
    Mxy = 0 * t + Mxy0;
elseif self.zeta == 0
    cos_t = cos(self.theta);
    sin_t = sin(self.theta);

    cos2 = cos_t * cos_t;
    sin2 = sin_t * sin_t;
    cossin = cos_t * sin_t;

    tks = 2 * kappa / self.sigma;

    Mxx0_r = cos2 * Mxx0 + sin2 * Myy0 + 2 * cos_t * sin_t * Mxy0;
    Myy0_r = sin2 * Mxx0 + cos2 * Myy0 - 2 * cos_t * sin_t * Mxy0;
    Mxy0_r = -cossin * Mxx0 + cossin * Myy0 + (cos2 - sin2) * Mxy0;

    Maa = (Mxx0_r + tks) * exp(self.sigma * t) - tks;
    Mbb = (Myy0_r - tks) * exp(-self.sigma * t) + tks;
    Mab = Mxy0_r;

    Mxx = cos2 * Maa + sin2 * Mbb - 2 * cos_t * sin_t * Mab;
    Myy = sin2 * Maa + cos2 * Mbb + 2 * cos_t * sin_t * Mab;
    Mxy = cossin * Maa - cossin * Mbb + (cos2 - sin2) * Mab;
elseif self.sigma == 0
    A = (Mxx0 + Myy0) / 2;
    B = -Mxy0;
    C = (Mxx0 - Myy0) / 2;

    Mxx = A + 2 * kappa * t + B * sin(self.zeta * t) + C * cos(self.zeta * t);
    Myy = A + 2 * kappa * t - B * sin(self.zeta * t) - C * cos(self.zeta * t);
    Mxy = -B * cos(self.zeta * t) + C * sin(self.zeta * t);
elseif self.zeta * self.zeta < self.sigma * self.sigma
    cos_t = cos(self.theta);
    sin_t = sin(self.theta);

    cos2 = cos_t * cos_t;
    sin2 = sin_t * sin_t;
    cossin = cos_t * sin_t;

    s = sqrt(self.sigma * self.sigma - self.zeta * self.zeta);
    tks = 2 * kappa * self.sigma / (s * s);

    Mxx1 = cos2 * Mxx0 + sin2 * Myy0 + 2 * cos_t * sin_t * Mxy0;
    Myy1 = sin2 * Mxx0 + cos2 * Myy0 - 2 * cos_t * sin_t * Mxy0;
    Mxy1 = -cossin * Mxx0 + cossin * Myy0 + (cos2 - sin2) * Mxy0;

    A = (1 + self.sigma / s) * Mxx1 / 2 - (self.zeta / s) * Mxy1 - (1 - self.sigma / s) * Myy1 / 2 + tks;
    B = (1 - self.sigma / s) * Mxx1 / 2 + (self.zeta / s) * Mxy1 - (1 + self.sigma / s) * Myy1 / 2 + tks;
    C = -(self.zeta / s) * Mxx1 + (2 * self.sigma / s) * Mxy1 - (self.zeta / s) * Myy1;

    Maa = (A / 2) * (1 + self.sigma / s) * exp(s * t) + (B / 2) * (1 - self.sigma / s) * exp(-s * t) + (self.zeta / s) * C / 2 - (2 * kappa / (s * s)) * (self.zeta * self.zeta * t + self.sigma);
    Mbb = -(A / 2) * (1 - self.sigma / s) * exp(s * t) - (B / 2) * (1 + self.sigma / s) * exp(-s * t) + (self.zeta / s) * C / 2 - (2 * kappa / (s * s)) * (self.zeta * self.zeta * t - self.sigma);
    Mab = (A / 2) * (self.zeta / s) * exp(s * t) - (B / 2) * (self.zeta / s) * exp(-s * t) + (self.sigma / s) * C / 2 - (2 * kappa / (s * s)) * (self.zeta * self.sigma * t);

    Mxx = cos2 * Maa + sin2 * Mbb - 2 * cos_t * sin_t * Mab;
    Myy = sin2 * Maa + cos2 * Mbb + 2 * cos_t * sin_t * Mab;
    Mxy = cossin * Maa - cossin * Mbb + (cos2 - sin2) * Mab;
elseif self.zeta * self.zeta == self.sigma * self.sigma
    cos_t = cos(self.theta);
    sin_t = sin(self.theta);

    cos2 = cos_t * cos_t;
    sin2 = sin_t * sin_t;
    cossin = cos_t * sin_t;

    Mxx1 = cos2 * Mxx0 + sin2 * Myy0 + 2 * cos_t * sin_t * Mxy0;
    Myy1 = sin2 * Mxx0 + cos2 * Myy0 - 2 * cos_t * sin_t * Mxy0;
    Mxy1 = -cossin * Mxx0 + cossin * Myy0 + (cos2 - sin2) * Mxy0;

    A = -2 * self.zeta * Mxy1 + self.sigma * (Mxx1 + Myy1);
    B = Mxx1 - Myy1;
    C = Mxx1 + Myy1;

    Maa = kappa * self.sigma^2 * t.^3 / 3 + kappa * self.sigma * t.^2 + 2 * kappa * t + A * (self.sigma * t.^2 + 2 * t) / 4 + B * (self.sigma * t + 1) / 2 + C / 2;
    Mbb = kappa * self.sigma^2 * t.^3 / 3 - kappa * self.sigma * t.^2 + 2 * kappa * t + A * (self.sigma * t.^2 - 2 * t) / 4 + B * (self.sigma * t - 1) / 2 + C / 2;
    Mab = kappa * self.sigma * self.zeta * t.^3 / 3 + A * (self.zeta * t.^2 / 4 - 1 / (2 * self.zeta)) + B * self.zeta * t / 2 + self.sigma * C / (2 * self.zeta);

    Mxx = cos2 * Maa + sin2 * Mbb - 2 * cos_t * sin_t * Mab;
    Myy = sin2 * Maa + cos2 * Mbb + 2 * cos_t * sin_t * Mab;
    Mxy = cossin * Maa - cossin * Mbb + (cos2 - sin2) * Mab;
else
    cos_t = cos(self.theta);
    sin_t = sin(self.theta);

    cos2 = cos_t * cos_t;
    sin2 = sin_t * sin_t;
    cossin = cos_t * sin_t;

    s = sqrt(self.zeta * self.zeta - self.sigma * self.sigma);

    Mxx1 = cos2 * Mxx0 + sin2 * Myy0 + 2 * cos_t * sin_t * Mxy0;
    Myy1 = sin2 * Mxx0 + cos2 * Myy0 - 2 * cos_t * sin_t * Mxy0;
    Mxy1 = -cossin * Mxx0 + cossin * Myy0 + (cos2 - sin2) * Mxy0;

    A = Mxx1 - Myy1 - 4 * kappa * self.sigma / s^2;
    B = (self.sigma / s) * Mxx1 - (2 * self.zeta / s) * Mxy1 + (self.sigma / s) * Myy1;
    C = (self.zeta^2 / s^2) * Mxx1 - (2 * self.zeta * self.sigma / s^2) * Mxy1 + (self.zeta^2 / s^2) * Myy1;

    Maa = (A / 2) * (cos(s * t) + (self.sigma / s) * sin(s * t)) + (B / 2) * (sin(s * t) - (self.sigma / s) * cos(s * t)) + C / 2 + (2 * kappa / s^2) * (self.sigma + self.zeta^2 * t);
    Mbb = -(B / 2) * (sin(s * t) + (self.sigma / s) * cos(s * t)) - (A / 2) * (cos(s * t) - (self.sigma / s) * sin(s * t)) + C / 2 + (2 * kappa / s^2) * (-self.sigma + self.zeta^2 * t);
    Mab = (A * self.zeta / (2 * s)) * sin(s * t) - (B * self.zeta / (2 * s)) * cos(s * t) + C * self.sigma / (2 * self.zeta) + 2 * kappa * self.sigma * self.zeta * t / s^2;

    Mxx = cos2 * Maa + sin2 * Mbb - 2 * cos_t * sin_t * Mab;
    Myy = sin2 * Maa + cos2 * Mbb + 2 * cos_t * sin_t * Mab;
    Mxy = cossin * Maa - cossin * Mbb + (cos2 - sin2) * Mab;
end
