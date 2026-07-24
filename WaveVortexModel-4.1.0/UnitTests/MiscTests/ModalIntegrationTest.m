Lxy = 750e3;
N = 16;
Nz = 60;

Lz = 4000;
N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
L_gm = 1300; % thermocline exponential scale, meters
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([Lxy, Lxy, Lz], [N, N, Nz], N2=N2,latitude=31);

% coeff = (wvt.Q ./ wvt.P);

%%
v = wvt.PF0inv(:,2);

IntZF = wvt.QG0inv*(squeeze(wvt.h_0 .* wvt.Q0 ./ wvt.P0).*wvt.PF0);
intv = IntZF * v;
intv = wvt.intZF(v);

figure
tiledlayout(1,2)
nexttile
plot(v,wvt.z)
nexttile
plot(intv,wvt.z), hold on
plot(cumtrapz(wvt.z,v),wvt.z)

%%

v = wvt.QG0inv(:,2);

PF0inv = wvt.PF0inv;
PF0inv(:,2:end) = PF0inv(:,2:end) - PF0inv(1,2:end);
% PF0inv = wvt.PF0inv - wvt.PF0inv(1,:);
IntZG = - wvt.g*PF0inv*(squeeze(wvt.P0./wvt.Q0).*(wvt.QG0 .* shiftdim(1./wvt.N2,-1)));
intv = IntZG * v;
intv = wvt.intZG(v);

figure
tiledlayout(1,2)
nexttile
plot(v,wvt.z)
nexttile
plot(intv,wvt.z), hold on
plot(cumtrapz(wvt.z,v),wvt.z)
