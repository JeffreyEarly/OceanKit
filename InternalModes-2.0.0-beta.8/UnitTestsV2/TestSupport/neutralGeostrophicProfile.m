function profile = neutralGeostrophicProfile()
% Return the prescribed neutral-cap profile and its physical joins.
profile.zDomain = [-1000 0];
profile.joins = [-1000 -400 -200 0];
profile.N2 = @evaluateN2;
profile.integratedN2 = 0.0175;
profile.description = "N2=0 above -200 m; quintic transition on [-400,-200] m; N=0.005 s^-1 below -400 m";
end

function values = evaluateN2(z)
t = min(1,max(0,(-200-z)/200));
values = (5e-3)^2*(10*t.^3-15*t.^4+6*t.^5);
end
