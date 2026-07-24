% (xe,ye) are both function handles, that take time as a parameter. So
% although the shape of the eddy is limited to a two parameter gaussian,
% it's propagation path is entirely unspecified.
test_eddy = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);

% Let's define a path analytically... but we could use an existing path
% using interp
x0 = 1500e3; y0 = 1500e3; vx = -2.0e-2; vy = -0.3e-2;
x_lin = @(t) x0+vx*t;
y_lin = @(t) y0+vy*t;

% Lets actually make an eddy with a chosen set of parameters.
my_eddy = @(x,y,t) test_eddy(x,y,t,0.15,80e3,x_lin,y_lin);

ats = AlongTrackSimulator();
currentMissions = ats.currentMissions;
ats.summarizeMissionWithName(currentMissions);

time = (0:10*86400).';

figure
tl = tiledlayout("flow",TileSpacing="tight");
for iMission = 1:length(ats.currentMissions)
    obs = ats.projectedPointsForMissionWithName(currentMissions(iMission),Lx=2000e3,Ly=2000e3,lat0=30,lon0=0,time=time);
    obs.ssh = my_eddy(obs.x,obs.y,obs.t);
    nexttile
    scatter3(obs.x/1e3,obs.y/1e3,obs.ssh*1e2,[],obs.ssh*1e2,'filled')
    repeatCycle = ats.repeatCycleForMissionWithName(currentMissions(iMission))/86400;
    longname = ats.missionParameters(currentMissions(iMission)).name;
    if isinf(repeatCycle)
        titleString = longname + " (" + currentMissions(iMission) + ") geodetic orbit";
    else
        titleString = longname + " (" + currentMissions(iMission) + ") " + string(repeatCycle) + " day repeat cycle";
    end
    title(titleString)
end
title(tl,"Ground tracks for " + string(max(time)/86400) + " days")

% % Now fit this to a simple model eddy...
% model_eddy = @(x,y,t,A,L,x0,y0,cx,cy) A.*exp(-((x-x0-cx*t).^2 + (y-y0-cy*t).^2)/L^2);
% penalty_function = @(p) mean((obs.ssh - model_eddy(obs.x,obs.y,obs.t,p(1),p(2),p(3),p(4),p(5),p(6))).^2);
% 
% % feed in slightly worse values. These values are terribly scaled, but it
% % seems to return exactly the right value anyway.
% pmin=fminsearch(penalty_function,[0.13, 85e3, 1450e3, 1450e3, -2e-2, -3e-3]);