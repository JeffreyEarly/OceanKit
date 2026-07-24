ats = AlongTrackSimulator();

figure
alongtrack = ats.projectedPointsForRepeatMissionWithName("s6a",Lx=2000e3,Ly=2000e3,lat0=30,lon0=1);
scatter(alongtrack.x, alongtrack.y,1.5^2,"blue","filled"), hold on
alongtrack = ats.projectedPointsForRepeatMissionWithName("j3n",Lx=2000e3,Ly=2000e3,lat0=30,lon0=1);
scatter(alongtrack.x, alongtrack.y,1.5^2,"red","filled")
%%
% figure
% scatter(x, y,1.5^2,"red","filled")

% 
% lag = 3000;
% 
% figure
% scatter(lon(1:lag),lat(1:lag),1.5^2,"red","filled"), hold on
% scatter(lon((end-lag):end),lat((end-lag):end),1.5^2,"blue","filled")
% legend('cycle start','cycle end')