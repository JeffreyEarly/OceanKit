ats = AlongTrackSimulator();
missionName = ats.missions{5}
% missionName = "s6a"

[lat,lon,time] = ats.repeatGroundTrackForMissionWithName(missionName);

lag = 3000;

figure
scatter(lon(1:lag),lat(1:lag),1.5^2,"red","filled"), hold on
scatter(lon((end-lag):end),lat((end-lag):end),1.5^2,"blue","filled")
legend('cycle start','cycle end')

% %%
% figure
% scatter(lon,lat,1.5^2,"red","filled"), hold on