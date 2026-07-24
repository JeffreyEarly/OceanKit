atse = AlongTrackSimulatorEmpirical();
missionName = atse.missions{24};
missionName = "s6a"

ats = AlongTrackSimulator();
T_orbit = ats.orbitalPeriodForMissionWithName(missionName);


[lat_e,lon_e,time_e] = atse.groundTrackForMissionWithName(missionName);
t_e = seconds(time_e - time_e(1));

first_orbit = 120;
last_orbit = 127;
t = first_orbit*T_orbit:1:last_orbit*T_orbit;
N_orbits = 20;

indices = t_e<last_orbit*T_orbit & t_e > first_orbit*T_orbit;

%%
figure
% scatter(lon_e,lat_e,1.5^2,"blue","filled"), hold on
scatter(lon_e(indices),lat_e(indices),1.5^2,"blue","filled"), hold on
[lat,lon,time] = ats.groundTrackForMissionWithName(missionName,time=t);
scatter(lon,lat,1.5^2,"red","filled")