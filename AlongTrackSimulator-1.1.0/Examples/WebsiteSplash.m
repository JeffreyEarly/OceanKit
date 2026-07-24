ats = AlongTrackSimulator();

fig = figure(Name="Ground tracks",Position=[100 100 500 500]);
[lat,lon,~] = ats.groundTrackForMissionWithName("s6a",N_orbits=5);
geoscatter(lat,lon,1.5^2,"blue","filled"), hold on
[lat,lon,~] = ats.groundTrackForMissionWithName("swon",N_orbits=5);
geoscatter(lat,lon,1.5^2,"red","filled")
legend(ats.missionParameters('s6a').name,ats.missionParameters('swon').name,Location="southeast")

geobasemap landcover

exportgraphics(fig,"groundtracks_splash.jpg",Resolution=250)