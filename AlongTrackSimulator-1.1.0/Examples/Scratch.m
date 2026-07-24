path = "./orbital_paths.nc";
info = ncinfo("@AlongTrackSimulator/" + path);
missions = {info.Groups.Name};

d = dictionary();
for iMission = 1:length(missions)
    mission = missions{iMission};
    var.startTime = days(ncread("@AlongTrackSimulator/" + path, mission + "/start_time")) + datetime(1950,01,01);
    var.endTime = ncread("@AlongTrackSimulator/" + path, mission + "/end_time") + datetime(1950,01,01);
    var.lon = ncread("@AlongTrackSimulator/" + path, mission + "/longitude");
    var.lat = ncread("@AlongTrackSimulator/" + path, mission + "/latitude");
    var.time = days(ncread("@AlongTrackSimulator/" + path, mission + "/time")) + datetime(1950,01,01);
    d(mission) = var;
end