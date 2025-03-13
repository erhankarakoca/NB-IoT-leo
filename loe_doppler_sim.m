% % Frekans: 2 GHz
% % 
% % Uydu Yörünge Verisi:
% % 
% % STARLINK-1728
% % 1 46582U 20070BC  24252.77330876  .00021462  00000-0  14554-2 0  9996
% % 2 46582  53.0530  91.8195 0001518  88.8264 271.2898 15.06398015217350
% % 
% % Gözlem Aralığı
% % 2024-09-09T07:27:00Z - 2024-09-09T07:38:00Z (GMT+0)
% % 
% % Referans Noktasının Konumu: (40.789549, 29.451292) - (Lat, Lon)
% % 
% % 4 farklı UE için konumlar (50 km beam yarıçapı için):
% % Kuzeydoğu (NE): (41.1081, 29.8720)
% % Kuzeybatı (NW): (41.1081, 29.0306)
% % Güneydoğu (SE): (40.4710, 29.8720)
% % Güneybatı (SW): (40.4710, 29.0306)

%%

% hs = 1500e3; % meters
% el = 45;     % degrees
% freq = 5e9;  % Hz
% hg = 0;      % meters
% 
% time = 0:6950; % seconds
% 
% shift = dopplerShiftCircularOrbit(el,hs,hg,freq,time);
% 
% figure
% plot(time,shift)
% title("Doppler Shift vs Time")
% xlabel("Time (seconds)")
% ylabel("Doppler Shift (Hz)")
% grid on

%%
% % Gözlem Aralığı
% % 2024-09-09T07:27:00Z - 2024-09-09T07:38:00Z (GMT+0)
% startTime = datetime(2024,9,9,7,27,0);
startTime = datetime(2024,9,9,7,28,32);
stopTime = startTime + seconds(500);
sampleTime = 1;                         % In seconds
sc = satelliteScenario( ...
    startTime,stopTime,sampleTime);

tleFile = "starlink-1728.tle";
sat = satellite(sc, tleFile);


% Referans Noktasının Konumu: (40.789549, 29.451292) - (Lat, Lon)
% % Kuzeydoğu (NE): (41.1081, 29.8720)
% % Kuzeybatı (NW): (41.1081, 29.0306)
% % Güneydoğu (SE): (40.4710, 29.8720)
% % Güneybatı (SW): (40.4710, 29.0306)
ref_ground = [40.789549, 29.451292]; 
ref_name = "REF";

NE_ground = [41.1081, 29.8720];
NE_name = "NE";

NW_ground = [41.1081, 29.0306];
NW_name = "NW";

SE_ground = [40.4710, 29.8720];
SE_name = "SE";

SW_ground = [40.4710, 29.0306];
SW_name = "SW";

reference_point = groundStation(sc, "Name", ref_name,...
                          "Latitude",  ref_ground(1), ...
                          "Longitude", ref_ground(2));
reference_point.MarkerColor = [0 0.4470 0.7410]	;


gs_NE = groundStation(sc, 'Name', NE_name,...
                          'Latitude',  NE_ground(1), ...
                          'Longitude', NE_ground(2));
gs_NE.MarkerColor = [0.8500 0.3250 0.0980];
gs_NE.ShowLabel = false;

gs_NW = groundStation(sc, "Name",NW_name,...
                          "Latitude",  NW_ground(1), ...
                          "Longitude", NW_ground(2));
gs_NW.MarkerColor = [0.9290 0.6940 0.1250];
gs_NW.ShowLabel = false;


gs_SE = groundStation(sc, "Name", SE_name,...
                          "Latitude",  SE_ground(1), ...
                          "Longitude", SE_ground(2));
gs_SE.MarkerColor = [0.4940 0.1840 0.5560];
gs_SE.ShowLabel = false;

gs_SW = groundStation(sc, "Name", SW_name,...
                          "Latitude",  SW_ground(1), ...
                          "Longitude", SW_ground(2));
gs_SW.MarkerColor = [0.4660 0.6740 0.1880];
gs_SW.ShowLabel = false;

% ac_reference = access(sat,reference_point);

% ac_NE = access(sat,gs_NE);

% ac_NW = access(sat,gs_NW);

% ac_SE = access(sat,gs_SE);

% ac_SW = access(sat,gs_SW);

play(sc);

%%
% orbitTime = ;
% 09-09-2024-07:31:34;
% 09-09-2024-07:32:37;
% 09-09-2024-07:33:41;
orbitTime_1 = datetime(2024,9,9,07,31,34);
orbitTime_2 = datetime(2024,9,9,07,32,37);
orbitTime_3 = datetime(2024,9,9,07,33,41);

tle = tleread("starlink-1728.tle");
[satPos,satVelocity] = propagateOrbit(orbitTime_1, ...
                                 tle, ...
                                 "OutputCoordinateFrame","geographic");

sat_altitude = satPos(3);
[az, el, ~] = aer(gs_NE, sat, orbitTime_1);
disp("Satellite altitude  : " + sat_altitude);
disp("Satellite velocity  : " + norm(satVelocity));
disp("Satellite elevation : " + el);

%%

[delay_reference, time_reference] = latency(sat,reference_point);
[delay_NE,time_NE] = latency(sat,gs_NE);
[delay_NW,time_NW] = latency(sat,gs_NW);
[delay_SE,time_SE] = latency(sat,gs_SE);
[delay_SW,time_SW] = latency(sat,gs_SW);

figure
plot(time_reference,delay_reference(1,:)*1000)                  % Plot in milliseconds
xlim([time_reference(1) time_reference(end)])
title("Latency")
xlabel("Simulation Time")
ylabel("Latency (ms)")
grid on


hold on 
plot(time_NE,delay_NE(1,:)*1000)                  % Plot in milliseconds
hold on 

plot(time_NW,delay_NW(1,:)*1000)                  % Plot in milliseconds
hold on 

plot(time_SE,delay_SE(1,:)*1000)                  % Plot in milliseconds
hold on 

plot(time_SW,delay_SW(1,:)*1000)                  % Plot in milliseconds
legend("Reference Point", "NE", "NW", "SE", "SW")

%%
% figure
% latencyRate = diff(delay,1,2)/sampleTime;
% plot(time(1:end-1),latencyRate(1,:)*1e6)          % Plot in microseconds/second
% xlim([time(1) time(end-1)])
% title("First Satellite's Latency Rate vs. Time")
% xlabel("Simulation Time")
% ylabel("Latency Rate (\mus/s)")
% grid on

%%
% Emitted carrier frequency in Hz
fc = 2e9;
% Calculate Doppler shift and get additional information about Doppler. The
% Doppler information output contains relative velocity and Doppler rate.
% The dopplershift function internally performs access analysis and returns
% NaN whenever there is no access.
[fShift_reference,time_reference,dopplerInfo_reference] = dopplershift(sat,reference_point,Frequency=fc);
[fShift_NE,time_NE,dopplerInfo_NE] = dopplershift(sat,gs_NE,Frequency=fc);
[fShift_NW,time_NW,dopplerInfo_NW] = dopplershift(sat,gs_NW,Frequency=fc);
[fShift_SE,time_SE,dopplerInfo_SE] = dopplershift(sat,gs_SE,Frequency=fc);
[fShift_SW,time_SW,dopplerInfo_SW] = dopplershift(sat,gs_SW,Frequency=fc);


figure
% plot(time_reference,fShift_reference(1,:)/1e3)                        % Plot in kHz

% grid on
% hold on 

plot(time_NE, (fShift_NE(1,:)-fShift_reference)/1e3)   
xlim([time_reference(1) time_reference(end)])
title("Residual Doppler Shift")
xlabel("Simulation Time")
ylabel("Doppler Shift (kHz)") % Plot in kHz
grid on

hold on 
plot(time_NW, (fShift_NW(1,:)-fShift_reference)/1e3)                        % Plot in kHz

hold on 
plot(time_SE, (fShift_SE(1,:)-fShift_reference)/1e3)                        % Plot in kHz

hold on 
plot(time_SW, (fShift_SW(1,:)-fShift_reference)/1e3) % Plot in kHz
legend("NE", "NW", "SE", "SW")




figure
plot(time_reference(1:end-1), dopplerInfo_reference.DopplerRate(1,:))  % Plot in Hz/s
xlim([time_reference(1) time_reference(end-1)])
title("First Satellite's Doppler Rate vs. Time")
xlabel("Simulation Time")
ylabel("Doppler Rate (Hz/s)")
grid on

hold on 

plot(time_NE(1:end-1),dopplerInfo_NE.DopplerRate(1,:))  % Plot in Hz/s
hold on 

plot(time_NW(1:end-1),dopplerInfo_NW.DopplerRate(1,:))  % Plot in Hz/s
hold on 

plot(time_SE(1:end-1),dopplerInfo_SE.DopplerRate(1,:))  % Plot in Hz/s
hold on 

plot(time_SW(1:end-1),dopplerInfo_SW.DopplerRate(1,:))  % Plot in Hz/s

legend("Reference Point", "NE", "NW", "SE", "SW")