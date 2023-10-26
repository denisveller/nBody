%% Celestial Body simualtor - DO NOT ADJUST
tStepCelestial = 1;
nDays = 365;

tMax = 24*60*60*nDays;

%Julian date: 2458399.1507421, Artemis 1 Post Disposal

%initial positions given in km, coordinate 0,0 is the solar system
%barycenter

moon = massiveBody("moon",7.349 * 10^22, [1.446521618838699E+08 , 3.745053632575661E+07 , 1.165587060644478E+04], [-7.865719625981020E+00 , 2.774043016334739E+01 , 5.913011174721206E-02]);
earth = massiveBody("earth",5.97219*10^24, [1.450173573606565E+08 , 3.741219391969249E+07 , -1.350422362502851E+04], [-7.735342008869714E+00 , 2.880155923976655E+01 ,-1.641112159793678E-03]);
sun = massiveBody("sun",1.988500*10^30, [-1.497177482546914E+04 , 1.081890812458717E+06 ,-1.103830272471352E+04], [-1.310634574107055E-02 , 4.585129640636250E-03 , 3.275255759649712E-04]);
objArr = [sun,earth,moon];

% determines whether regenration is needed. Save and load data to avoid
% needing to run this. any changes to any above configuration parameters
% will triger a rerun, as will 'rawOutput' being missing from the workspace
if ~(exist('rawOutput','var') && length(rawOutput) == length(objArr)*3*(tMax/tStepCelestial + 1) && sun.position(1)==rawOutput(1))
    [rawOutput, rawOutputVelocity] = nBodyMatrixGen(tMax,tStepCelestial,objArr);
end

shaped = reshape(rawOutput,3,[])';
sunArr = shaped(1:3:end,:);
earthArr = shaped(2:3:end,:);
moonArr = shaped(3:3:end,:);

shapedV = reshape(rawOutputVelocity,3,[])';
sunVelArr = shapedV(1:3:end,:);
earthVelArr = shapedV(2:3:end,:);
moonVelArr = shapedV(3:3:end,:);

%% Satelite simulator
%configuration for running the satelite through existing body
%positions
tStep = 10; %must be an integer multiple of tStep from genFresh
nDays = 154.5; %must be less than or equal to number of days used to generate data

%configuration for satelite properties
mass = 24; %kg
%starting_orbit_alt = 200; %km
solar_rad_pressure_mag = 0; %N
solar_rad_pressure = [1 0 0]* solar_rad_pressure_mag;

tMax = 24*60*60*nDays; % turns input t into number of seconds


sunMass = sun.mass;
earthMass = earth.mass;
moonMass = moon.mass;

%% Legacy Orbit Setup
%this creates a perfectly circular orbit in the plane of the moon. at t=0,
%the satelite is directly under the moon

%earthRadius = 6371; %km
%orbit_radius = starting_orbit_alt + earthRadius; %km
%vcs = sqrt((earth.G*earth.mass)/(orbit_radius*1000))/1000; %km

%moonUnitPos = (moon.position-earth.position)/norm(moon.position-earth.position);
%moonUnitVel = (moon.velocity-earth.velocity)/norm(moon.velocity-earth.velocity);

% axisAngle = vrrotvec(moonUnitPos,moonUnitVel);
% axisAngle(4) = pi/2;
% rotMat = vrrotvec2mat(axisAngle);
% tangentialUnit = rotMat*moonUnitPos';
% tangentialUnit = tangentialUnit';

% starting_position = moonUnitPos*orbit_radius + earth.position/1000;
% starting_velocity = tangentialUnit*vcs + earth.velocity/1000;

%% Setup for a 2018 artemis 1 trajectory per https://ntrs.nasa.gov/api/citations/20205004731/downloads/SLS%20Earth-Moon%20Departure%20Trajectory%20Analysis_Draft4.pdf
a = 206959.1154;
e = 0.966351726;
i = 28.26180851; %degrees
RAAN = 35.37723043;%degrees
w = 41.23445507; %degrees
f = 118.0730252; %degrees
u = 398600.435436;

i = i * (pi/180);
RAAN =  RAAN * (pi/180);
w = w * (pi/180);
f = f * (pi/180);

tilt = -23.4;

[r,v] = elem_to_eci(a,e,i,RAAN,w,f,u);

R = [1 0 0; 0 cosd(tilt) -1*sind(tilt); 0 sind(tilt) cosd(tilt)];

r = (R*r')';
v = (R*v')';

% E = (e + cos(f))/(1+e*cos(f));
% M= E - e*sin(E);
% 
% M*180/pi


starting_position = r*1000 + earth.position;
starting_velocity = v*1000 + earth.velocity;

sat = secondaryBody('Sat',mass,starting_position/1000,starting_velocity/1000);

t = 0;
n = 1;

satPosVec = zeros(tMax/tStep,3);
satPosVec(1,:) = sat.position;

earthArrTstepChange = earthArr(1:tStep/tStepCelestial:end,:);
moonArrTstepChange = moonArr(1:tStep/tStepCelestial:end,:);
sunArrTstepChange = sunArr(1:tStep/tStepCelestial:end,:);

earthVelArrTstepChange = earthVelArr(1:tStep/tStepCelestial:end,:);
moonVelArrTstepChange = moonVelArr(1:tStep/tStepCelestial:end,:);
sunVelArrTstepChange = sunVelArr(1:tStep/tStepCelestial:end,:);



%% trajectory tuning

t_tcm_1 = 0*24*60*60 + 12*60*60 + 0*60;
t_tcm_2 = 35*24*60*60 + 0*60*60 + 0*60;
t_tcm_3 = 142*24*60*60 + 0*60*60 + 0*60;

t_tcm_4 = 96*24*60*60 + 12*60*60 + 0*60;
t_tcm_5 = 98*24*60*60 + 18*60*60 + 0*60;
t_tcm_6 = 136*24*60*60 + 0*60*60 + 0*60;

dv_tcm_1_prograde = 11.55;
dv_tcm_1_normal = 0;
dv_tcm_1_binormal = 0;

dv_tcm_2_prograde = 2.8;
dv_tcm_2_normal = 0;
dv_tcm_2_binormal = 0;

dv_tcm_3_prograde = 0;
dv_tcm_3_normal = -15.1;
dv_tcm_3_binormal = 0;

t_EDL_start_looking = 13300000 - 3600;
EDL_DONE = false;

dist = 10^10;
distE = 10^10;

% transient heating of entry vehicle
T_internal = 288.5457;
T_vec = [];

while t < tMax
    %computing distance info for timing perilune burns
    newDist = norm(sat.position - moonArrTstepChange(n,:)); %Insertion at closest aproach
    delta = dist - newDist;

    newDistE = norm(sat.position - earthArrTstepChange(n,:)); 
    deltaE = distE - newDistE;
    distEsurf = distE-6371000;

    %computing directions

    sunDir = (sat.position - sunArrTstepChange(n,:))/norm(sat.position-sunArrTstepChange(n,:));
    solar_rad_pressure = sunDir * solar_rad_pressure_mag;

    earthRelPrograde = (sat.velocity - earthVelArrTstepChange(n,:))/norm(sat.velocity - earthVelArrTstepChange(n,:));

    earthDir = (sat.position-earthArrTstepChange(n,:))/norm(sat.position-earthArrTstepChange(n,:));
    axisAngle = vrrotvec(earthRelPrograde,earthDir);
    axisAngle(4) = pi/2;
    rotMat = vrrotvec2mat(axisAngle);

    earthRelNormal = (rotMat*earthRelPrograde')';
    earthRelBinormal = cross(earthRelPrograde,earthRelNormal);

    moonRelPrograde = (sat.velocity - moonVelArrTstepChange(n,:))/norm(sat.velocity - moonVelArrTstepChange(n,:));

    moonDir = (sat.position - moonArrTstepChange(n,:))/norm(sat.position-moonArrTstepChange(n,:));
    axisAngle = vrrotvec(moonRelPrograde,moonDir);
    axisAngle(4) = pi/2;
    rotMat = vrrotvec2mat(axisAngle);

    moonRelNormal = (rotMat*moonRelPrograde')';
    moonRelBinormal = cross(moonRelPrograde,moonRelNormal);


    if (t == t_tcm_1) % TCM1
        vPrograde = dv_tcm_1_prograde*earthRelPrograde;
        vNormal = dv_tcm_1_normal*earthRelNormal;
        vBinormal = dv_tcm_1_binormal*earthRelBinormal;
        sat.velocity = sat.velocity + vPrograde + vNormal + vBinormal;
    end

    if (t == t_tcm_2) % TCM2
        vPrograde = dv_tcm_2_prograde*earthRelPrograde;
        vNormal = dv_tcm_2_normal*earthRelNormal;
        vBinormal = dv_tcm_2_binormal*earthRelBinormal;
        sat.velocity = sat.velocity + vPrograde + vNormal + vBinormal;
    end

    if (t == t_tcm_3) % TCM3
        vPrograde = dv_tcm_3_prograde*earthRelPrograde;
        vNormal = dv_tcm_3_normal*earthRelNormal;
        vBinormal = dv_tcm_3_binormal*earthRelBinormal;
        sat.velocity = sat.velocity + vPrograde + vNormal + vBinormal;
    end

    if (t > t_EDL_start_looking) && (distEsurf < 150000) && (EDL_DONE == false) % EDL time reporting
        EDL_DONE = true;
%         Spits out data
        t
        distEsurf
        norm(sat.velocity - earthVelArrTstepChange(n,:))
        axisAngle = vrrotvec(earthRelPrograde,earthDir);
        axisAngle(4)*(180/pi)-90
    end

    if (t > 13309860 - 3600)&&(t < 13309860)

        %pre EDL transient thermal analysis

        sigma = 5.670374419e-8;
        Se = 1366.1;

        k = 0.167; %W/mk, PICA,https://ntrs.nasa.gov/api/citations/20110014590/downloads/20110014590.pdf
        thickness = 0.015; %m,
        heatSheildDensity = 0.280 * (1/1000) * (100^3);
        heatSheildCp = 900; %https://www.jstage.jst.go.jp/article/tjsass/56/3/56_T-12-17/_pdf

        alt = distEsurf/1000;

        alpha_TPS = 0.85;
        epsilon_TPS = 0.85; %see EDL code
        a_tps = (12884.998121 + 26235.796147)/(1000^2);
        
        a_sup = 0.0;
        
        alpha_backside = 0.13;
        epsilon_backside = 0.9; %white paint per slides
        a_backside = (30673.258889 + 15160.15423)/(1000^2) + a_sup;
        
        Q_internal = 2; %arbitrary 2 watt heat load from basic electronics - likley high
        
        earth_vf = 0.5*(1-cos(asin(6371/(6371+alt))));
        T_earth = 273.15 + 15;
        is_day = 1;
        
        T_EV = (Q_internal+ is_day*(a_tps*alpha_TPS*Se*0.5 + a_backside*alpha_backside*Se*0.5) + sigma*(T_earth^4)*earth_vf*(epsilon_TPS*a_tps + epsilon_backside*a_backside))/(sigma*(epsilon_backside*a_backside + epsilon_TPS*a_tps));
        T_EV = (T_EV^(1/4));
        
        q_dot_cond = conductionOut(T_EV,T_internal,k);
        T_internal = T_internal + (2*q_dot_cond*tStep)/(thickness*heatSheildDensity*heatSheildCp);
        T_vec(end+1,:) = [T_internal, T_EV];
    end




    sat = sat.netAcceleration([earthArrTstepChange(n,:) earthMass; moonArrTstepChange(n,:) moonMass; sunArrTstepChange(n,:) sunMass],solar_rad_pressure);
    sat = sat.integrate(tStep);
    n = n + 1;

    satPosVec(n,:) = sat.position;
    t = t + tStep;

    dist = newDist;
    distE=newDistE;
end

%% Display code
factor = 6;

%1/factor data points are used - default is 6, coresponding to plotting
%data 1x a minute. does not affect underlying computation

satArrL = satPosVec(1:factor:end,:);
moonArrTstepChange = moonArrTstepChange(1:factor:end,:);
earthArrTstepChange = earthArrTstepChange(1:factor:end,:);
sunArrTstepChange = sunArrTstepChange(1:factor:end,:);

moonArrTstepChange = moonArrTstepChange(1:length(satArrL),:);
earthArrTstepChange = earthArrTstepChange(1:length(satArrL),:);
sunArrTstepChange = sunArrTstepChange(1:length(satArrL),:);

%gateway
gateway_perilune = 3*10^6 + 1.74*10^6;
gateway_apolune = 70*10^6 + 1.74*10^6;

a = (gateway_apolune + gateway_perilune)/1000;
e = (gateway_apolune - gateway_perilune)/(gateway_apolune+gateway_perilune);
i = pi/2;
RAAN = 0;
w = pi/2;
u = 4902.800066;
gateway_pos_arr = [];
num = 100;
for j = 1:num
    f = (2*pi*j)/num;
    [gateway_pos_arr(j,:),~] = elem_to_eci(a,e,i,RAAN,w,f,u);
end
gateway_pos_arr = gateway_pos_arr * 1000;
Xg = gateway_pos_arr(:,1);
Yg = gateway_pos_arr(:,2);
Zg = gateway_pos_arr(:,3);


%ECI view
moonArrP1 = moonArrTstepChange - earthArrTstepChange;
satArrP1 = satArrL - earthArrTstepChange;
sunArrP1 = sunArrTstepChange - earthArrTstepChange;


X = moonArrP1(:,1);
Y = moonArrP1(:,2);
Z = moonArrP1(:,3);


X1 = satArrP1(:,1);
Y1 = satArrP1(:,2);
Z1 = satArrP1(:,3);

X_sun = sunArrP1(:,1);
Y_sun = sunArrP1(:,2);
Z_sun = sunArrP1(:,3);

[Xs, Ys, Zs] = sphere;


figure
hold on
plot3(X,Y,Z)
plot3(X1,Y1,Z1)
plot3(X_sun,Y_sun,Z_sun)
surf(Xs*6371000,Ys*6371000,Zs*6371000)
title("ECI")
axis([-2*10^9, 2*10^9,-2*10^9,2*10^9,-2*10^9, 2*10^9])


% T_time = 1:length(T_vec(:,1));
% T_time = T_time * tStep;
% figure
% hold on
% plot(T_time, T_vec(:,1))
% plot(T_time, T_vec(:,2))
% hold off
% title("Transient Heating of Entry Vehicle")


%MCI view
% satArrP2 = satArrL - moonArrTstepChange;
% 
% X1 = satArrP2(:,1);
% Y1 = satArrP2(:,2);
% Z1 = satArrP2(:,3);

% figure
% hold on
% plot3(X1,Y1,Z1)
% %plot3(Xg,Yg,Zg)
% axis([-2*10^8, 2*10^8,-2*10^8, 2*10^8,-2*10^8, 2*10^8])
% 
% surf(Xs*1740000,Ys*1740000,Zs*1740000)
% title("MCI")

function [qDotCond] = conductionOut(T, Tinterior, k)
    qDotCond = k*(T-Tinterior);
end
function T2_new = thermalEstimator(T2,T1,emmisivity,k,q_dot_in,tStep,radius,thickness,heatSheildDensity,heatSheildCp)
    q_dot_cond = conductionOut(T2,T1,k);
    q_dot_rad = radiationOut(emmisivity,T2);
    if q_dot_in < 1e3
        T2_new = T2 - (2*q_dot_rad*tStep)/(thickness*heatSheildDensity*heatSheildCp); %uses same linear assumption as T1 heating - probably really bad to use it twice but hey.
    else
        qRatio = q_dot_in/(q_dot_cond+q_dot_rad);
        if abs(1-qRatio) > 0.01
            T2 = T2*((qRatio)^0.1);
            T2_new = thermalEstimator(T2,T1,emmisivity,k,q_dot_in,tStep,radius,thickness,heatSheildDensity,heatSheildCp);
        else
            T2_new = T2;
        end
    end
end

