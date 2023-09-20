%% Celestial Body simualtor - DO NOT ADJUST
tStepCelestial = 1;
nDays = 180;

tMax = 24*60*60*nDays;

%initial positions given in km, coordinate 0,0 is the solar system
%barycenter

moon = massiveBody("moon",7.349 * 10^22, [-2.649907633791349E+07 , 1.451115329759790E+08 , 1.920803579456359E+04], [-3.031955180389345E+01 ,-4.429909321012257E+00 , 8.719919857097347E-02]);
earth = massiveBody("earth",5.97219*10^24, [-2.682456095823074E+07 , 1.448975704171931E+08 , 2.413001115289330E+04], [-2.981482351984199E+01 ,-5.295918918751036E+00 ,-5.493700746179719E-04]);
sun = massiveBody("sun",1.988500*10^30, [-1.354630417531006E+06, 1.420546299731827E+04, 3.143940434895225E+04], [1.636302711538109E-03, -1.558286309774412E-02, 8.961626340290262E-05]);
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
nDays = 120; %must be less than or equal to number of days used to generate data

%configuration for satelite properties
mass = 20; %kg
starting_orbit_alt = 200; %km
solar_rad_pressure = 0; %N

tMax = 24*60*60*nDays; % turns input t into number of seconds


sunMass = sun.mass;
earthMass = earth.mass;
moonMass = moon.mass;

%this creates a perfectly circular orbit in the plane of the moon. at t=0,
%the satelite is directly under the moon

earthRadius = 6371; %km
orbit_radius = starting_orbit_alt + earthRadius; %km
vcs = sqrt((earth.G*earth.mass)/(orbit_radius*1000))/1000; %km

moonUnitPos = (moon.position-earth.position)/norm(moon.position-earth.position);
moonUnitVel = (moon.velocity-earth.velocity)/norm(moon.velocity-earth.velocity);

axisAngle = vrrotvec(moonUnitPos,moonUnitVel);
axisAngle(4) = pi/2;
rotMat = vrrotvec2mat(axisAngle);
tangentialUnit = rotMat*moonUnitPos';
tangentialUnit = tangentialUnit';

starting_position = moonUnitPos*orbit_radius + earth.position/1000;
starting_velocity = tangentialUnit*vcs + earth.velocity/1000;



sat = secondaryBody('Sat',mass,starting_position,starting_velocity);

t = 0;
n = 1;

satPosVec = zeros(tMax/tStep,3);
satPosVec(1,:) = sat.position;

earthArrL = earthArr(1:tStep/tStepCelestial:end,:);
moonArrL = moonArr(1:tStep/tStepCelestial:end,:);
sunArrL = sunArr(1:tStep/tStepCelestial:end,:);

% trajectory tuning
t_tli = 52*60;    
dv_tli = 3197.650085; %m/s

t_moon_insertion = 110*24*60*60 + 3*60*60 + 15*60;%s
dv_mi = 120; %m/s
vDir = [1 0.1 0.3];
vDir = vDir/norm(vDir);
vAdd = vDir*dv_mi;
while t < tMax
    if (t == t_tli) % TLI
        vStart = sat.velocity - earthVelArr(n,:);
        vDir = vStart/norm(vStart);
        vMagNew = dv_tli + norm(vStart);
        vNew = vMagNew * vDir;
        sat.velocity = earthVelArr(n,:) + vNew;
    end
    if (t == t_moon_insertion) % NRHO Insertion
        sat.velocity = sat.velocity + vAdd;
    end
    sat = sat.netAcceleration([earthArrL(n,:) earthMass; moonArrL(n,:) moonMass; sunArrL(n,:) sunMass],solar_rad_pressure);
    sat = sat.integrate(tStep);
    
    n = n + 1;

    satPosVec(n,:) = sat.position;
    t = t + tStep;
end

%% Display code
factor = 6;

%1/factor data points are used - default is 6, coresponding to plotting
%data 1x a minute. does not affect underlying computation

satArrL = satPosVec(1:factor:end,:);
moonArrL = moonArrL(1:factor:end,:);
earthArrL = earthArrL(1:factor:end,:);
moonArrL = moonArrL(1:length(satArrL),:);
earthArrL = earthArrL(1:length(satArrL),:);


%ECI view
moonArrP1 = moonArrL - earthArrL;
satArrP1 = satArrL - earthArrL;

X = moonArrP1(:,1);
Y = moonArrP1(:,2);
Z = moonArrP1(:,3);


X1 = satArrP1(:,1);
Y1 = satArrP1(:,2);
Z1 = satArrP1(:,3);

[Xs, Ys, Zs] = sphere;


figure
hold on
plot3(X,Y,Z)
plot3(X1,Y1,Z1)
surf(Xs*6371000,Ys*6371000,Zs*6371000)
title("ECI")
axis([-3*10^9, 3*10^9,-3*10^9,3*10^9,-3*10^9, 3*10^9])

%MCI view
satArrP2 = satArrL - moonArrL;

X1 = satArrP2(:,1);
Y1 = satArrP2(:,2);
Z1 = satArrP2(:,3);

figure
hold on
plot3(X1,Y1,Z1)
axis([-.75*10^8, .75*10^8,-.75*10^8, .75*10^8,-.75*10^8, .75*10^8])

surf(Xs*1740000,Ys*1740000,Zs*1740000)
title("MCI")


