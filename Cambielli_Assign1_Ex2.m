%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     1
% Exercise:       2 - Impulsive Guidance
% Author:         Alessandro Cambielli

%% REQUEST 1

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 2 Request 1  \n \n')

% integration frame and center
center = 'SSB'; % solar system barycenter
frame = 'ECLIPJ2000'; % ecliptic plane

% define close approach time window
closeapproach = datetime(2029,1,1):hours(1):datetime(2029,7,31);
date = cspice_str2et(char(closeapproach));

% asteroid_label = '99942 Apophis';
initial_epoch = '2029-Jan-01 00:00:00.0000 TDB';
final_epoch = '2029-Jul-31 00:00:00.0000 TDB';

pos_apophis = cspice_spkpos('20099942',date,frame,'NONE',center);
pos_earth = cspice_spkpos('Earth',date,frame,'NONE',center);
pos_moon = cspice_spkpos('Moon',date,frame,'NONE',center);
pos_sun = cspice_spkpos('Sun',date,frame,'NONE',center);

% Earth-Apophis distance
rr_earth_apophis = pos_apophis - pos_earth;
dist_km1 = sqrt(sum(rr_earth_apophis.^2,1));
dist_AU1 = cspice_convrt(dist_km1,'km','au');

% 'Relative distance Earth-Apophis (@ECLIPJ2000)'
figure(1)
subplot(1,3,1);
hold on
grid on
plot(closeapproach, dist_AU1,'r','LineWidth', 1.5)
% Plot settings
set(gca,'FontSize',12)
legend('Relative distance Earth-Apophis','Location','northeast','FontSize',12);
xlabel('$Epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$Distance$ [AU]','Interpreter','latex','FontSize', 20)

% Moon-Apophis distance
rr_moon_apophis = pos_apophis - pos_moon;
dist_km2 = sqrt(sum(rr_moon_apophis.^2,1));
dist_AU2 = cspice_convrt(dist_km2,'km','au');

% 'Relative distance Moon-Apophis (@ECLIPJ2000)'
subplot(1,3,2);
hold on
grid on
plot(closeapproach, dist_AU2,'r','LineWidth', 1.5)
% Plot settings
set(gca,'FontSize',12)
legend('Relative distance Moon-Apophis','Location','northeast','FontSize',12);
xlabel('$Epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$Distance$ [AU]','Interpreter','latex','FontSize', 20)

% Sun-Apophis distance
rr_sun_apophis = pos_apophis - pos_sun;
dist_km3 = sqrt(sum(rr_sun_apophis.^2,1));
dist_AU3 = cspice_convrt(dist_km3,'km','au');

% 'Relative distance Sun-Apophis (@ECLIPJ2000)'
subplot(1,3,3);
hold on
grid on
plot(closeapproach, dist_AU3,'r','LineWidth', 1.5)
ylim([0.86 1.15])
% Plot settings
set(gca,'FontSize',12)
legend('Relative distance Sun-Apophis','Location','northeast','FontSize',12);
xlabel('$Epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$Distance$ [AU]','Interpreter','latex','FontSize', 20)

% Earth-Apophis-Sun angle
apophis_earth = pos_apophis - pos_earth;
apophis_sun = pos_apophis - pos_sun;
earth_apophis_sun_angle = acosd(dot(apophis_earth,apophis_sun)./...
    (vecnorm(apophis_sun,2,1).*vecnorm(apophis_earth,2,1) ));

% 'Earth-Apophis-Sun angle @ECLIPJ2000'
figure(2)
hold on
grid on
plot(closeapproach,earth_apophis_sun_angle,'r','LineWidth', 2)
% Plot settings
set(gca,'FontSize',12)
legend('Earth-Apophis-Sun angle','Location','northeast','FontSize',12);
xlabel('$Epoch$ [date]','Interpreter','latex','FontSize', 20)
ylabel('$Angle$ [deg]','Interpreter','latex','FontSize', 20)


% ground-tracks computation over the closest approach window
min_time = fmincon(@apophisEarthclosest,date(round(length(date)/2)),...
    [],[],[],[],date(1),date(end),[]);
date_min = cspice_et2utc(min_time,'C',1);
fprintf('   closest approach date: %s UTC \n',date_min)
timewindow = (min_time-6*3600):(min_time+6*3600);
rectcoord = cspice_spkpos('20099942',timewindow,'IAU_EARTH','NONE','EARTH');

% retrieve Earth's radii
Re = cspice_bodvrd('399','RADII',3);

% compute flatness of the Earth 
fe = (Re(1)-Re(3))/Re(1);

% obtain longitudes and latitudes
[lon,lat,~] = cspice_recgeo(rectcoord, Re(1), fe);

% conversion from radiants to degrees
lon_deg = rad2deg(lon); 
lat_deg = rad2deg(lat);
wrap = wrapTo180(lon_deg);

% 'Aphopis Ground-tracks'
figure(3)
img = imread('earthsurface.jpg');  % Load a sample image
img = flipud(img);
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap,lat_deg,'r','LineStyle', 'none', 'Marker', '.')
hold on
plot(wrap(1,1),lat_deg(1,1),'o',wrap(end),lat_deg(end),'o','LineWidth',2)
xlim([-180 180]); ylim([-90 90]); grid on;
xticks(-180:30:180);
yticks(-90:30:90);
% Plot settings
set(gca,'FontSize',12)
legend('Ground-Track', 'Start', 'End','Location','northeast',...
    'FontSize',12);
xlabel('$Longitude$ [deg]','Interpreter','latex','FontSize', 20)
ylabel('$Latitude$ [deg]','Interpreter','latex','FontSize', 20)

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REQUEST 3

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 2 Request 3  \n \n')

% If the solver stops without finding a proper solution, run again.
% Maximum number of function evaluations or maximum number of iterations
% has probably been reached.
% That means that the initial condition generated is awful.

% Setup for the plots:
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% select integration frame string (SPICE naming convention)
center = 'SSB';
frame = 'ECLIPJ2000'; % ecliptic plane, different from J2000 (equatorial)!

bodies = nbody_init(labels);

% retrieve Earth's radii
Re = cspice_bodvrd('399','RADII',3);

% INITIAL CONDITIONS
% launch window
LWO = cspice_str2et(char(datetime(2024,10,1))); % Launch Window Open
LWC = cspice_str2et(char(datetime(2025,2,1))); % Launch Window Close

% deep space manoeuvre window
DSMO = cspice_str2et(char(datetime(2024,10,1)+calmonths(6))); % DSM Open
DSMC = cspice_str2et(char(datetime(2025,2,1)+calmonths(18))); % DSM Close

% impact window
IMPO = cspice_str2et(char(datetime(2028,8,1))); % Impact Window Open
IMPC = cspice_str2et(char(datetime(2029,2,28))); % Impact Window Close

lb = [-1e12*ones(18,1); LWO; DSMO; IMPO];
ub = [1e12*ones(18,1); LWC; DSMC; IMPC];

% define the options for fmincon and the ode
optoptions = optimoptions('fmincon','Display','iter', 'MaxIter', 3000,...
    'MaxFunEvals', 3000);
odeoptions = odeset('RelTol',1e-11,'AbsTol',1e-12);

rng("shuffle")
gs = GlobalSearch;

% define the initial conditions
x0 = zeros(21,1);
x0(19) = LWO+(LWC-LWO)*rand(1,1); 
x0(20) = DSMO+(DSMC-DSMO)*rand(1,1); 
x0(21) = IMPO+(IMPC-IMPO)*rand(1,1); 

x0(1:6) = cspice_spkezr('Earth', x0(19), frame, 'NONE', center);

[x2guess, ~, ~] = keplerian_propagator(x0(19), x0(1:6), x0(20), 'Sun');
x0(7:12) = x2guess;

x3guess = cspice_spkezr('20099942', x0(21), frame, 'NONE', center);
x0(13:18) = x3guess;

% optimization
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',@(x)objectiveF(x,frame,center,bodies),'lb',lb,...
    'ub',ub,'options',optoptions,'nonlcon',@(x)nonlcon(x));

[states,fval] = run(gs,problem);

% RESULTS
% retrieve the state of the Earth
state_E = cspice_spkezr('Earth',states(19),frame,'NONE',center);

% compute the dvs of the mission phases 
dv0 = states(4:6) - state_E(4:6);
dv0_norm = norm(dv0);

[state_SC, ~, ~] = keplerian_propagator(states(19), states(1:6), states(20), 'Sun');

dv_dsm = states(10:12) - state_SC(4:6);
dv_dsm_norm = norm(dv_dsm);

% compute the total dv of the mission (max = 5 km/s)
dvtot = dv0_norm + dv_dsm_norm;

% retrieve the launch, dsm, impact and closest approach dates
t0 = cspice_et2utc(states(19),'C',1);
tdsm = cspice_et2utc(states(20),'C',1);
timp = cspice_et2utc(states(21),'C',1);

state_afterimp = cspice_spkezr('20099942',states(21),frame,'NONE',center) +...
    [zeros(3,1); states(16:18)*5e-5];
opt = odeset('RelTol',1e-8,'AbsTol',1e-9,'Events',@(t,s) myEvent(t,s));
[t,s] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame),[states(21) states(21)*10],...
    state_afterimp,opt);
tclosest = cspice_et2utc(t(end),'C',1);

% retrieve the closest approach distance
mindistance = abs(fval)/Re(1);

% Display results
fprintf('Solutions:\n\n')
fprintf('   departure date: %s UTC \n',t0)
fprintf('   dsm date: %s UTC \n', tdsm)
fprintf('   impact date: %s UTC \n', timp)
fprintf('   closest approach date: %s UTC \n', tclosest)
fprintf('   minimum distance: %f earth radii \n', mindistance)
fprintf('   deltav: %f km/s \n', dvtot)
fprintf('------------------------------------------------------------------------------------\n')

GM = cspice_bodvrd('Sun', 'GM', 1);
points = 100000;

t1 = linspace(states(19),states(20),points);
[tt1, xx1] = ode78(@(t,x) keplerian_rhs(t,x,GM), t1, states(1:6), odeoptions);

t2 = linspace(states(20),states(21),points);
[tt2, xx2] = ode78(@(t,x) keplerian_rhs(t,x,GM), t2, states(7:12), odeoptions);

tapophis = linspace(states(19),states(21),points);
posApophis1 = cspice_spkpos('20099942',tapophis,frame,'NONE',center);

t3 = linspace(states(21),t(end),points);
[tt3, xx3] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame), t3, state_afterimp,...
    odeoptions);

tearth = linspace(states(19),states(21),points);
posEarth = cspice_spkpos('Earth',tearth,frame,'NONE',center);

% 'Asteroid impacting mission (@Sun ECLIPTIC J2000)'
figure(4)
plot3(xx1(:,1), xx1(:,2), xx1(:,3), 'b', 'LineWidth',2)
hold on ; grid on ; 
plot3(xx2(:,1), xx2(:,2), xx2(:,3), 'c', 'LineWidth',2)
plot3(posApophis1(1,:), posApophis1(2,:), posApophis1(3,:), 'k', 'LineWidth',2)
plot3(xx3(:,1), xx3(:,2), xx3(:,3), 'm', 'LineWidth',2)
plot3(posEarth(1,:), posEarth(2,:), posEarth(3,:), 'MarkerFaceColor',...
    '#77AC30', 'LineWidth',1)
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#FFFF00') % Sun
scatter3(xx1(1,1), xx1(1,2), xx1(1,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#0072BD') % launch
scatter3(xx1(end,1), xx1(end,2), xx1(end,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#4DBEEE') % dsm
scatter3(xx3(1,1), xx3(1,2), xx3(1,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#7E2F8E') % impact
scatter3(xx3(end,1), xx3(end,2), xx3(end,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#A2142F') % closest approach
ylim([-2e8 2e8])
view(2)
% Plot settings
set(gca,'FontSize',12)
legend('first segment', 'second segment', 'apophis before',...
    'apophis after','Earth orbit','Sun','Launch','DSM','Impact',...
    'Closest encounter','Location','bestoutside','FontSize',12);
xlabel('$x$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [km]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [km]','Interpreter','latex','FontSize', 20)

% 'Asteroid impacting mission (@Sun ECLIPTIC J2000)'
figure(5)
plot3(xx1(:,1), xx1(:,2), xx1(:,3), 'b', 'LineWidth',2)
hold on ; grid on ; 
plot3(xx2(:,1), xx2(:,2), xx2(:,3), 'c', 'LineWidth',2)
plot3(posApophis1(1,:), posApophis1(2,:), posApophis1(3,:), 'k', 'LineWidth',2)
plot3(xx3(:,1), xx3(:,2), xx3(:,3), 'm', 'LineWidth',2)
plot3(posEarth(1,:), posEarth(2,:), posEarth(3,:), 'MarkerFaceColor',...
    '#77AC30', 'LineWidth',1)
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#FFFF00') % Sun
scatter3(xx1(1,1), xx1(1,2), xx1(1,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#0072BD') % launch
scatter3(xx1(end,1), xx1(end,2), xx1(end,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#4DBEEE') % dsm
scatter3(xx3(1,1), xx3(1,2), xx3(1,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#7E2F8E') % impact
scatter3(xx3(end,1), xx3(end,2), xx3(end,3),'filled', 'SizeData', 120,...
    'MarkerFaceColor','#A2142F') % closest approach
view(3)
% Plot settings
set(gca,'FontSize',12)
xlabel('$x$ [km]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [km]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [km]','Interpreter','latex','FontSize', 20)

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS

function [bodies] = nbody_init(labels)
% The function returns a cell array populated with structures containing 
% the body label and the associated gravitational constant.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   labels : cell-array with object labels [1,n] 
%
% Outputs:
%   bodies : cell-array with struct elements containing the following [1,n] 
%                  bodies{i}.name -> body label
%                  bodies{i}.GM -> gravitational constant [km^3/s^2]
%

% Initialize cell array bodies
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)

    bodies{i}.name = labels{i};
    bodies{i}.GM = cspice_bodvrd(labels{i}, 'GM', 1);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [distance] = apophisEarthclosest(t)
% The function is the target to find the time of Apophis' minimum approach
% with respect to the Earth.
% 
% Author: Alessandro Cambielli
%
% Inputs:
%   t : times [1,1] 
%
% Outputs:
%   distance : Apophis-Earth distance [1,1]
%

distance = cspice_spkpos('20099942',t,'IAU_EARTH','NONE','EARTH');
distance = vecnorm(distance);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = keplerian_rhs(~, x, GM)
%KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
% Author
%   Name: ALESSANDRO 
%   Surname: CAMBIELLI
%
% Inputs:
%   t   : [ 1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [ 1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6); 

% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xf, tt, xx] = keplerian_propagator(et0, x0, etf, attractor)
% The function propagates the keplerian motion around an attractor body
%
% Author: Alessandro Cambielli
%
% Inputs:
%   et0 : initial time [1,1] 
%   x0 : initial state conditions [6,1] 
%   etf : final time [1,1] 
%   attractor : name of the attractor body [string] 
%
% Outputs:
%   xf : final state of the body [6,1]
%   tt : integration time vector [n]
%   xx : state of the body at all points [n]
%

if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

options = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [et0 etf], x0, options);

% Extract final state vector
xf = xx(end,1:6)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = objectiveF(x,frame,center,bodies)
% The function sets the function to minimize in fmincon
%
% Author: Alessandro Cambielli
%
% Inputs:
%   x : variables [21,1] 
%   frame : reference frame [string] 
%   center : considered barycenter [string] 
%   bodies : bodies to be considered in the analysis [string] 
%
% Outputs:
%   J : value to be minimized [1,1]
%

t_imp = x(21);

state_afterimp = cspice_spkezr('20099942',t_imp,frame,'NONE',center) +...
    [zeros(3,1); x(16:18)*5e-5];
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@(t,s) myEvent(t,s));
[t,s] = ode78(@(t,s) nbody_rhs(t,s,bodies,frame),[t_imp t_imp*20],...
    state_afterimp,options);

% to maximize a function, minimize its negative
J = -norm(s(end,1:3)-cspice_spkpos('EARTH',t(end),frame,'NONE',center)');  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, Ceq] = nonlcon(x)
% The function sets the equality and inequality constraints for fmincon
%
% Author: Alessandro Cambielli
%
% Inputs:
%   x : variables [21,1]  
%
% Outputs:
%   C : inequality constraints [1,1]
%   Ceq : equality constraints [15,1]
%

t0 = x(19);
tdsm = x(20);
timp = x(21);
x1 = x(1:6);
x2 = x(7:12);
x3 = x(13:18);

dvmax = 5;

[xf1, ~, ~] = keplerian_propagator(t0, x1, tdsm, 'Sun');
[xf2, ~, ~] = keplerian_propagator(tdsm, x2, timp, 'Sun');

s_earth = cspice_spkezr('EARTH',t0,'ECLIPJ2000','NONE','SSB');
s_apophis = cspice_spkpos('20099942',timp,'ECLIPJ2000','NONE','SSB');

r0 = s_earth(1:3);
v0 = s_earth(4:6);
dv1 = x1(4:6) - v0;
dv2 = x2(4:6) - xf1(4:6);

Ceq(1:3) = x1(1:3) - r0; %psi initial
Ceq(4:6) = xf1(1:3) - x2(1:3); % z1
Ceq(7:12) = xf2 - x3; % z2
Ceq(13:15) = x3(1:3) - s_apophis; %psi final

% inequality constraints
C = norm(dv1) + norm(dv2) - dvmax;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dxdt] = nbody_rhs(t, x, bodies, frame)
% The function evaluates the right-hand-side of a newtonian N-body
% propagator.
% The integration centre is the Solar-System-Barycentre (SSB) and only
% Newtonian gravitational accelerations are considered.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   t : ephemeris time (ET SPICE), seconds past J2000 (TDB) [1,1]
%   x : cartesian state vector wrt Solar-System-Barycentre [6,1] 
%   bodies : cell array of bodies [string]
%   frame : reference frame [string] 
%
% Outputs:
%   dxdt : RHS, newtonian gravitational acceleration only [6,1]
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position from state x
rr_ssb_obj = x(1:3);

% Loop over all bodies
for i = 1:length(bodies)

    % Retrieve position and velocity of i-th celestial body wrt Solar
    % System Barycentre in inertial frame
    rv_ssb_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'SSB');

    % Extract object position wrt. i-th celestial body
    rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3);

    % Compute square distance and distance
    dist2 = dot(rr_body_obj, rr_body_obj);
    dist = sqrt(dist2);

    % Compute the gravitational acceleration using Newton's law
    aa_grav = - bodies{i}.GM * rr_body_obj / (dist*dist2);

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + aa_grav;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction] = myEvent(t,x)
% The event function detects the closest approach between Apophis and the
% Earth after the impact with the spacecraft.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   x : state vector [6,1]  
%
% Outputs:
%   value : expression describing the event [1,1]
%   isterminal : 1 if the integration is to terminate when the event occurs
%   direction : 1  locates only zeros where the event function is increasing
%

x_e = cspice_spkezr('EARTH',t,'ECLIPJ2000','NONE','SSB');
distance = x(1:3) - x_e(1:3);
velocity = x(4:6) - x_e(4:6);
value = 2*dot(distance,velocity);
isterminal = 1;         
direction  = 1;         

end

















