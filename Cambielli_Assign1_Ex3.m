%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     1
% Exercise:       3 - Continuous Guidance
% Author:         Alessandro Cambielli

%% REQUEST 2

% Clear memory workspace and set path
clearvars; close all; clc 

fprintf('Exercise 3 Request 2  \n \n')

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

% Data 
t0 = '2023-05-28-14:13:09.000 UTC'; 
et0 = cspice_str2et(t0); % Ephemeris time [s] 
xx0 = cspice_spkezr('Earth',et0,'ECLIPJ2000','NONE','Sun'); % Initial state vector [km,km/s]
rr0 = xx0(1:3); % Earth's initial position [km]                                                   
vv0 = xx0(4:6); % Earth's initial velocity [km/s]                                                      

m0 = 1000;  % Spacecraft mass [kg]
Tmax = 800e-6;  % Maximum thrust  [kN]
Isp = 3120;  % Specific impulse [s]
mu = cspice_bodvrd('Sun','GM',1);  % Sun gravitational parameter [km^3/s^2]
g0 = 9.80665e-3;  % Gravitational constant [km/s^2]

% Adimensionalizing unit
LU = cspice_convrt(1,'AU','km'); % Length unit (Astronomical unit) [km]
MU = m0; % Mass unit [kg]
TU = sqrt(LU^3/mu); % Time unit [s] 

% Adimensional parameters
data.mu = 1; % Adimensional Sun parameter[-]
data.rr0 = rr0/LU; % Adimensional initial position [-]
data.vv0 = vv0/LU*TU; % Adimensional initial velocity [-] 
data.Tmax = Tmax*TU^2/MU/LU; % Adimensional thrust [-]
data.Isp = Isp/TU; % Adimensional specific impulse [-]
data.g0 = g0*TU^2/LU; % Adimensional sea level acceleration constant [-]
data.m0 = m0/MU; % Adimensional mass [-]
data.et0 = et0/TU; % Adimensional initial time [-]

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REQUEST 3

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 3 Request 3  \n \n')

% Data 
t0 = '2023-05-28-14:13:09.000 UTC'; 
et0 = cspice_str2et(t0); % Ephemeris time [s] 
xx0 = cspice_spkezr('Earth',et0,'ECLIPJ2000','NONE','Sun'); % Initial state vector [km,km/s]
rr0 = xx0(1:3); % Earth's initial position [km]                                                   
vv0 = xx0(4:6); % Earth's initial velocity [km/s]                                                      

m0 = 1000;  % Spacecraft mass [kg]
Tmax = 800e-6;  % Maximum thrust  [kN]
Isp = 3120;  % Specific impulse [s]
mu = cspice_bodvrd('Sun','GM',1);  % Sun gravitational parameter [km^3/s^2]
g0 = 9.80665e-3;  % Gravitational constant [km/s^2]

% Adimensionalizing unit
LU = cspice_convrt(1,'AU','km'); % Length unit (Astronomical unit) [km]
MU = m0; % Mass unit [kg]
TU = sqrt(LU^3/mu); % Time unit [s] 

% Adimensional parameters
data.mu = 1; % Adimensional Sun parameter[-]
data.rr0 = rr0/LU; % Adimensional initial position [-]
data.vv0 = vv0/LU*TU; % Adimensional initial velocity [-] 
data.Tmax = Tmax*TU^2/MU/LU; % Adimensional thrust [-]
data.Isp = Isp/TU; % Adimensional specific impulse [-]
data.g0 = g0*TU^2/LU; % Adimensional sea level acceleration constant [-]
data.m0 = m0/MU; % Adimensional mass [-]
data.et0 = et0/TU; % Adimensional initial time [-]

% If the solver stops without finding a proper solution, run again.
% Maximum number of function evaluations or maximum number of iterations
% has probably been reached.
% That means that the initial condition generated is awful.

error = 1;
tol = 1e-9; % Set error tolerance for exit criteria

% Initial guess for the fsolve problem 
tf_guess = (2*pi-0).*rand(1) + 0;
lambda0 = zeros(7,1);  
lambda0(1:7) = (20-(-20)).*rand(7,1) + (-20);
etf = tf_guess + data.et0;

% Optimal lambda and tf finder
while error > tol

    % shooting algorithm
    options = optimset('Algorithm', 'levenberg-marquardt','Display', ...
                        'iter','TolFun',1e-15, 'MaxIter',2000, 'MaxFunEvals',2000); 
    [Y_opt,fval] = fsolve(@(Y) shooting(Y,data), [lambda0; etf], options);

    % Check for the error 
    error = norm(fval);

    lambda0 = Y_opt(1:7);
    etf = Y_opt(8);

end

optimal_lambda0 = lambda0;
optimal_etf = etf;

% Time of Flight [days]
ToF = (optimal_etf - data.et0)*TU/cspice_spd;

% Final state error with respect to Venus 
residual = fval;

err_rr = norm(residual(1:3)) * LU ; % Final potision error[km]
err_vv = norm(residual(4:6)) * LU/TU*1e3 ;  % Final velocity error [m/s]

% Final time string
Final_date = cspice_et2utc(optimal_etf*TU,'C',1); 

% Display results
fprintf('Solutions of the time-optimal Earth-Venus transfer considering 800 mN continuous thrust are:\n\n')
fprintf('   λ = [');
fprintf('%f ',optimal_lambda0(1:7))
fprintf(']''\n')
fprintf('   Arrival date: %s UTC \n',Final_date)
fprintf('   Time of flight: %f days\n', ToF)
fprintf('   Position error: %f km\n', err_rr)
fprintf('   Velocity error: %e m/s\n', err_vv)
fprintf('------------------------------------------------------------------------------------\n')

% Integrating state and the costate using the optimal initial condition  
options = odeset('AbsTol',1e-13 , 'RelTol', 1e-13); 
tspan = linspace(data.et0,optimal_etf,1000); 
Y0 = [data.rr0; data.vv0 ; data.m0; optimal_lambda0]; 
[time, Y] = ode113(@(t,Y) TPBVP(t,Y,data), tspan, Y0, options);

% Settings for the plots
points = 100000;

% Propagation in the useful tspan
tspan_E1 = TU*linspace(data.et0,optimal_etf,points); 
tspan_V1 = TU*linspace(data.et0,optimal_etf,points);

tspan_E2 = TU*linspace(optimal_etf,1.03*optimal_etf,points); 
tspan_V2 = TU*linspace(optimal_etf,1.02*optimal_etf,points);

rr_E_prop1 = cspice_spkpos('Earth', tspan_E1,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_V_prop1 = cspice_spkpos('Venus', tspan_V1,'ECLIPJ2000', 'NONE', 'SSB');

rr_E_prop2 = cspice_spkpos('Earth', tspan_E2,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_V_prop2 = cspice_spkpos('Venus', tspan_V2,'ECLIPJ2000', 'NONE', 'SSB'); 

% 'Time optimal Earth-Venus transfer orbit (@Sun ECLIPTIC J2000)'
figure(1)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
xlim([-1.3 1.3])
ylim([-1.5 1.5])
view(2)
% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Venus orbit','','', 'Sun',...
    'Departure from Earth','Arrival at Venus', 'Location',...
    'bestoutside','FontSize',14) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)

% 'Time optimal Earth-Venus transfer orbit (@Sun ECLIPTIC J2000)'
figure(2)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
view(3)
% Plot settings
set(gca,'FontSize',12)
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)


% Computing the thrust pointing unit vector  
alpha = zeros(length(time),3); 

for i = 1:size(Y,1)
    alpha(i,:) = -Y(i,11:13)./norm(Y(i,11:13)); 
end

% 'Thrust pointing unit vector emphasized' 
figure(3)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
quiver3(Y(1:4:end,1), Y(1:4:end,2), Y(1:4:end,3), alpha(1:4:end,1), ...
    alpha(1:4:end,2), alpha(1:4:end,3), 'Color', '#A3A3A3', 'LineWidth', 0.7);
view(3)
% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Venus orbit','','','Sun','Departure from Earth', ...
        'Arrival at Venus','Thrust pointing vector','Location','bestoutside',...
        'FontSize',14) 
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
zlabel('$z$ [-]', 'Interpreter', 'latex', 'FontSize', 20)


% Hamiltonian computation 
hamilton = zeros(size(Y,1),1); 

for i = 1:size(Y,1)
    hamilton(i) = 1 + dot(Y(i,8:10),Y(i,4:6)) - (data.mu/norm(Y(i,1:3))^3)*dot(Y(i,11:13), ...
                Y(i,1:3)) + (data.Tmax/(data.Isp*data.g0))*(-norm(Y(i,11:13)) * ...
                data.Isp*data.g0/Y(i,7) - Y(i,14));
end

% 'Time evolution of the Hamiltonian along the time-optimal trajectory'
figure(4)
plot((time-data.et0)*TU/cspice_spd, hamilton,'r','LineWidth', 2); 
grid on ; axis equal 
% Plot settings
set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('H [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Hamiltonian','FontSize',14)


% Switching function computation 
St = zeros(length(time),1);

for i = 1 : length(time)
    St(i) = -norm(Y(i, 11:13))*data.Isp*data.g0/Y(i,7) - Y(i,14); 
end

% 'Time Evolution of the switching function'
figure(5)
plot((time-data.et0)*TU/cspice_spd, St,'r','LineWidth',2)
grid on ;  
% Plot settings
set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('$S_t$ [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Switching function','FontSize',14)
xlim([0,ToF]) 
ylim([-14,2]) 

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REQUEST 4

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 3 Request 4  \n \n')

% Data 
t0 = '2023-05-28-14:13:09.000 UTC'; 
et0 = cspice_str2et(t0); % Ephemeris time [s] 
xx0 = cspice_spkezr('Earth',et0,'ECLIPJ2000','NONE','Sun'); % Initial state vector [km,km/s]
rr0 = xx0(1:3); % Earth's initial position [km]                                                   
vv0 = xx0(4:6); % Earth's initial velocity [km/s]                                                      

m0 = 1000;  % Spacecraft mass [kg]
Tmax = 500e-6;  % Maximum thrust  [kN]
Isp = 3120;  % Specific impulse [s]
mu = cspice_bodvrd('Sun','GM',1);  % Sun gravitational parameter [km^3/s^2]
g0 = 9.80665e-3;  % Gravitational constant [km/s^2]

% Adimensionalizing unit
LU = cspice_convrt(1,'AU','km'); % Length unit (Astronomical unit) [km]
MU = m0; % Mass unit [kg]
TU = sqrt(LU^3/mu); % Time unit [s] 

% Adimensional parameters
data.mu = 1; % Adimensional Sun parameter[-]
data.rr0 = rr0/LU; % Adimensional initial position [-]
data.vv0 = vv0/LU*TU; % Adimensional initial velocity [-] 
data.Tmax = Tmax*TU^2/MU/LU; % Adimensional thrust [-]
data.Isp = Isp/TU; % Adimensional specific impulse [-]
data.g0 = g0*TU^2/LU; % Adimensional sea level acceleration constant [-]
data.m0 = m0/MU; % Adimensional mass [-]
data.et0 = et0/TU; % Adimensional initial time [-]

error = 1;
tol = 1e-9; % Set error tolerance for exit criteria

% Initial guess for the fsolve problem using numerical continuation
tf_guess = 4.669239;
lambda0 = [0.509268; -13.311608; 0.115802; 5.046448;...
    -10.018109; 1.456822; 1.709593];
etf = tf_guess + data.et0;

% Optimal lambda and tf finder
while error > tol

    %Finding the roots
    options = optimset('Algorithm', 'levenberg-marquardt', 'Display', ...
                        'iter','TolFun',1e-15, 'MaxIter', 2000, 'MaxFunEvals', 2000); 
    [Y_opt,fval] = fsolve(@(Y) shooting(Y,data), [lambda0; etf], options);

    % Check for the error
    error = norm(fval);

    lambda0 = Y_opt(1:7);
    etf = Y_opt(8);

end

optimal_lambda0 = lambda0;
optimal_etf = etf;

% Time of Flight [days]
ToF = (optimal_etf - data.et0)*TU/cspice_spd;

% Final state error with respect to Venus 
residual = fval;

err_rr = norm(residual(1:3))*LU ; % Final potision error [km]
err_vv = norm(residual(4:6))*LU/TU*1e3 ;  % Final velocity error [m/s]

% Final time string
Final_date = cspice_et2utc(optimal_etf*TU,'C',1); 

% % Display results
fprintf('Solutions of the time-optimal Earth-Venus transfer considering 500 mN continuous thrust are:\n\n')
fprintf('   λ = [');
fprintf('%f ',optimal_lambda0(1:7))
fprintf(']''\n')
fprintf('   Arrival date: %s UTC \n',Final_date)
fprintf('   Time of flight: %f days\n', ToF)
fprintf('   Position error: %f km\n', err_rr)
fprintf('   Velocity error: %e m/s\n', err_vv)
fprintf('------------------------------------------------------------------------------------\n')

% Integrating state and the costate with the optimal initial conditions 
options = odeset('AbsTol',1e-13 , 'RelTol', 1e-13); 
tspan = linspace(data.et0,optimal_etf,1000); 
Y0 = [data.rr0; data.vv0 ; data.m0; optimal_lambda0]; 
[time, Y] = ode113(@(t,Y) TPBVP(t,Y,data),tspan, Y0, options);

points = 100000;

% Propagation in the useful tspan
tspan_E1 = TU*linspace(data.et0,optimal_etf,points); 
tspan_V1 = TU*linspace(data.et0,optimal_etf,points);

tspan_E2 = TU*linspace(optimal_etf,1.03*optimal_etf,points); 
tspan_V2 = TU*linspace(optimal_etf,1.02*optimal_etf,points);

rr_E_prop1 = cspice_spkpos('Earth', tspan_E1,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_V_prop1 = cspice_spkpos('Venus', tspan_V1,'ECLIPJ2000', 'NONE', 'SSB');

rr_E_prop2 = cspice_spkpos('Earth', tspan_E2,'ECLIPJ2000', 'NONE', 'SSB'); 
rr_V_prop2 = cspice_spkpos('Venus', tspan_V2,'ECLIPJ2000', 'NONE', 'SSB'); 

% 'Time optimal Earth-Venus transfer orbit (@Sun ECLIPTIC J2000)'
figure(6)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
xlim([-1.3 1.3])
ylim([-1.5 1.5])
view(2)
% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Venus orbit','','', 'Sun',...
    'Departure from Earth','Arrival at Venus', 'FontSize',14,'Location',...
    'bestoutside','FontSize',12) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)

% 'Time optimal Earth-Venus transfer orbit (@Sun ECLIPTIC J2000)'
figure(7)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
view(3)
% Plot settings
set(gca,'FontSize',12)
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)


% Computing the thrust pointing unit vector  
alpha = zeros(length(time), 3);

for i = 1:size(Y,1)
    alpha(i,:) = -Y(i,11:13)./norm(Y(i,11:13)); 
end

% 'Thrust pointing unit vector emphasized'
figure(8)
plot3(Y(:,1), Y(:,2), Y(:,3), 'color', '#FF0000', 'LineWidth',2)
hold on ; grid on ; 
plot3(rr_E_prop1(1,:)/LU, rr_E_prop1(2,:)/LU, rr_E_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#0072BD')
plot3(rr_V_prop1(1,:)/LU, rr_V_prop1(2,:)/LU, rr_V_prop1(3,:)/LU, '-',...
    'LineWidth', 1.5, 'Color', '#EDB120') 
plot3(rr_E_prop2(1,:)/LU, rr_E_prop2(2,:)/LU, rr_E_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#0072BD')
plot3(rr_V_prop2(1,:)/LU, rr_V_prop2(2,:)/LU, rr_V_prop2(3,:)/LU, '--',...
    'LineWidth', 1, 'Color', '#EDB120') 
scatter3(0, 0, 0, 'filled', 'SizeData', 300,'MarkerFaceColor','#D95319') % Sun
scatter3(rr_E_prop1(1,1)/LU, rr_E_prop1(2,1)/LU, rr_E_prop1(3,1)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#0072BD') % Earth 
scatter3(rr_V_prop1(1, end)/LU, rr_V_prop1(2, end)/LU, rr_V_prop1(3, end)/LU,...
    'filled', 'SizeData', 140,'MarkerFaceColor','#EDB120') % Venus
quiver3(Y(1:4:end,1), Y(1:4:end,2), Y(1:4:end,3), alpha(1:4:end,1), ...
    alpha(1:4:end,2), alpha(1:4:end,3), 'Color', '#A3A3A3', 'LineWidth', 0.7);
xlim([-1.4,1.4])
ylim([-1.4,1.4])
view(3)
% Plot settings
set(gca,'FontSize',12)
legend('Transfer orbit','Earth orbit','Venus orbit','','','Sun','Departure from Earth', ...
        'Arrival at Venus','Thrust pointing vector','Location','bestoutside',...
        'FontSize',14) 
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 20)
zlabel('$z$ [-]', 'Interpreter', 'latex', 'FontSize', 20)


% Hamiltonian computation and plot
hamilton = zeros(size(Y,1),1); 

for i = 1:size(Y,1)
    hamilton(i) = 1 + dot(Y(i,8:10),Y(i,4:6)) -  (data.mu/norm(Y(i,1:3))^3)*dot(Y(i,11:13), ...
                Y(i,1:3)) + (data.Tmax/(data.Isp*data.g0))*(-norm(Y(i,11:13)) * ...
                data.Isp*data.g0/Y(i,7) - Y(i,14));
end

% 'Time evolution of the Hamiltonian along the time-optimal trajectory'
figure(9)
plot((time-data.et0)*TU/cspice_spd, hamilton,'r','LineWidth', 2) ; 
grid on ; axis equal 
% plot settings
set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('H [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Hamiltonian','FontSize',14)

% Switching function computation 
St = zeros(length(time),1); 

for i = 1 : length(time)
    St(i) = -norm(Y(i, 11:13))*data.Isp*data.g0/Y(i,7) - Y(i,14); 
end

% 'Time Evolution of the switching function'
figure(10)
plot((time-data.et0)*TU/cspice_spd, St,'r','LineWidth',2)
grid on ;  
% plot settings
set(gca,'FontSize',12)
xlabel('time [days]','FontSize',20, 'Interpreter', 'Latex')
ylabel('$S_t$ [-]','FontSize', 20, 'Interpreter', 'Latex')
legend('Switching function','FontSize',14)
xlim([0,ToF]) 
ylim([-30,4]) 

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS

function [dydt] = TPBVP(~,Y,data)
% The function provides the odefun for the two-point boundary value
% problem
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   Y : state and costate vector:  [14,1]             
%         -State:
%          r     Position     [3,1]    
%          v     Velocity     [3,1]    
%          m     Mass     [1,1]     
%         -Costate:
%          lambda_r     Position costate     [3,1]    
%          lambda_v     Velocity costate     [3,1]    
%          lambda_m     Mass costate     [1,1]   
%   data : data of the problem   [struct]    
%
% Outputs:
%   dydt : derivative of state and costate vector   [14,1]
%

% data extraction
mu = data.mu;
Tmax = data.Tmax; 
Isp= data.Isp; 
g0 = data.g0;

% Extract the state vectors
rr = Y(1:3);   
vv = Y(4:6);

% Extract the mass
m = Y(7);  

% Extract the costate vectors
lambda_rr = Y(8:10);
lambda_vv = Y(11:13);

% Compute the norms
r_norm = norm(rr); 
lambda_vv_norm = norm(lambda_vv);

% Time-optimal throttling factor
u_star = 1 ;   

% State and costate derivatives
dydt = [vv(1); vv(2); vv(3);
    -mu*rr(1)/r_norm^3 - (u_star*Tmax/m)*(lambda_vv(1)/lambda_vv_norm);...
    -mu*rr(2)/r_norm^3 - (u_star*Tmax/m)*(lambda_vv(2)/lambda_vv_norm);...
    -mu*rr(3)/r_norm^3 - (u_star*Tmax/m)*(lambda_vv(3)/lambda_vv_norm);...
    -u_star*Tmax/(Isp*g0);...
    -3*mu*dot(rr,lambda_vv)*rr/r_norm^5 + mu*lambda_vv/r_norm^3;...
    -lambda_rr;...
    -u_star*lambda_vv_norm*Tmax/(m^2)];

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = shooting(yy,data) 

% The function provides the shooting function for the time optimal 
% TPBVP for continous low-thrust guidance in the two-body problem.
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   yy : unknowns array:  [8,1]             
%          lambda_r     [3,1]    
%          lambda_v     [3,1]    
%          lambda_m     [1,1]     
%          tf     [1,1]       
%   data : data of the problem   [struct]    
%
% Outputs:
%   F : shooting function   [8,1]
%

% data extraction 
mu = data.mu;
Tmax = data.Tmax; 
Isp = data.Isp; 
et0 = data.et0; 
rr0 = data.rr0;
vv0 = data.vv0; 
g0 = data.g0;
m0 = data.m0;
u_star = 1; 

% Adimensional parameters 
LU = 149597870700e-3;
TU = 5.0226e6; 

%Propagation
y0 = [rr0; vv0; m0; yy(1:7)] ;
options = odeset('RelTol', 2.5e-13, 'AbsTol', 2.5e-13);
[~, YY] = ode113(@(t,Y) TPBVP(t,Y,data), [et0 yy(8)], y0, options);

% Retrieving the final conditions 
YY_f = YY(end,:)'; 

rr_f = YY_f(1:3); % position
vv_f = YY_f(4:6);  % velocity 
m_f = YY_f(7);  % mass

lambda_rr_f = YY_f(8:10); % lambda position                
lambda_vv_f = YY_f(11:13); % lambda velocity          
lambda_m_f = YY_f(14); % lambda mass   

% Venus' state vector at final time
state_venus = cspice_spkezr('Venus', yy(8)*TU, 'ECLIPJ2000', 'NONE', 'Sun');

% Pre-allocating for speed
rr_venus = zeros(3,1); % Position 
vv_venus = zeros(3,1); % Velocity 

% Extract final position, velocity and acceleration of Venus
rr_venus(1:3) = state_venus(1:3)/LU;  % Venus position vector
vv_venus(1:3) = state_venus(4:6)/LU*TU;  % Venus velocity vector
aa_venus = -(mu/norm(rr_venus)^3)*rr_venus;  % Venus acceleration vector

% Hamiltonian at final time 
H = 1 + dot(lambda_rr_f,vv_f) - mu/(norm(rr_f)^3)*dot(rr_f, lambda_vv_f)...
    - u_star*Tmax*norm(lambda_vv_f)/m_f - lambda_m_f*u_star*Tmax/(Isp*g0);

% zero-conditions vector
F = zeros(8,1);

F(1:3) = rr_f - rr_venus;
F(4:6) = vv_f - vv_venus;
F(7) = lambda_m_f;
F(8) = H - dot(lambda_rr_f,vv_venus) - dot(lambda_vv_f,aa_venus);

end



















