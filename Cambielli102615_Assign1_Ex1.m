%% Spacecraft Guidance and Navigation - Periodic Orbit (2023/2024)
% Assignment:     1
% Exercise:       1 - Period orbits
% Author:         Alessandro Cambielli

% ID number:      10619322
%                 102615

%% REQUEST 1

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 1 Request 1  \n \n')

% Initialization
mu = 0.012150; % Moon-Earth CR3BP

% L1 position
f = @(x) x-(1-mu)*(x+mu)./abs(x+mu).^3-mu*(x+mu-1)./abs(x+mu-1).^3;

options = optimset('Display','iter');
L_1 = fzero(f,[0.01, 0.97],options);

fprintf('L1=%.10f \n \n', L_1)

% 'The Lagrange L1 point in the Earth-Moon system'
figure(1)
grid on
hold on
fplot(f,'k')
plot(-mu,0,'mo','MarkerFaceColor','m')
plot(1-mu,0,'co','MarkerFaceColor','c')
plot(L_1,0,'r*','MarkerFaceColor','r')
axis([-2 2 -30 30])
% Plot settings
set(gca,'FontSize',12)
legend('dU/dx','m1','m2','L1', 'Location',...
    'northeast','FontSize',10) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$dU/dx$ [-]','Interpreter','latex','FontSize', 20)

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REQUEST 2

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 1 Request 2  \n \n')

% Initialization
mu = 0.012150; % Moon-Earth CR3BP
t0 = 0;
tf = 2; % Propagation time

%Initial conditions
x0 = [1.08892819445324; 0; 0.0591799623455459;...
       0; 0.257888699435051; 0];

options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14,'Events',@x_axis_crossing);
optionsSTM = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11],'Events',@x_axis_crossing);

% L2 position
f = @(x) x-(1-mu)*(x+mu)./abs(x+mu).^3-mu*(x+mu-1)./abs(x+mu-1).^3;
L_2 = fzero(f,[1, 1.5]);

[~, xx_initial] = ode113(@(t,X) CR3BP(t,X,mu), [t0 tf], x0, options);

% 'Initial L2 halo orbit in the Earth-Moon System'
figure(2)
grid on
hold on
plot3(xx_initial(:,1),xx_initial(:,2),xx_initial(:,3),'r','LineWidth',2)
plot(L_2,0,'ko','MarkerFaceColor','k')
xlim([1.07 1.17])
% Plot settings
set(gca,'FontSize',12)
legend('initial halo orbit', 'L2', 'Location', 'northeast','FontSize',14) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)

% Differential correction scheme

xx = ones(1,6);
counter = 0;

while abs(xx(end,4)) > 1e-8 && abs(xx(end,6)) > 1e-8

    %Integration
    [~, xx] = ode113(@(t,X) CR3BP_STM(t,X,mu), [t0 tf],...
        [x0; reshape(eye(6),[],1)], options);

    % Retrieve the final STM
    STMf = (reshape(xx(end,7:42),6,6))';

    %Correction
    deltas = [-xx(end,4); -xx(end,6)];
    el1 = [STMf(4,1) STMf(4,5); STMf(6,1) STMf(6,5)];
    corr = el1 \ deltas;
    x0(1) = x0(1) + corr(1);
    x0(5) = x0(5) + corr(2);

    counter = counter + 1;

end

fprintf('convergence is reached in %d iterations \n', counter)
fprintf('new initial conditions found \n \n')

% Final integration
tff = 5;  % Propagation time
options_per = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
[tt, xx_periodic] = ode113(@(t,X) CR3BP(t,X,mu), [t0 tff], x0, options_per);

% 'Periodic L2 halo orbit in the Earth-Moon System'
figure(3)
grid on
hold on
plot3(xx_periodic(:,1),xx_periodic(:,2),xx_periodic(:,3),'r','LineWidth',2)
plot(L_2,0,'ko','MarkerFaceColor','k')
xlim([1.07 1.19])
% Plot settings
set(gca,'FontSize',12)
legend('periodic halo orbit', 'L2', 'Location', 'northeast','FontSize',14) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 20)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 20)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 20)

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% REQUEST 3

% Clear memory workspace and set path
clearvars; close all; clc 

% Load SPICE kernels: 
cspice_furnsh('ex02.tm');
format long g

fprintf('Exercise 1 Request 3  \n \n')

% Initialization
mu = 0.012150; % Moon-Earth CR3BP
t0 = 0;
tf = 2; % Propagation time
tff = 5; % Final propagation time

options_per = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-14,'Events',@x_axis_crossing);

x0 = [1.09027817963494; 0; 0.0591799623455459;...
       0; 0.260348873992119; 0];

z0_in = 0.0591799623455459;
z0_fin = 0.034;

values = linspace(z0_in, z0_fin, 12);

cc = hsv(length(values));

for iter = 1:length(values)  

    xx = ones(1,6);

    %Initialization
    x0(3) = values(iter);

    while abs(xx(end,4)) > 1e-8 && abs(xx(end,6)) > 1e-8

        %Integration
        [~, xx] = ode113(@(t,X) CR3BP_STM(t,X,mu), [t0 tf],...
        [x0; reshape(eye(6),[],1)], options);

        STMf = (reshape(xx(end,7:42),6,6))';
    
        %Correction
        deltas = [-xx(end,4); -xx(end,6)];
        el1 = [STMf(4,1) STMf(4,5); STMf(6,1) STMf(6,5)];
        corr = el1 \ deltas;
        x0(1) = x0(1) + corr(1);
        x0(5) = x0(5) + corr(2);

    end

% L2 position
f = @(x) x-(1-mu)*(x+mu)./abs(x+mu).^3-mu*(x+mu-1)./abs(x+mu-1).^3;
L_2 = fzero(f,[1, 1.5]);

[tt, xx_periodic] = ode113(@(t,X) CR3BP(t,X,mu), [t0 tff], x0, options_per); 

% 'Family of L2 periodic halo orbits in the Earth-Moon system'
figure(4)
grid on
hold on
plot3(xx_periodic(:,1),xx_periodic(:,2),xx_periodic(:,3),'color',cc(iter,:),'LineWidth',2)
plot(L_2,0,'ko','MarkerFaceColor','k')
xlim([1.08 1.185])
xticks = (1.08:6:1.185);
yticks = (-0.15:6:0.15);
zticks = (-0.1:6:0.6);
view(2)
% Plot settings
set(gca,'FontSize',12)
legend('periodic halo orbits', 'L2', 'Location',...
    'northeast','FontSize',16) 
xlabel('$x$ [-]','Interpreter','latex','FontSize', 24)
ylabel('$y$ [-]','Interpreter','latex','FontSize', 24)
zlabel('$z$ [-]','Interpreter','latex','FontSize', 24)

end

% Clear kernels
cspice_kclear(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS


function [dxdt] = CR3BP(~, X, mu)
% The function provides the odefun for the circular restricted three-boy
% problem
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   X : cartesian state vector wrt Solar-System-Barycentre [6,1]
%   mu : gravitational parameter [1,1]    
%
% Outputs:
%   dxdt : [6,1] 
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract variables
x  = X(1);
y  = X(2);
z  = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

% Compute distances from bodies 1 and 2
r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

 % Compute derivative of the potential
dUdx = x-(1-mu)*(x+mu)/r1^3-mu*(x+mu-1)/r2^3;
dUdy = y-(1-mu)*y/r1^3-mu*y/r2^3;
dUdz = -(1-mu)*z/r1^3-mu*z/r2^3;

dxdt(1:3) = [vx;vy;vz]; % Position detivative is object's velocity
dxdt(4:6) = [2*vy + dUdx; -2*vx + dUdy; dUdz]; % Sum up acceleration to right-hand-side

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dstate = CR3BP_STM(~, state, mu)
% The function provides the vectorial field for the circular restricted
% three-body problem merging the state with the computation of the state
% transition matrix
%
% Author: Alessandro Cambielli
%  
% Inputs:
%   state : cartesian state and STM [42,1]
%   mu : gravitational parameter [1,1]    
%
% Outputs:
%   dstate : derivative of cartesian state and STM [42,1] 
%

%Initialize
dstate = zeros(42,1);

x = state(1);
y = state(2);
z = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

STM = reshape(state(7:42),6,6);   %From equations to STM

% CRTBP dynamics
r1 = sqrt((x+mu)^2+y^2+z^2);
r2 = sqrt((x+mu-1)^2+y^2+z^2);

% Compute the derivatives of the flow
df4dx = 1-(1-mu)/r1^3+3*(1-mu)*(x+mu)^2/r1^5-mu/r2^3+3*mu*(x+mu-1)^2/r2^5;
df4dy = 3*(1-mu)*(x+mu)*y/r1^5+3*mu*(x+mu-1)*y/r2^5;
df4dz = 3*(1-mu)*(x+mu)*z/r1^5+3*mu*(x+mu-1)*z/r2^5;
df5dy = 1-(1-mu)/r1^3+3*(1-mu)*y^2/r1^5-mu/r2^3+3*mu*y^2/r2^5;
df5dz = 3*(1-mu)*y*z/r1^5+3*mu*y*z/r2^5;
df6dz = -(1-mu)/r1^3+3*(1-mu)*z^2/r1^5-mu/r2^3+3*mu*z^2/r2^5;

% Assemble the matrix A(t) = dfdx
A = [0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 1;...
    df4dx, df4dy, df4dz, 0, 2, 0;...
    df4dy, df5dy, df5dz, -2, 0, 0;...
    df4dz, df5dz, df6dz, 0, 0, 0];

% STM derivative
dSTM = A*STM;

% assemble derivatives of cartesian state and STM
dstate(1:3) = [vx; vy; vz];
dstate(4:6) = [2*vy + x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;...
    -2*vx + y - (1-mu)*y/r1^3 - mu*y/r2^3;...
    - (1-mu)*z/r1^3 - mu*z/r2^3];
dstate(7:12) = dSTM(1,1:6)';
dstate(13:18) = dSTM(2,1:6)';
dstate(19:24) = dSTM(3,1:6)';
dstate(25:30) = dSTM(4,1:6)';
dstate(31:36) = dSTM(5,1:6)';
dstate(37:42) = dSTM(6,1:6)';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = x_axis_crossing(~,xx,~)
% x_axis_crossing is an event function that detects when the second
% component (y) of the state reaches 0.
%
% Author: Alessandro Cambielli
%
% Inputs:
%   xx : state vector [6,1]  
%
% Outputs:
%   value : expression describing the event [1,1]
%   isterminal : 1 if the integration is to terminate when the event occurs
%   direction : -1 locates only zeros where the event function is decreasing
%

value = xx(2);
isterminal = 1;
direction = -1;

end


