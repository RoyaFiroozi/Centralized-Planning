% Author: Roya Firoozi
% Year: 2020
% Motion Planning: Centralized Approach with 3 vehicles
yalmip('clear')
clear all
close all
%% Parameters

g = 9.81; % gravitational acceleration
delta_t = 0.1; % controller sampling time
d_min = 0.5; % desing parameter chosen by user, minimum desired distance between the vehicles
v_max = 15; % speed limit of the road
T_traj = 100; % simulation time
% road
noLane = 2; %number of lanes 
laneWidth = 3.7; % standard lane width 
road_right = 0;
road_center = road_right + laneWidth;
road_left = road_right + noLane * laneWidth;
 
% initial condition 
x1init = 7; y1init = laneWidth/2;
x2init = 0.5; y2init = road_center + laneWidth/2;
x3init = 14; y3init = road_center + laneWidth/2;

steer_limit = 0.3; % bound on steering angle, unit is radian
accel_limit = 4; % bound on acceleration, unit is m/s^2
change_steer_limit = 0.2; % bound on change of steering 
change_accel_limit = 0.5; % bound on change of acceleration 

% vehicle size: ioniq
l_f = 4.47/2; 
l_r = 4.47/2;
h = l_f + l_r; % length of vehicle
w = 1.82; % width of vehicle
%% Model and MPC Data and Variables
nz = 4; % Number of states: states are position (x,y), and heading psi, velocity, v.  z = [x,y,psi,v]
nu = 2; % Number of inputs: acceleration and steering are control inputs  u = [a,delta_f]
ny = 4;

% polytope dimensions
nlambda = 4;
nmu = 4;
ns = 2;
nA = 4;
nb = 4;


% MPC 
N = 7;  %MPC horizon
if N == 1
    u1 = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1)); 
    u1 = u1(1);
    u2 = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1)); 
    u2 = u2(1);
    u3 = sdpvar(repmat(nu,1,N+1),repmat(1,1,N+1)); 
    u3 = u3(1);
    beta1 = sdpvar(repmat(1,1,N+1+1), repmat(1,1,N+1+1));
    beta1 = beta1(1:2);
    beta2 = sdpvar(repmat(1,1,N+1+1), repmat(1,1,N+1+1));
    beta2 = beta2(1:2);
    beta3 = sdpvar(repmat(1,1,N+1+1), repmat(1,1,N+1+1));
    beta3 = beta3(1:2);
else
    u1 = sdpvar(repmat(nu,1,N),repmat(1,1,N));
    u1_prev = sdpvar(2,1);
    u2 = sdpvar(repmat(nu,1,N),repmat(1,1,N));
    u2_prev = sdpvar(2,1);
    u3 = sdpvar(repmat(nu,1,N),repmat(1,1,N));
    u3_prev = sdpvar(2,1);
    beta1 = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));
    beta2 = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));
    beta3 = sdpvar(repmat(1,1,N+1),repmat(1,1,N+1));

end
z1 = sdpvar(repmat(nz,1,N+1),repmat(1,1,N+1));
z2 = sdpvar(repmat(nz,1,N+1),repmat(1,1,N+1));
z3 = sdpvar(repmat(nz,1,N+1),repmat(1,1,N+1));
r1 = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1)); 
r2 = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));
r3 = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));



lambda12 = sdpvar(repmat(nlambda,1,N+1),repmat(1,1,N+1));
mu12 = sdpvar(repmat(nmu,1,N+1),repmat(1,1,N+1)); 
s12 = sdpvar(repmat(ns,1,N+1),repmat(1,1,N+1));

lambda13 = sdpvar(repmat(nlambda,1,N+1),repmat(1,1,N+1));
mu13 = sdpvar(repmat(nmu,1,N+1),repmat(1,1,N+1)); 
s13 = sdpvar(repmat(ns,1,N+1),repmat(1,1,N+1));

lambda23 = sdpvar(repmat(nlambda,1,N+1),repmat(1,1,N+1));
mu23 = sdpvar(repmat(nmu,1,N+1),repmat(1,1,N+1)); 
s23 = sdpvar(repmat(ns,1,N+1),repmat(1,1,N+1));


A1 = sdpvar(repmat(nA,1,N+1),repmat(2,1,N+1));
b1 = sdpvar(repmat(nb,1,N+1),repmat(1,1,N+1)); 

A2 = sdpvar(repmat(nA,1,N+1),repmat(2,1,N+1));
b2 = sdpvar(repmat(nb,1,N+1),repmat(1,1,N+1)); 

A3 = sdpvar(repmat(nA,1,N+1),repmat(2,1,N+1));
b3 = sdpvar(repmat(nb,1,N+1),repmat(1,1,N+1)); 

constraints = [];

objective = 0;

Q = 0.1*diag([1,100,1,1]);
R = 0.1*diag([1,1]);

%%
for k = 1:N
    objective = objective + (z1{k}-r1{k})'*Q*(z1{k}-r1{k}) + u1{k}'*R*u1{k} + ...
                            (z2{k}-r2{k})'*Q*(z2{k}-r2{k}) + u2{k}'*R*u2{k} + ...
                            (z3{k}-r3{k})'*Q*(z3{k}-r3{k}) + u3{k}'*R*u3{k};
    % Vehicle 1
    beta1{k} = atan((l_r./(l_f+l_r))*tan(u1{k}(2)));
    constraints = [constraints, z1{k+1}(1) == z1{k}(1) + delta_t*z1{k}(4)*cos(z1{k}(3) + beta1{k})]; %nonlinear MPC
    constraints = [constraints, z1{k+1}(2) == z1{k}(2) + delta_t*z1{k}(4)*sin(z1{k}(3) + beta1{k})]; % kinematic bicycle model
    constraints = [constraints, z1{k+1}(3) == z1{k}(3) + delta_t*(z1{k}(4)./l_r)*sin(beta1{k})];
    constraints = [constraints, z1{k+1}(4) == z1{k}(4) + delta_t*u1{k}(1)];
    constraints = [constraints,  z1{k}(2) + w/2 <= road_left]; % road boundary constraint 
    constraints = [constraints, -z1{k}(2) + w/2 <= -road_right]; % road boundary constraint
    constraints = [constraints, 0.0 <= z1{k}(4) <= v_max+4];

    constraints = [constraints, -accel_limit <= u1{k}(1) <= accel_limit]; % acceleration input constraint
    constraints = [constraints, -steer_limit <= u1{k}(2) <= steer_limit]; % steering input constraint 
    
    if k ~= 1
        constraints = [constraints, -change_accel_limit <= u1{k}(1)-u1{k-1}(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u1{k}(2)-u1{k-1}(2) <= change_steer_limit]; % change of steering input constraint 
    end
    if k == 1
        constraints = [constraints, -change_accel_limit <= u1{k}(1)-u1_prev(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u1{k}(2)-u1_prev(2) <= change_steer_limit]; % change of steering input constraint 
    end        
    
    % Vehicle 2
    beta2{k} = atan((l_r./(l_f+l_r))*tan(u2{k}(2)));
    constraints = [constraints, z2{k+1}(1) == z2{k}(1) + delta_t*z2{k}(4)*cos(z2{k}(3) + beta2{k})]; %nonlinear MPC
    constraints = [constraints, z2{k+1}(2) == z2{k}(2) + delta_t*z2{k}(4)*sin(z2{k}(3) + beta2{k})]; % kinematic bicycle model
    constraints = [constraints, z2{k+1}(3) == z2{k}(3) + delta_t*(z2{k}(4)./l_r)*sin(beta2{k})];
    constraints = [constraints, z2{k+1}(4) == z2{k}(4) + delta_t*u2{k}(1)];
    constraints = [constraints,  z2{k}(2) + w/2 <= road_left]; % road boundary constraint 
    constraints = [constraints, -z2{k}(2) + w/2 <= -road_right]; % road boundary constraint 
    constraints = [constraints, 0.0 <= z2{k}(4) <= v_max+4];
    
    constraints = [constraints, -accel_limit <= u2{k}(1) <= accel_limit]; % acceleration input constraint
    constraints = [constraints, -steer_limit <= u2{k}(2) <= steer_limit]; % steering input constraint 
    if k ~= 1
        constraints = [constraints, -change_accel_limit <= u2{k}(1)-u2{k-1}(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u2{k}(2)-u2{k-1}(2) <= change_steer_limit]; % change of steering input constraint 
    end
    if k == 1
        constraints = [constraints, -change_accel_limit <= u2{k}(1)-u2_prev(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u2{k}(2)-u2_prev(2) <= change_steer_limit]; % change of steering input constraint
    end
    
    % Vehicle 3
    beta3{k} = atan((l_r./(l_f+l_r))*tan(u3{k}(2)));
    constraints = [constraints, z3{k+1}(1) == z3{k}(1) + delta_t*z3{k}(4)*cos(z3{k}(3) + beta3{k})]; %nonlinear MPC
    constraints = [constraints, z3{k+1}(2) == z3{k}(2) + delta_t*z3{k}(4)*sin(z3{k}(3) + beta3{k})]; % kinematic bicycle model
    constraints = [constraints, z3{k+1}(3) == z3{k}(3) + delta_t*(z3{k}(4)./l_r)*sin(beta3{k})];
    constraints = [constraints, z3{k+1}(4) == z3{k}(4) + delta_t*u3{k}(1)];
    constraints = [constraints,  z3{k}(2) + w/2 <= road_left]; % road boundary constraint
    constraints = [constraints, -z3{k}(2) + w/2 <= -road_right]; % road boundary constraint
    constraints = [constraints, 0.0 <= z3{k}(4) <= v_max+4];
    
    constraints = [constraints, -accel_limit <= u3{k}(1) <= accel_limit]; % acceleration input constraint
    constraints = [constraints, -steer_limit <= u3{k}(2) <= steer_limit]; % steering input constraint
    if k ~= 1
        constraints = [constraints, -change_accel_limit <= u3{k}(1)-u3{k-1}(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u3{k}(2)-u3{k-1}(2) <= change_steer_limit]; % change of steering input constraint
    end
    if k == 1
        constraints = [constraints, -change_accel_limit <= u3{k}(1)-u3_prev(1) <= change_accel_limit]; % change of acceleration input constraint
        constraints = [constraints, -change_steer_limit <= u3{k}(2)-u3_prev(2) <= change_steer_limit]; % change of steering input constraint
    end
    
    % polytope constraints:
    [A1{k+1}, b1{k+1}] = rotation_translation([z1{k+1}(1); z1{k+1}(2)], z1{k+1}(3), h, w);
    [A2{k+1}, b2{k+1}] = rotation_translation([z2{k+1}(1); z2{k+1}(2)], z2{k+1}(3), h, w);
    [A3{k+1}, b3{k+1}] = rotation_translation([z3{k+1}(1); z3{k+1}(2)], z3{k+1}(3), h, w);
    
    % vehicle 1 w/ vehicle 2
    constraints = [constraints, b1{k}'*lambda12{k} + b2{k}'*mu12{k} <= -d_min];
    constraints = [constraints, A1{k}'*lambda12{k} + s12{k} == 0];
    constraints = [constraints, A2{k}'*mu12{k} - s12{k} == 0];
    constraints = [constraints, -lambda12{k} <= zeros(nlambda, 1)];
    constraints = [constraints, -mu12{k} <= zeros(nmu, 1)];
    constraints = [constraints, (s12{k}(1)^2 + s12{k}(2)^2)^(1/2) <= 1]; % norm one is used.
       
        % vehicle 1 w/ vehicle 3
    constraints = [constraints, b1{k}'*lambda13{k} + b3{k}'*mu13{k} <= -d_min];
    constraints = [constraints, A1{k}'*lambda13{k} + s13{k} == 0];
    constraints = [constraints, A3{k}'*mu13{k} - s13{k} == 0];
    constraints = [constraints, -lambda13{k} <= zeros(nlambda, 1)];
    constraints = [constraints, -mu13{k} <= zeros(nmu, 1)];
    constraints = [constraints, (s13{k}(1)^2 + s13{k}(2)^2)^(1/2) <= 1]; % norm one is used.
        
    % vehicle 2 w/ vehicle 3
    constraints = [constraints, b2{k}'*lambda23{k} + b3{k}'*mu23{k} <= -d_min];
    constraints = [constraints, A2{k}'*lambda23{k} + s23{k} == 0];
    constraints = [constraints, A3{k}'*mu23{k} - s23{k} == 0];
    constraints = [constraints, -lambda23{k} <= zeros(nlambda, 1)];
    constraints = [constraints, -mu23{k} <= zeros(nmu, 1)];
    constraints = [constraints, (s23{k}(1)^2 + s23{k}(2)^2)^(1/2) <= 1]; % norm one is used.
    
end

% Vehicle 1:
constraints = [constraints,  z1{N+1}(2) + w/2 <= road_left]; % road boundary constraint 
constraints = [constraints, -z1{N+1}(2) + w/2 <= -road_right]; % road boundary constraint
constraints = [constraints, 0.0 <= z1{N+1}(4) <= v_max+4];

% Vehicle 2:
constraints = [constraints,  z2{N+1}(2) + w/2 <= road_left]; % road boundary constraint 
constraints = [constraints, -z2{N+1}(2) + w/2 <= -road_right]; % road boundary constraint
constraints = [constraints, 0.0 <= z2{N+1}(4) <= v_max+4];

% Vehicle 3:
constraints = [constraints,  z3{N+1}(2) + w/2 <= road_left]; % road boundary constraint 
constraints = [constraints, -z3{N+1}(2) + w/2 <= -road_right]; % road boundary constraint
constraints = [constraints, 0.0 <= z3{N+1}(4) <= v_max+4];

% vehicle 1 w/ vehicle 2
constraints = [constraints, b1{N+1}'*lambda12{N+1} + b2{N+1}'*mu12{N+1} <= -d_min];
constraints = [constraints, A1{N+1}'*lambda12{N+1} + s12{N+1} == 0];
constraints = [constraints, A2{N+1}'*mu12{N+1} - s12{N+1} == 0];
constraints = [constraints, -lambda12{N+1} <= zeros(nlambda, 1)];
constraints = [constraints, -mu12{N+1} <= zeros(nmu, 1)];
constraints = [constraints, (s12{N+1}(1)^2 + s12{N+1}(2)^2)^(1/2) <= 1]; % norm one is used.

% vehicle 1 w/ vehicle 3
constraints = [constraints, b1{N+1}'*lambda13{N+1} + b3{N+1}'*mu13{N+1} <= -d_min];
constraints = [constraints, A1{N+1}'*lambda13{N+1} + s13{N+1} == 0];
constraints = [constraints, A3{N+1}'*mu13{N+1} - s13{N+1} == 0];
constraints = [constraints, -lambda13{N+1} <= zeros(nlambda, 1)];
constraints = [constraints, -mu13{N+1} <= zeros(nmu, 1)];
constraints = [constraints, (s13{N+1}(1)^2 + s13{N+1}(2)^2)^(1/2) <= 1]; % norm one is used.

% vehicle 2 w/ vehicle 3
constraints = [constraints, b2{N+1}'*lambda23{N+1} + b3{N+1}'*mu23{N+1} <= -d_min];
constraints = [constraints, A2{N+1}'*lambda23{N+1} + s23{N+1} == 0];
constraints = [constraints, A3{N+1}'*mu23{N+1} - s23{N+1} == 0];
constraints = [constraints, -lambda23{N+1} <= zeros(nlambda, 1)];
constraints = [constraints, -mu23{N+1} <= zeros(nmu, 1)];
constraints = [constraints, (s23{N+1}(1)^2 + s23{N+1}(2)^2)^(1/2) <= 1]; % norm one is used.

objective = objective + (z1{N+1}-r1{N+1})'*Q*(z1{N+1}-r1{N+1}) + ...
                        (z2{N+1}-r2{N+1})'*Q*(z2{N+1}-r2{N+1}) + ...
                        (z3{N+1}-r3{N+1})'*Q*(z3{N+1}-r3{N+1});

parameters_in = {z1{1}, z2{1}, z3{1}, [r1{:}], [r2{:}], [r3{:}], ... 
                 A1{1}, A2{1}, A3{1}, b1{1}, b2{1}, b3{1}, ...
                 u1_prev, u2_prev, u3_prev};
solutions_out = {[u1{:}], [u2{:}], [u3{:}], ...
                 [z1{:}], [z2{:}], [z3{:}], ...
                 [lambda12{:}], [mu12{:}], [s12{:}], ...
                 [lambda13{:}], [mu13{:}], [s13{:}], ...
                 [lambda23{:}], [mu23{:}], [s23{:}]};
% controller with 4 states, 2 inputs
controller = optimizer(constraints, objective, sdpsettings('solver','ipopt','verbos',2), parameters_in, solutions_out);

U1(:,1) = [0;0];
U2(:,1) = [0;0];
U3(:,1) = [0;0];
%% Generating Reference Trajectory

x_des1 = x1init; y_des1 = y1init; psi_des1 = 0; v_des1 = v_max;
ref1(:,1) = [x_des1 y_des1 psi_des1 v_des1]';

x_des2 = x2init; y_des2 = y2init; psi_des2 = 0; v_des2 = v_max;
ref2(:,1) = [x_des2 y_des2 psi_des2 v_des2]';

x_des3 = x3init; y_des3 = y3init; psi_des3 = 0; v_des3 = v_max;
ref3(:,1) = [x_des3 y_des3 psi_des3 v_des3]';

for i = 1:T_traj/4
    x_des_old1 = x_des1;
    y_des_old1 = y_des1;
    x_des1 = v_max*delta_t + x_des1;
    xdot1 = (x_des1 - x_des_old1)/delta_t;
    ydot1 = (y_des1 - y_des_old1)/delta_t;
    v_des1 = sqrt(xdot1^2+ydot1^2);
    ref1(:,i+1) = [x_des1 y_des1 psi_des1 v_des1]';
    
    x_des_old2 = x_des2;
    y_des_old2 = y_des2;
    x_des2 = v_max*delta_t + x_des2;
    xdot2 = (x_des2 - x_des_old2)/delta_t;
    ydot2 = (y_des2 - y_des_old2)/delta_t;
    v_des2 = sqrt(xdot2^2+ydot2^2);
    ref2(:,i+1) = [x_des2 y_des2 psi_des2 v_des2+0.5]';
    
    x_des_old3 = x_des3;
    y_des_old3 = y_des3;
    x_des3 = v_max*delta_t + x_des3;
    xdot3 = (x_des3 - x_des_old3)/delta_t;
    ydot3 = (y_des3 - y_des_old3)/delta_t;
    v_des3 = sqrt(xdot3^2+ydot3^2);
    ref3(:,i+1) = [x_des3 y_des3 psi_des3 v_des3]';
end

for i = ((T_traj/4)+1):(T_traj/2)
    x_des_old1 = x_des1;
    y_des_old1 = y_des1;
    y_des1 = y_des1 + laneWidth/((T_traj/2 - ((T_traj/4))));
    x_des1 = v_max*delta_t + x_des1;
    xdot1 = (x_des1 - x_des_old1)/delta_t;
    ydot1 = (y_des1 - y_des_old1)/delta_t;
    v_des1 = sqrt(xdot1^2+ydot1^2);
    ref1(:,i+1) = [x_des1 y_des1 psi_des1 v_des1]';
    
    x_des_old2 = x_des2;
    y_des_old2 = y_des2;
    x_des2 = v_max*delta_t + x_des2;
    xdot2 = (x_des2 - x_des_old2)/delta_t;
    ydot2 = (y_des2 - y_des_old2)/delta_t;
    v_des2 = sqrt(xdot2^2+ydot2^2);
    ref2(:,i+1) = [x_des2 y_des2 psi_des2 v_des2+0.5]';
    
    x_des_old3 = x_des3;
    y_des_old3 = y_des3;
    x_des3 = v_max*delta_t + x_des3;
    xdot3 = (x_des3 - x_des_old3)/delta_t;
    ydot3 = (y_des3 - y_des_old3)/delta_t;
    v_des3 = sqrt(xdot3^2+ydot3^2);
    ref3(:,i+1) = [x_des3 y_des3 psi_des3 v_des3]';
end

for i = ((T_traj/2)+1):T_traj
    x_des_old1 = x_des1;
    y_des_old1 = y_des1;
    x_des1 = v_max*delta_t + x_des1;
    xdot1 = (x_des1 - x_des_old1)/delta_t;
    ydot1 = (y_des1 - y_des_old1)/delta_t;
    v_des1 = sqrt(xdot1^2+ydot1^2);
    ref1(:,i+1) = [x_des1 y_des2 psi_des1 v_des1]';
    
    x_des_old2 = x_des2;
    y_des_old2 = y_des2;
    x_des2 = v_max*delta_t + x_des2;
    xdot2 = (x_des2 - x_des_old2)/delta_t;
    ydot2 = (y_des2 - y_des_old2)/delta_t;
    v_des2 = sqrt(xdot2^2+ydot2^2);
    ref2(:,i+1) = [x_des2 y_des2 psi_des2 v_des2]';
    
    x_des_old3 = x_des3;
    y_des_old3 = y_des3;
    x_des3 = v_max*delta_t + x_des3;
    xdot3 = (x_des3 - x_des_old3)/delta_t;
    ydot3 = (y_des3 - y_des_old3)/delta_t;
    v_des3 = sqrt(xdot3^2+ydot3^2);
    ref3(:,i+1) = [x_des3 y_des3 psi_des3 v_des3]';
end

T = length(ref1);

%% Simulation

% Initial Conditions
x_traj1 = zeros(nz,T);
u_traj1 = zeros(nu,T);
z1 = ref1(:,1);
x_traj1(:,1) = z1;

x_traj2 = zeros(nz,T);
u_traj2 = zeros(nu,T);
z2 = ref2(:,1);
x_traj2(:,1) = z2;

x_traj3 = zeros(nz,T);
u_traj3 = zeros(nu,T);
z3 = ref3(:,1);
x_traj3(:,1) = z3;

a1 = 0; delta_f1 = 0;
a2 = 0; delta_f2 = 0;
a3 = 0; delta_f3 = 0;

figure
hold on
scatter(z1(1), z1(2), 'or');
scatter(z2(1), z2(2), 'ob');
scatter(z3(1), z3(2), 'og');
axis equal
xlim([0.0, 100])
ylim([-0.5, noLane*road_center+0.5])
pause(0.1)

for i = 1:T-N
    [A1_0, b1_0] = rotation_translation([z1(1); z1(2)], z1(3), h, w);
    [A2_0, b2_0] = rotation_translation([z2(1); z2(2)], z2(3), h, w);
    [A3_0, b3_0] = rotation_translation([z3(1); z3(2)], z3(3), h, w);
    inputs = {z1, z2, z3, ref1(:,i:i+N), ref2(:,i:i+N), ref3(:,i:i+N), A1_0, A2_0, A3_0, b1_0, b2_0, b3_0, [a1; delta_f1], [a2; delta_f2], [a3; delta_f3]};
    [solutions, diagnostics,~,~,~,diag] = controller{inputs};
    timer(i) = diag.solvertime;
    U1 = solutions{1};
    U2 = solutions{2};
    U3 = solutions{3};
    X1 = solutions{4};
    X2 = solutions{5};
    X3 = solutions{6};
    lambda_12(:,:,i) = solutions{7};
    lambda_21(:,:,i) = solutions{8};
    s_12(:,:,i) = solutions{9};
    lambda_13(:,:,i) = solutions{10};
    lambda_31(:,:,i) = solutions{11};
    s_13(:,:,i) = solutions{12};
    lambda_23(:,:,i) = solutions{13};
    lambda_32(:,:,i) = solutions{14};
    s_23(:,:,i) = solutions{15};
    if diagnostics ~= 0
        error('The problem is infeasible');
    end
    a1 = U1(1,1);
    delta_f1 = U1(2,1); % steering angle
    x1 = z1(1); y1 = z1(2); psi1 = z1(3); v1 = z1(4);
    [x1,y1,psi1,v1] = kinematic_bicycle_model(x1, y1, psi1, v1, delta_t, l_f, l_r, a1, delta_f1);
    z1 = [x1,y1,psi1,v1]';
    x_traj1(:,i+1) = z1;
    u_traj1(:,i) = U1(:,1);
    
    a2 = U2(1,1); 
    delta_f2 = U2(2,1); % steering angle
    x2 = z2(1); y2 = z2(2); psi2 = z2(3); v2 = z2(4);
    [x2,y2,psi2,v2] = kinematic_bicycle_model(x2, y2, psi2, v2, delta_t, l_f, l_r, a2, delta_f2);
    z2 = [x2,y2,psi2,v2]';
    x_traj2(:,i+1) = z2;
    u_traj2(:,i) = U2(:,1);  
    
    a3 = U3(1,1); 
    delta_f3 = U3(2,1); % steering angle
    x3 = z3(1); y3 = z3(2); psi3 = z3(3); v3 = z3(4);
    [x3,y3,psi3,v3] = kinematic_bicycle_model(x3, y3, psi3, v3, delta_t, l_f, l_r, a3, delta_f3);
    z3 = [x3,y3,psi3,v3]';
    x_traj3(:,i+1) = z3;
    u_traj3(:,i) = U3(:,1);
    
    for j = 1:N+1
        [A_1(:,:,j), b_1(:,:,j)] = rotation_translation([X1(1,j); X1(2,j)], X1(3,j), h, w);
        [A_2(:,:,j), b_2(:,:,j)] = rotation_translation([X2(1,j); X2(2,j)], X2(3,j), h, w);
        [A_3(:,:,j), b_3(:,:,j)] = rotation_translation([X3(1,j); X3(2,j)], X3(3,j), h, w);
    end
    
    if i == 1
        A_1(:,:,1) = A1_0;
        A_2(:,:,1) = A2_0;
        A_3(:,:,1) = A3_0;
        b_1(:,:,1) = b1_0;
        b_2(:,:,1) = b2_0;
        b_3(:,:,1) = b3_0;
    end    
    diagnos(i) = diagnostics;
    
    % plot:
    scatter(x1, y1, 'or');
    scatter(x2, y2, 'ob');
    scatter(x3, y3, 'og');
    drawnow
    pause(0.1)
end


figure
hold on
plot(x_traj1(4,1:end-N))
plot(x_traj2(4,1:end-N))
plot(x_traj3(4,1:end-N),'--')
title('Velocity')
figure
hold on
plot(u_traj1(1,1:end-N))
plot(u_traj2(1,1:end-N))
plot(u_traj3(1,1:end-N),'--')
title('Acceleration')
disp(mean(timer))