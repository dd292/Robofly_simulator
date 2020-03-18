function [t, q, p, Q0, qd, LOG] = robofly_3dmatlab_simulator()
% 6DOF robofly simulator and renderer 
% Sawyer B. Fuller, 2014.07.28 
% 
% Simulates robofly dynamics including stroke-averaged aerodynamic drag.
% The simulation supports plugging in matlab transfer functions as
% regulators. Includes a few examples. The controllers are updated at 500
% Hz and then the dynamics are integrated with a fix-step integrator at 
% at 10 kHz. A variable-step integrator can also be enabled for more
% precision but is slower.
%
% The states that are computed are q = [theta, omegabody, posworld,
% vbody], that is, a 6DOF state vector with 12 components. The
% position posworld is given in world coordinates and the angles theta 
% are in XYZ Euler Angle convention (see function rot_matrix for details). 
% Velocity vectors (xdot_body and omega) are
% calculated in the coordinates of the body frame (but note that, as is
% convention, computed velocities are not _relative_ to the moving body
% frame, just given with respect the current orientation of its coordinate
% frame). To compute velocities in the world frame, use vworld = R *
% v_body, where R is the rotation matrix from the rot_matrix function.)
% This simulation loosely follows the treatment given in RB McGhee,
% ER Bachmann, MJ Zyda, "Rigid body dynamics, inertial reference frames,
% and graphics coordinate systems: A resolution of conflicting conventions
% and terminology," MOVES Academic Group Technical Report NPS-MV-01-002,
% Naval Postgraduate School, Monterey, CA (2000).
%
% useful simulation parameters: 
% * TFINAL: time duration of simulation  
% * ANIMATE: whether to animate while computing (small time cost)
% * PLOT: whether to plot trajectory results
% * RECALC: set to zero to plot previous results (eg. for a long simulation)
% * set the controller under 'controllers' comment. There 
%   a separate controller for each of the the altitude and torque dynamnics

% clear persistent vars in functions, otherwise problems with local state
clear all
%saved_db = dbstatus;
%clear <insert function name>
%dbstop(saved_db)

global CTRL_DT LOG Q0

RECALC = 1; % set to zero to plot previous results (eg. for a long simulation)
TFINAL = 4; % time duration of simulation
ANIMATE = 1; % whether to animate while computing (small time cost)
PLOT = 1; % whether to plot trajectory results
VARIABLE_STEP = 0; 
if ANIMATE
    fig = figure(22); clf()
    set(fig,'Name','Robofly simulator', ...
        'NumberTitle','off', ...
        'Position',[10,950,900,800], ...
        'Color',[0.1 0.1 0.2]);
         %'Menubar','none', ...
        %   %'Position',[10,10,1920, 1080], ... %1024,768], ...
    ax = axes('Position',[0 0 1 1]);
    axis equal
    light('Position',[-1 1 1]);
    light('Position',[-1 -1 1]);
    alpha(.8)
    view(30,30);
    set(ax,'Visible','off')
end
CTRL_HZ = 10000; %  update rate

% robofly params
p.winglength = .6 * 2.54 / 100; 
p.l = .013; % length of fly
p.h = .0025; % thickness of body
p.J = 1.5e-9; % 1d version
p.Jmat = diag([p.J, p.J, .5e-9]); %1.4e-9; % from science paper/solidworks
p.b_w = 2.0e-4; % aero drag on wings from wind tunnel tests, Ns/m
p.c_w = (p.h/2 + p.winglength * 2./3)^2 * p.b_w; % rot drag coeff around z [Nsm]
p.r_w = .007; %p.l; % z-dist from ctr of wings to ctr of mass 
p.m = 111e-6; %80e-6; %mass of fly
p.g = 9.81; 
p.max_f_l = 1.5 * p.m * p.g; 
p.max_torque = 2e-6;  
p.ks = 0;%0.5e-6; % tether stiffness
p.force_bias_x = 0; %0.0001; %N
p.torque_bias_y = 0; %-0.1e-6; %Nm
p.gyro_bias_y = 0; %0.1; % rad/s
% leave the following three zero for this simulation:
p.force_bias_y = 0; %N
p.torque_bias_x = 0; %Nm
p.gyro_bias_x = 0; % rad/s

% rendering
diagnl = 1/2*sqrt(p.l^2+p.h^2);
angl = atan(p.l/p.h);
plotbox = inline('plot([x+diag*cos(th-angl) x+diag*cos(th+angl) x+diag*cos(th-angl+pi) x+diag*cos(th+angl+pi) x+diag*cos(th-angl)], [y+diag*sin(th-angl) y+diag*sin(th+angl) y+diag*sin(th-angl+pi) y+diag*sin(th+angl+pi) y+diag*sin(th-angl)])','x','y','th','angl','diag');

CTRL_DT = 1/CTRL_HZ; 
t = 0:CTRL_DT:TFINAL; 

if RECALC
    %state: theta, omegabody, posworld, vbody
    Q0 = [.2,-.2,0,  -1,0,1, .04,.04,.01,  .1,-.3,0]';
    %Q0 = [0,0,0,  0,0,0, 0,0,.07,  0,0,0]';
    qd = [0,0,0,  0,0,0, 0,0,.08,  0,0,0]'; % desired state
        
    % initialize state
    LOG.q = zeros(length(t), length(Q0)); 
    LOG.u = zeros(length(t), 4);
    qout = Q0; 
 
    % controllers
    tauc_controller_step = @damping_controller_step;
    %tauc_controller_step = @open_loop_phase_shifter;
    f_l_controller_step = @altitude_controller_step;
    
    for idx = 1:length(t)
        tau_c = tauc_controller_step(t(idx), qout, qd, p);
        f_l = f_l_controller_step(t(idx), qout, qd, p); 
        %[f_l, tau_c] = wing_force_step(t(idx), qout, qd, p); 
        u = [f_l; tau_c]; 
        if VARIABLE_STEP
            [tode,qode] = ode45(@fly_dynamics, [t(idx), t(idx)+CTRL_DT], ...
                qout', [], u, p); 
            if size(qode, 1) > 100
                disp('unstable simulation detected, aborting');
                break
            end
            qout = qode(end,:)';
        else
            qout = qout + CTRL_DT * fly_dynamics(t(idx), qout, u, p); 
        end
            
        LOG.q(idx,:) = qout; 
        LOG.u(idx,:) = u'; 
        
        if ANIMATE
            if mod(idx, CTRL_HZ/300) == 0
                if 0 % 2d
                    plot([-.05, .05], [0, 0])
                    hold on
                    plotbox(LOG.q(idx, 8), LOG.q(idx, 9), LOG.q(idx, 1), angl, diagnl);
                    axis equal
                    plot(qd(8), qd(9), '+')
                    hold off
                    title(['t=' num2str(t(idx),'%1.3f')])
                    xlim([-.06, .06])
                else
                    render_body(LOG.q(idx,1:3)', LOG.q(idx,7:9)')
                    %title(['t=' num2str(t(idx),'%1.3f')], 'Color', [1,1,1])
                    h = text(0,0,.1, sprintf('t=%1.3f',t(idx)), 'Color', [1,1,1]);
                end
                drawnow
                delete(h)
            end
        end
    end
    save simulation_results.mat
else
    load simulation_results.mat
end

if PLOT
    fig = 1; iax = 1; 
    if 1
        figure(fig); fig = fig + 1; 
        % floor
        plot([-.05, .05], [0, 0])
        hold on
        %    trajectory of cm
        plot(LOG.q(:,8), LOG.q(:,9), 'k')
        for idx = 1:10:length(t)
            plotbox(LOG.q(idx, 8), LOG.q(idx, 9), LOG.q(idx, 1), angl, diagnl);
        end
        axis equal
        % setpoint location
        plot(qd(8), qd(9), '+')
        xlim([-.06, .06])
        xlabel('y position (m)')
        ylabel('z position (m)')
        hold off
    end 
    
    figure(fig); fig = fig + 1; 
    idx = 1;
    ax(iax) = subplot(4,3,idx); iax = iax+1;
    plot(t, LOG.q(:,0+idx), 'r'); ylim([-1, 1])
    ylabel(['theta', num2str(idx), ' (rad)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,idx); iax = iax+1;
    plot(t, LOG.q(:,0+idx), 'r'); ylim([-1, 1])
    ylabel(['theta', num2str(idx), ' (rad)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,idx); iax = iax+1;
    plot(t, LOG.q(:,0+idx), 'r'); ylim([-1, 1])
    ylabel(['theta', num2str(idx), ' (rad)'])
    
    idx = 1; 
    offset = 3; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:, offset+idx), 'r'); ylim([-5, 5])
    ylabel(['omega', num2str(idx), ' (rad/s)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:, offset+idx), 'r'); ylim([-5, 5])
    ylabel(['omega', num2str(idx), ' (rad/s)'])
    idx = idx + 1;
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:, offset+idx), 'r'); ylim([-5, 5])
    ylabel(['omega', num2str(idx), ' (rad/s)'])
    
    idx = 1;
    offset = 6; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.1, .1])
    ylabel(['pos', num2str(idx), ' (m)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.1, .1])
    ylabel(['pos', num2str(idx), ' (m)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.1, .1])
    ylabel(['pos', num2str(idx), ' (m)'])
    
    idx = 1;
    offset = 9; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.5, .5])
    ylabel(['v', num2str(idx), ' (m/s)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.5, .5])
    ylabel(['v', num2str(idx), ' (m/s)'])
    idx = idx + 1; 
    ax(iax) = subplot(4,3,offset+idx); iax = iax+1;
    plot(t, LOG.q(:,offset+idx), 'r'); ylim([-.5, .5])
    ylabel(['v', num2str(idx), ' (m/s)'])
        
    figure(fig); fig = fig + 1;  
    ax(iax) = subplot(211); iax = iax+1;
    plot(t, LOG.u(:,1))
    ylabel('f_l (N)')
    ax(iax) = subplot(212); iax = iax+1;
    plot(t, LOG.u(:,2:4))
    ylabel('torques (Nm)')
% 
%     figure(2)
%     subplot(221)
%     plot(t, LOG.q(:,4))
%     xlabel('time (s)')
%     ylabel('z (m)')
%     subplot(223)
%     plot(t, LOG.u(:, 1))
%     xlabel('time (s)')
%     ylabel('fl (N)')
%     subplot(222)
%     plot(t, LOG.q(:,1), 'r')
%     hold on 
%     plot(t, LOG.q(:,3))
%     hold off
%     xlabel('time (s)')
%     ylabel('theta (rad, red), x (m, blue)')
%     subplot(224)
%     plot(t, LOG.u(:, 2))
%     xlabel('time (s)')
%     ylabel('tauc (Nm)')
%     
    % bias corrections
    if isfield(LOG, 'correction_torque')
        figure(3)
        subplot(311)
        plot(t, LOG.correction_torque); 
        ylabel('correction torque (Nm)')
        subplot(312)
        plot(t, LOG.correction_theta); 
        ylabel('theta correction (rad)')
        subplot(313)
        plot(t, LOG.thetad); 
        ylabel('theta desired (rad)')
        xlabel('time (s)')
    end
    linkaxes(ax, 'x')
    xlim([0, TFINAL])
end

end

function [f_l, tau_c] = wing_force_step(t, q, ~, p)
% simulation to examine if phase offset can produce yaw rotation
flapping_f = 120; % rad/s
flapping_w = flapping_f * 2 * pi; 
right_wing_phase_offset = 60 * pi / 180; % rad
lcpwing = .01; % length to center of pressure of wing; assume 2/3 lwing. 
fl_amplitude_wing = p.m * p.g / 2; % assume lift goes as a sinusoid at 2f
fd_amplitude = fl_amplitude_wing; % assume drag = lift
if 0
    left_wing_drag_force = fd_amplitude * sin(flapping_w * t + 0);
    right_wing_drag_force = fd_amplitude * sin(flapping_w * t + right_wing_phase_offset);

    left_wing_yaw_torque =  lcpwing * left_wing_drag_force; 
    right_wing_yaw_torque = -lcpwing * right_wing_drag_force;
    yaw_torque = left_wing_yaw_torque + right_wing_yaw_torque;

    left_wing_pitch_torque =  -p.r_w * left_wing_drag_force; 
    right_wing_pitch_torque = -p.r_w * right_wing_drag_force;
    pitch_torque = left_wing_pitch_torque + right_wing_pitch_torque;


    left_wing_lift_force = fl_amplitude_wing * ...
        (1 +  sin(2*flapping_w * t + 0));
    right_wing_lift_force = fl_amplitude_wing * ...
        (1 + sin(2*flapping_w * t + right_wing_phase_offset));
    f_l = left_wing_lift_force + right_wing_lift_force; 

    left_wing_roll_torque = left_wing_lift_force * lcpwing; 
    right_wing_roll_torque = right_wing_lift_force * (-lcpwing); 
    roll_torque = left_wing_roll_torque + right_wing_roll_torque;
else
    f_l = p.m * p.g; 
    pitch_torque = fd_amplitude * lcpwing * sin(flapping_w * t); 
    %roll_torque = fd_amplitude * lcpwing * cos(flapping_w * t); 
    roll_torque = fd_amplitude * lcpwing * sin(2*flapping_w * t + -right_wing_phase_offset); 
    yaw_torque = fd_amplitude * lcpwing * sin(2*flapping_w * t + right_wing_phase_offset); 
end

% add stabilizing term
omegabody = q(4:6); 
tau_c = -1e-7 * omegabody; 
tau_c = clip(tau_c, p.max_torque); 
tau_c = tau_c + [roll_torque; pitch_torque; yaw_torque];
% tau_c = tau_c + [roll_torque; 0; yaw_torque];
% tau_c = [roll_torque; 0; yaw_torque];

end

function tau_c = lateral_nestedloop_controller_step(~, q, qd, p)
% inner-outer loop config based on lqr. ignores qd for theta and thetadot 
persistent initialized Ktheta Kx Kitheta Kix integral_of_ex1 integral_of_ex2 integral_of_etheta1 integral_of_etheta2 idx
global CTRL_DT LOG 
if isempty(initialized)
    Kx = 2*[4, .6]; % xworld, vworld
    Kix = 0; % integral component 
    Ktheta = 2*[2e-6, 1e-7]; % theta, omegabody
    Kitheta = 2e-6; % integral component
    integral_of_ex1 = 0; %-.22 / Kix; %0; 
    integral_of_etheta1 = 0; %2.99e-7 / Kitheta; %0;
    integral_of_ex2 = 0; 
    integral_of_etheta2 = 0; 
    idx = 1; 
    initialized = 1; 
end
% plant states: theta, omegabody, posworld, vbody    
ex1 = [qd(3) - q(3), qd(5) - q(5)]'; % x, v
ex1 = [qd(3) - q(3), qd(5) - q(5)]'; % x, v
integral_of_ex = integral_of_ex + CTRL_DT * (qd(3) - q(3)); 

%etheta = [qd(1) - q(1), qd(2) - q(2)]';% theta, thetadot
thetad = Kx * ex + Kix * integral_of_ex; 
thetad = min(1, max(-1, thetad));  
theta_error = 0;%.2; 
etheta = [thetad - q(1) + theta_error, -q(2)]';% theta, thetadot
integral_of_etheta = integral_of_etheta + CTRL_DT * (thetad - qd(1) + theta_error); 
disturbance = 0;%-.3e-6; 
tau_c = Ktheta * etheta + Kitheta * integral_of_etheta + disturbance;
LOG.correction_torque(idx) = Kitheta * integral_of_etheta;
LOG.correction_theta(idx) = Kix * integral_of_ex; 
LOG.thetad(idx) = thetad; 
idx = idx + 1; 
end


function tau_c = damping_controller_step(~, q, ~, p)
omegabody = q(4:6); 
tau_c = -2e-7 * omegabody; 
tau_c = clip(tau_c, p.max_torque); 
tau_c(3) = 0; 
end

function f_l = altitude_controller_step(~, q, qd, p)
% hand-tuned altitude controller
global CTRL_DT
persistent initialized A B C D x 
e = qd(9) - q(9);
if isempty(initialized)
    %controller_tf = tf([.0003, .001],[.0002, 1]) + tf(.5e-3, [1, .1]);
    controller_tf = tf([.01, .01, .02], [.001, 1, 0]); 
    [A, B, C, D] = ssdata(c2d(controller_tf, CTRL_DT, 'foh')); %first order hold
    x = zeros(size(A, 1), 1); 
    f_l = 0; % initial output  
    initialized = 1; 
else
    x = A * x + B * e; % update dynamics by dt
    f_l = C * x + D * e; 
end
f_l = f_l + .8*p.m * p.g; % feedforward gravity cancellation
f_l = clip(f_l, p.max_f_l); 
end

function f_l = constant_force_step(t, q, qd, p)
% for testing
if t < .1
    f_l = p.m * p.g;
else
    f_l = 0; 
end
f_l = f_l + p.m * p.g; % feedforward gravity cancellation
end


function clipped = clip(x, absmax)
    clipped = max(min(x, absmax), -absmax);
end

function A = estimate_state_jacobian(dynamics, t, q, varargin)
n_states = length(q); 
A = zeros(n_states); 
dq = sqrt(eps) * q; % optimal dq from wikipedia on numerical diff. 
for idx = 1:n_states
    dqvec = zeros(n_states, 1); 
    if q(idx) == 0
        dqvec(idx) = 1e-15; 
    else
        dqvec(idx) = dq(idx); 
    end
    A(:, idx) = (dynamics(t, q+dqvec, varargin{:}) - dynamics(t, q-dqvec, varargin{:}))/...
        (2 * dqvec(idx));
end
end

function qdot = fly_dynamics_wrapper(t, q, u, p)
% break out parts
torque_bias_estimate = q(13:15); 
qdot = zeros(18, 1); % assume no bias changes, update the rest

qdot(1:12) = fly_dynamics(t, q(1:12), u, p, torque_bias_estimate); 
% qdot(13:18) = 0; % assume no bias changes
end

function y = measurement_model(~, q, p)
gyro_bias = q(16:18); 
y = [q(4:6) + gyro_bias; q(10:12)];
end

function qdot = fly_dynamics(t, q, u, p, tau_disturb, f_disturb)
% given state and inputs u from controller, calculate derivatives for ode45
% integrator. 
theta = q(1:3);
omegabody = q(4:6);
posworld = q(7:9); % not used in calculation
vbody = q(10:12); 
f_l = u(1);
tau_c = u(2:4); 

R = rot_matrix(theta);
W = omega2thetadot_matrix(theta);

if nargin < 6, f_disturb = zeros(3,1); end
if nargin < 5, tau_disturb = zeros(3,1); end

% aerodynamic forces and torques
[f_d, tau_d] = fly_aerodynamics(vbody, omegabody, p);

% forces (f) (body coords)
f_g = R' * [0, 0, - p.g * p.m]'; % gravity
f_l = [0, 0, f_l]'; % control force
%f_disturb = [p.force_bias_x, p.force_bias_y, 0]'; 
f = f_l + f_g + f_disturb + f_d; 

% moments/torues (tau) (body coords)
%tau_disturb = [p.torque_bias_x, p.torque_bias_y, 0]'; 
tau = tau_c + tau_d + tau_disturb - p.ks * [theta(1), theta(2), 0]'; 

fictitious_f = p.m * cross(omegabody, vbody); %[0; omega; 0], [q(4); 0; q(6)]); % omega x v
fictitious_tau = cross(omegabody, p.Jmat * omegabody); 

% geometric
xdotworld = R * vbody; 
vdotbody = 1/p.m * (f - fictitious_f); 
thetadot = W * omegabody; 
omegadotbody = p.Jmat \ (tau - fictitious_tau);

qdot = [thetadot; omegadotbody; xdotworld; vdotbody]; 
end

function c = cross(a,b) % fast cross because default cross func is slow
c = [a(2)*b(3)-a(3)*b(2); a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)];
end


function R = rot_matrix(eulerAngles)
% R matrix to convert 3-vectors in body coords to world coords
% v = Rv' where v is in world frame and v' is in body frame


% XYZ/zyx (airplane) convention
% theta1 = thetaz = thetaX; theta2 = thetay = thetaY; theta3 = thetax =
% thetaZ. got it? : )
% no that is dumb. the order is different, but make theta1 always the
% rotation around x. 
cz = cos(eulerAngles(3)); 
cy = cos(eulerAngles(2));
cx = cos(eulerAngles(1));
sz = sin(eulerAngles(3));
sy = sin(eulerAngles(2));
sx = sin(eulerAngles(1));
R = [cz*cy,   cz*sy*sx - cx*sz,   sz*sx + cz*cx*sy; 
     cy*sz,   cz*cx + sz*sy*sx,   cx*sz*sy - cz*sx;
     -sy,    cy*sx,              cy*cx];  

% ZYX (vicon) convention
% convention here is ZYX convention: the three coordinates in euler_theta  
% mean rotate around world Z, then world Y, then world X 
% alternately, rot around body x, then new body y, then new body z
% (the latter two are equivalent)
% theta1 = thetax = thetaZ; theta2 = thetay = thetaY; theta3 = thetaz =
% thetaX    
% c1 = cos(eulerAngles(1)); 
% c2 = cos(eulerAngles(2));
% c3 = cos(eulerAngles(3));
% s1 = sin(eulerAngles(1));
% s2 = sin(eulerAngles(2));
% s3 = sin(eulerAngles(3));
% R = [c2.*c3, -c2.*s3, s2;
%      c1.*s3 + c3.*s1.*s2, c1.*c3 - s1.*s2.*s3, -c2.*s1;
%      s1.*s3 - c1.*c3.*s2, c3.*s1 + c1.*s2.*s3, c1.*c2];
end

function rotatedVectors = rotVectors(eulerAngles, vectors)
% rotate vectors. 
% angles is 3xn or 3x1, vectors is 3xn 
% rotates a vector v' (vectors) given in body-frame coordinates to 
% v (rotatedVectors) given in world-frame coords according to 
% v = Rv'

dim = size(eulerAngles);
n = dim(2);
dim = size(vectors); 
nv = dim(2);
if n == 1
    eulerAngles = repmat(eulerAngles, [1, nv]); 
end
cz = cos(eulerAngles(3,:)); 
cy = cos(eulerAngles(2,:));
cx = cos(eulerAngles(1,:));
sz = sin(eulerAngles(3,:));
sy = sin(eulerAngles(2,:));
sx = sin(eulerAngles(1,:));
R1 = [cz.*cy;   cz.*sy.*sx - cx.*sz;   sz.*sx + cz.*cx.*sy];
R2 = [cy.*sz;   cz.*cx + sz.*sy.*sx;   cx.*sz.*sy - cz.*sx];
R3 = [-sy;    cy.*sx;              cy.*cx];  


% ZYX (vicon) convention
% rows of R
% R1 = [c2.*c3, -c2.*s3, s2];
% R2 = [c1.*s3 + c3.*s1.*s2, c1.*c3 - s1.*s2.*s3, -c2.*s1];
% R3 = [s1.*s3 - c1.*c3.*s2, c3.*s1 + c1.*s2.*s3, c1.*c2];

% do a dot product by summing vertically (along 1st dimension)
rotatedVectors = [
    sum(R1.*vectors,1); ...
    sum(R2.*vectors,1); ...
    sum(R3.*vectors,1)]; 
end


function transformedVectors = transformVectors(eulerAngles, translation, vectors)
% rotate and translate vectors. 
% eulerAngles is 3xn or 3x1, translation is 3xn or 3x1, vectors is 3xn 
% rotates a vector v' (vectors) given in body-frame coordinates 
% a vector v to world-frame coords  
% and translates it by t (translation)
% v = Rv' + t

dim = size(translation);
n = dim(2);
if n == 1
    dim = size(vectors); 
    nv = dim(2);
    translation = repmat(translation, [1, nv]); 
end
transformedVectors = rotVectors(eulerAngles, vectors) + translation; 
end

function W = omega2thetadot_matrix(euler_theta)
% transform euler angle rates to angular rot rate vector omega
% thetadot = W * omega
% this needs tobe corrected for 3d case for the current euler angle convention. 
st1 = sin(euler_theta(1));
ct1 = cos(euler_theta(1));
tt2 = tan(euler_theta(2));
ct2 = cos(euler_theta(2));

% XYZ (airplane) convention
W = [1, st1*tt2, ct1*tt2; 
    0, ct1, -st1;
    0, st1/ct2, ct1/ct2];

% ZYX (vicon) convention
% needs fixing
% W = [c3/c2,   -s3/c2,     0; 
%      s3,      c3,         0; 
%      -c3*t2,  s3*t2,      1];  

end


function [f_d, tau_d] = fly_aerodynamics(v, omega, p)
% calculate stroke-averaged forces due to aerodynamic drag on flapping 
% wings. assumes force is applied at point r_w away from center of mass in
% z-direction in body-coordinates. v_w is the vel at that point. 
% assumes drag force f_d = -b_w * v_w. this has been tested in a wind
% tunnel for x- and y-directions but not z-direction
r_w = [0, 0, p.r_w]';
v_w = v + cross(omega, r_w); % vel at midpoint of wings
f_d = - p.b_w * v_w;  
tau_d = cross(r_w, f_d); %p.r_w * f_d(1); % r_w cross f_d
end

function render_body(eulerAngles, position, and_wings)

if nargin < 3
    and_wings = 1;
end
if nargin < 2
    position = zeros(3,1); 
end
if nargin < 1
    eulerAngles = zeros(3,1); 
end

persistent moving_geometry fixed_geometry 
if isempty(moving_geometry)
    % initialize persistent vars
    b.l = 3e-3; 
    b.w = 2e-3; 
    b.h = 13e-3; 
    moving_geometry.body.patch = patch('EraseMode','normal');
    moving_geometry.body.verts = [
        b.l/2 * [-1, -1, 1, 1, -1, -1, 1, 1]; % length
        b.w/2 * [-1, 1, 1, -1, -1, 1, 1, -1]; % width
        b.h/2 * [-1, -1, -1, -1, 1, 1, 1, 1]]; % height
    moving_geometry.body.faces = [1 2 3 4; 2 6 7 3;  4 3 7 8; ...
                                  1 5 8 4; 1 2 6 5; 5 6 7 8];
    moving_geometry.shadow.patch = patch('EraseMode','normal');
    moving_geometry.shadow.verts = [
        -b.l/2, -b.w/2, .001; 
        b.l/2, -b.w/2, .001; 
        b.l/2, b.w/2, .001; 
        -b.l/2, b.w/2, .001]';
    moving_geometry.shadow.faces = 4:-1:1; 
    if and_wings
        vertr = 1e-3*[14.9           25         27.5         28.4           29         29.2         29.3         29.3         29.2         28.8         27.7         25.7         22.6         21.1         20.3         19.3         18.4         17.5         14.9;
                           35.2         35.2         35.2         35.1         34.9         34.7         34.6         34.5         34.3         33.8         33.1           32         30.7         30.2         30.1         30.1         30.4         31.1         33.6];
        % put into x-y-z coordinates
        vertr = [zeros(1, length(vertr)); vertr];
        % recenter and position relative to fly correctly
        vertr = vertr ...
            - repmat(vertr(:,1), [1, length(vertr)]) ...
            + repmat([0; b.w/2 + 0.05e-3; b.h/2], [1, length(vertr)]); 
        vertl = [vertr(1,:); -vertr(2,:); vertr(3,:)]; 
        moving_geometry.rightwing.patch = patch('EraseMode','normal');
        moving_geometry.rightwing.verts = vertr; 
        moving_geometry.rightwing.faces = length(vertr):-1:1;
        
        moving_geometry.leftwing.patch = patch('EraseMode','normal');
        moving_geometry.leftwing.verts = vertl;
        moving_geometry.leftwing.faces = length(vertl):-1:1;
        
    end
    fixed_geometry.floor.patch = patch('EraseMode','normal');
    fs = 50e-3; % floor size m
    fixed_geometry.floor.verts = [-fs, -fs, 0; fs, -fs, 0; fs, fs, 0; -fs, fs, 0]'; 
    fixed_geometry.floor.faces = 4:-1:1; 
end
    
fieldnames_ = fieldnames(moving_geometry); 
for idx = 1:length(fieldnames_)
    fieldname = fieldnames_{idx};
    verts = moving_geometry.(fieldname).verts; 
    if strcmp(fieldname, 'shadow')
        % special case for shadow, only rotate aroud z and don't shift z
        shadow_eulerAngles = [0; 0; eulerAngles(3)];
        shadow_position = [position(1:2); 0]; 
        transformed_verts = transformVectors(shadow_eulerAngles, shadow_position, verts);
        color = 'black'; 
    else
        transformed_verts = transformVectors(eulerAngles, position, verts); 
        color = 'white'; 
    end
    set(moving_geometry.(fieldname).patch, ...
        'Faces', moving_geometry.(fieldname).faces, ...
        'Vertices',transformed_verts', ...
        'FaceColor',color)
end

fieldnames_ = fieldnames(fixed_geometry); 
for idx = 1:length(fieldnames_)
    fieldname = fieldnames_{idx};
    set(fixed_geometry.(fieldname).patch, ...
        'Faces', fixed_geometry.(fieldname).faces, ...
        'Vertices',fixed_geometry.(fieldname).verts', ...
        'FaceColor','white')
end
     
%axis tight
axis([-.06 .06 -.06 .06 0 .09])
end

    