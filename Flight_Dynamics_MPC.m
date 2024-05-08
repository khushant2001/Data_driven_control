%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Khushant Khurana %
   % Flight Dynamics MPC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Note: Heavility Inspired by Dr. Brunton's Code!
The dataset can be found here: https://c3.ndc.nasa.gov/dashlink/resources/294/
%}


clc
close all 
clear

%% Importing Dataset!
[data_import,~,~] = xlsread("flight_data.xlsx");
time = data_import(:,1);
dt = time(end)/length(data_import);
data = data_import(:,5:10);
u = data_import(1:end-1,11:end);

%% Constructing functions for Sindy
data_dot = diff(data,1,1)/dt;
data(end,:) = [];

%% Build library and compute sparse regression

n=6;
Theta = [ones(length(data),1),data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,4).*data(:,5),data(:,5).*data(:,6),data(:,4).*data(:,6),sin(data(:,1)),sin(data(:,2)),sin(data(:,3)),cos(data(:,1)),cos(data(:,2)),cos(data(:,3)),u(:,1),u(:,2),u(:,3),u(:,4),u(:,5),u(:,5).^2];
lambda = .1; % lambda is our sparsificat ion knob.
Xi = sparsifyDynamics(Theta,data_dot, lambda,n);
functions = ["1","phi","theta","psi","p","q","r","pq","qr","pr","sin(phi)","sin(theta)","sin(psi)","cos(phi)","cos(theta)","cos(psi)","delta_a","delta_e","delta_r","delta_t","V_a","V_a^2"]';
heading = ["Basis Functions","phi_dot","theta_dot","psi_dot","p_dot","q_dot","r_dot"];
coefficients = [functions,Xi];
coefficients = vertcat(heading,coefficients);
disp(coefficients)

%% Plotting data to optimize the lambda. 

new_x_dot = Theta*Xi;
for i = 1:6
    y_new = cumtrapz(.01,new_x_dot(:,i));
    y_new = y_new + data(1,i);
    subplot(3,3,i)
    plot(time(1:end-1),y_new, lineWidth = 2)
    hold on 
    plot(time(1:end-1), data(:,i), LineWidth = 2)
    legend('SINDY Approximation', 'Original')
    title([functions(i+1), " approximation using Sindy"])
    ylabel("Attitude info [deg]")
    xlabel("Time [sec]")
    grid on
end

%% Designing MPC 
%{
Ts = dt; % Sampling time (adjust as needed)
N = 10; % Prediction horizon
M = 2; % Control horizon

% Create MPC object
mpcobj = mpc(new_x_dot, Ts, N, M);

% Specify constraints, weights, and other MPC parameters
mpcobj.MV = struct('Min', MinValue, 'Max', MaxValue);
mpcobj.OutputVariables = 'omega'; % Assuming you want to control angular velocity

% Step 4: Simulation and Validation
% Simulate MPC controller with your system model and data
simulator = sim(mpcobj, omega_data); % Simulate with angular velocity data

% Plot results
plot(simulator);
%}
%% Function Definitions!

function Xi = sparsifyDynamics(Theta,dXdt , lambda,n)
    % Compute Sparse regression: sequent ial least squares
    Xi = Theta\dXdt ; % Initial guess: Least-squares
    % Lambda is our sparsi f icat ion knob.
    for k=1:10
        smallinds = (abs(Xi )<lambda) ; % Find smal l coeff icients
        Xi(smallinds)=0; % and threshold
        for ind = 1:n % n is state dimension
            biginds = ~smallinds( : , ind) ;
            % Regress dynamics onto remaining terms to f ind sparse Xi
            Xi (biginds, ind) = Theta( : ,biginds) \dXdt ( : , ind) ;
        end
    end
end
