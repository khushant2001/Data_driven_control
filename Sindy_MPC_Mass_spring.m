%{
Doing Sindy on a Mass Spring System
%}

clear
clc

%% Parameters
k = 10;
b = 20;
m = 1;

%% State space equations
A = [0,1;-k/m,-b/m];
B = [0;1/m];
C = [1,0;0,1];
D = 0;

%% Lsim response to sinusoidal force input
sys = ss(A,B,C,D);
t = linspace(0,10,100);
dt = 10/100;
u = sin(t);
y = lsim(sys,u,t);

%% Constructing functions for Sindy
y_dot = diff(y)/dt;
y(end,:) = [];

%% Build library and compute sparse regression

n=2;
Theta = [ones(length(y),1),y,u(1:end-1)',y(:,1).^2,y(:,2).^2,u(1:end-1)'.^2];
lambda = .01; % lambda is our sparsificat ion knob.
Xi = sparsifyDynamics(Theta,y_dot, lambda,n);
functions = ["1","x","x_dot","u","x^2","x_dot^2","u^2"]';
heading = ["Basis Functions","Velocity","Acceleration"];
coefficients = [functions,Xi];
coefficients = vertcat(heading,coefficients);
disp(coefficients)

%% Constructing the new differential equations:
x_dot_coeff = Xi(:,1);
x_ddot_coeff = Xi(:,2);

%{
%% LQR Controller Design
Q = [10,0;0,1];
R = [10];
k = lqr(A_id,B_id,Q,R);
%}
%% Sindy Function Definitions

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