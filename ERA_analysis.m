%%%%%%%%%%%%%%%%%%%%
% Khushant Khurana %
       % ERA %
%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Initializing parameters for the Atomic Force Microscope TF
zeta1 = 0.03;
zeta2 = 0.03;
zeta3 = 0.042;
zeta4 = 0.025;
zeta5 = 0.032;
k = 5;
tau = 10^(-4);
w1 = 2*pi*2.4*1000;
w2 = 2*pi*2.6*1000;
w3 = 2*pi*6.5*1000;
w4 = 2*pi*8.3*1000;
w5 = 2*pi*9.3*1000;
nin = 1; % Number of inputs
nout = 1; % Number of outputs

%% Declaring the tf
s = tf('s');
num = k*w2^2*w3^2*w5^2*(s^2 + 2*zeta1*w1*s + w1^2)*(s^2 + 2*zeta4*w4*s + w4^2)*exp(-s*tau);
denom = w1^2*w4^2*(s^2 + 2*zeta2*w2*s + w2^2)*(s^2 + 2*zeta3*w3*s + w3^2)*(s^2 + 2*zeta5*w5*s + w5^2);

%% Impulse respsone to the tf
g = num/denom;
[y,t] = impulse(g);

% Determining the delta_t for the impulse response. 
% Would be used for the discrete case later!
dt = t(2) - t(1);

%% Plotting the impulse response!
figure
plot(t,y,LineWidth=2)
xlabel("Time [sec]")
ylabel("Y")
title("Impulse response of Atomic Force Microscope")
grid on
box off

%% Constructing the hankel matrices
H = hankel(y(1:end-1));
H2 = hankel(y(2:end)); % with one time shift

%% Performing svd on the Hankel matrix!
[U,S,V] = svd(H,'econ');

%% Performing ERA!
r = 100; % the reduced size. Can be altered!

Sigma = S(1:r,1:r);
Ur = U(:,1:r);
Vr = V(:,1:r);

% Deriving the A,B,C matrices
Ar = Sigma^(-.5)*Ur'*H2*Vr*Sigma^(-.5);
Br = Sigma^(-.5)*Ur'*H(:,1:nin);
Cr = H(1:nout,:)*Vr*Sigma^(-.5);

%% Plotting frequency repsonse comparison between original and reconstructed system

% Utility to allow for phase wrapping in bode plots
opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';

% Frequency response of reconstructed system!
figure
subplot(1,2,1)
bode(dt*ss(Ar,Br,Cr,0,ts = dt),opts)
title("ERA approximation")
grid on

% Frequency response of original system
subplot(1,2,2)
bode(g,opts)
title("Original")
grid on
