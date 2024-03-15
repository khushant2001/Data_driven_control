%{
%%%%%%%%%%%%%%%%%%%%%%%%%
    % Khushant Khurana %
%%%%%%%%%%%%%%%%%%%%%%%%%
%}

clear;
clc;
close all;

% Import files
filenameGrid = fullfile('airfoilDNS_grid.h5');
filename_parameters = fullfile('airfoilDNS_a25f0p05.h5');

% Get the velocities and time_parameters for all snapshots taken. 
ux = h5read(filename_parameters, '/ux');    
uy = h5read(filename_parameters, '/uy'); 
t_field = h5read(filename_parameters,'/t_field');

% Getting the grid points for the flow regime. 
x = h5read(filenameGrid,'/x');
y = h5read(filenameGrid,'/y');

% Getting the lengths of the arrays which will be used later for
% reshaping purposes. 
nx = length(x);
ny = length(y);
nt = length(t_field);

% Reshaping the 2d velocities for different snapshots into column vectors. 
uxreshape = reshape(ux,nx*ny,nt);
uyreshape = reshape(uy,nx*ny,nt);
data = [uxreshape;uyreshape];


%x = linspace(0,10,1000);
%data = x.^2;
% Subtracting the average mean, which serves as the purpose of
% normalization. 
meanSub = 1;
if meanSub
    dataMean = mean(data,2);
    data = data-dataMean*ones(1,nt);
end

% Performing SVD!
[U,Z,V] = svd(data,'econ');
temporal_matrix = Z * V';
X_data = temporal_matrix(1:6,1:length(Z));

%% Calculating derivative of data
%X_dot = gradient(X_data);
X_dot = diff(X_data,1,2)/.01;
%X_dot(end,:) = [];
X_data(:,end) = [];
%% Build library and compute sparse regression
n = 3;
Theta = poolData(X_data,n,2); % up to third order polynomials
lambda = .0139; % lambda is our sparsificat ion knob.
Xi = sparsifyDynamics(Theta,X_dot, lambda,n);
yout=poolDataLIST({'x','y','z'},Xi,n,3);

new_x_dot = Theta*Xi;
figure

for i = 1:6
    y_new = cumtrapz(.01,new_x_dot(i,:));
    y_new = y_new + X_data(i,1);
    subplot(3,3,i)
    plot(t_field(1:end-1),y_new, lineWidth = 2)
    hold on 
    plot(t_field(1:end-1), X_data(i,:), LineWidth = 2)
    legend('SINDY Approximation', 'Original')
    title(["Temporal amplitude - ", num2str(i), " approximation using Sindy"])
    ylabel("Temporal Amplitude")
    xlabel("Snapshots")
    grid on
end

%% Function definitions
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
function yout = poolData(yin,nVars,polyorder)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

ind = 1;
% poly order 0
yout(:,ind) = ones(n,1);
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = yin(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = yin(:,i).*yin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end
end

function yout = poolDataLIST(yin,ahat,nVars,polyorder)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);

ind = 1;
% poly order 0
yout{ind,1} = ['1'];
ind = ind+1;

% poly order 1
for i=1:nVars
    yout(ind,1) = yin(i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout{ind,1} = [yin{i},yin{j}];
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout{ind,1} = [yin{i},yin{j},yin{k}];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout{ind,1} = [yin{i},yin{j},yin{k},yin{l}];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout{ind,1} = [yin{i},yin{j},yin{k},yin{l},yin{m}];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end




output = yout;
newout(1) = {''};
for k=1:length(yin)
    newout{1,1+k} = [yin{k},'dot'];
end
% newout = {'','xdot','ydot','udot'};
for k=1:size(ahat,1)
    newout(k+1,1) = output(k);
    for j=1:length(yin)
        newout{k+1,1+j} = ahat(k,j);
    end
end
newout;
end
