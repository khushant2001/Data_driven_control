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


data = importdata("AIL_test_M02.mat");
data_matrix = [];

for i = 1:25
    column = reshape(data{i}, 180000,1);
    data_matrix(:,i) = column;
end

meanSub = 0;
if meanSub
    dataMean = mean(data_matrix,2);
    data = data_matrix-dataMean*ones(1,nt);
end

%% Calculating derivative of data
X_dot = diff(data_matrix,1,1)/.01;
data_matrix(end,:) = [];

%% Build library and compute sparse regression
n = 4;
Theta = poolData(data_matrix,n,2); % up to third order polynomials
lambda = .0139; % lambda is our sparsificat ion knob.
Xi = sparsifyDynamics(Theta,X_dot, lambda,n);
new_x_dot = Theta*Xi;

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