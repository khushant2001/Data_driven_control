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
meanSub = 1;

% Subtracting the average mean, which serves as the purpose of
% normalization. 
if meanSub
    dataMean = mean(data,2);
    data = data-dataMean*ones(1,nt);
end

% Performing SVD!
[U,Z,V] = svd(data,'econ');
modes = linspace(0,401,401);

% Finding the importance of individual eigenvalues for different snapshots.
figure;
semilogy(modes,(diag(Z).^2), linewidth = 2);
xlabel("Modes");
ylabel("Eigenvalues");
title("Lin log scale");

U_ux = U(1:length(uxreshape), :);
U_uy = U(length(uyreshape)+1:end, :);

MM = 0.01; % defining contour levels on the plot. 
v = -1:0.1:1;
v(11)=[];

% Plotting the spatial modes for both velocities. 
figure;
for k = 1:6
    subplot(2, 3, k);
    contourf(x, y, transpose(reshape(U_ux(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('ux spatial modes')

figure;
for k = 1:12
    subplot(3, 4, k);
    contourf(x, y, transpose(reshape(U_uy(:,k), nx, ny)), MM * v);
    caxis([-MM MM]);
    colorbar;
end
sgtitle('uy spatial modes')

temporal_matrix = Z * V'; % This matrix contains information for the temporal amplitude. 

% Plotting the ux temporal nodes

x_axis_2 = linspace(0,length(Z),length(Z));
% Repeat for the next five subplots
figure;
for i = 1:6
    subplot(2,3,i);
    % Placeholder plot
    plot(x_axis_2, temporal_matrix(i,1:length(x_axis_2)));
    xlim([0,400]);
    title(['ux temporal amplitude ', num2str(i)]);
end

% Plotting the uy temporal nodes

x_axis_2 = linspace(0,length(Z),length(Z));
figure;
% Repeat for the next five subplots
for i = 1:6
    subplot(2,3,i);
    % Placeholder plot
    plot(x_axis_2, temporal_matrix(i+201,1:length(x_axis_2)));
    xlim([0,400]);
    title(['uy temporal amplitude ', num2str(i)]);
end

cycles = 1.5;
temporal_freq = 1.5/400;
disp("Frequency of the temporal amplitude = ", num2str(temporal_freq));
% Lets reconstruct the flow regime.
r = 20; % Parameter used for reconstruction. Higher the value of r, more
% information is used for reconstruction. 
new_matrix = U(:,1:r) * Z(1:r,1:r) * V(:,1:r)' + repmat(dataMean, 1, nt);

% Plotting the new reconstructed flow regimes. 
figure;
z = transpose(reshape(real(new_matrix(end/2+1:end,1)), nx, ny));
contourf(x, y, z);
title("Reconstruction of the y-velocity");
colorbar;

figure;
contourf(x, y, transpose(reshape(new_matrix(1:end/2,1), nx, ny)));
title("Reconstruction of the x-velocity");
colorbar;

% The following code can be commented out to get the video. IMP! DOING SO
% WILL TAKE ABOUT 25 MIN. 
%{
%% Making a video
% Setup for video creation
videoFilename = 'FluidFlowReconstruction.avi';
v = VideoWriter(videoFilename, 'Uncompressed AVI'); % Choose your desired format
v.FrameRate = 10; % Adjust frame rate as needed
open(v);

% Loop through each time step
for t = 1:nt
    % Create a new figure for each frame to ensure consistent sizing
    fig = figure('Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Visible', 'off');

    % Original flow plot
    subplot(1, 5, 1);
    imagesc(x, y, squeeze(uy(:,:,t))');
    axis equal tight;
    title('Original Flow');
    colorbar;

    % Ranks for reconstruction
    ranks = [2, 4, 10, 20];

    for i = 1:length(ranks)
        r = ranks(i);
        % Reconstruct flow using the first 'r' modes
        data_approx = U(:,1:r) * Z(1:r,1:r) * V(:,1:r)' + repmat(dataMean, 1, nt); % Add back the mean
        uy_approx = reshape(real(new_matrix(end/2+1:end,t)), nx, ny);

        % Plot reconstruction
        subplot(1, 5, i+1);
        imagesc(x, y, uy_approx');
                                                                                                                                                                                                                                axis equal tight;
        title(['Reconstruction r=', num2str(r)]);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        colorbar;
    end

    drawnow;

    % Capture and write frame
    frame = getframe(fig);   
    writeVideo(v, frame);

    % Close the figure to free up memory
    close(fig);
end

close(v);
%}