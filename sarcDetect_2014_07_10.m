% Performs the detection of the sarcomeres in an image.
%
% The outputs are:
%  - Binary skeleton image of restored sarcomeres
%  - *.sarcDetect.Settings.mat file containing sarcomere alignment data and
%  parameters needed to run sarcDetectMulti
%
% Original version by Adam W. Feinberg, Feb 2009
% School of Egineering and Applied Sciences
% Havard University, Cambridge, MA 02138
%
% Last updated October 2019 by Davi M. Lyra Leite
% Viterbi School of Engineering
% University of Southern California, Los Angeles, CA 90089


%% ========================================================================
% Clearing workspace
clear all
close all
clc

% Changing the path to add the finger print detection toolbox
addpath(genpath('Finger Print Detection/'));


%% ========================================================================
% User inputs file:
[file,path] = uigetfile({'*.TIF';'*.tif';'*.bmp';'*.jpg';'*.*'},...
    'Select Image File...');
filename = [path file];
disp(filename)
disp(' ')

% Read the image
im = imread(filename);
info = imfinfo(filename);
width = info.Width;
height = info.Height;
disp(['width = ',num2str(width)])
disp(['height = ',num2str(height)])
disp(' ')


%% ========================================================================
% Microscope selection routine to define appropriate sarcomere spacing
micselection;


%% ========================================================================
% Identify ridge-like regions and normalize image
index = 0;
while index < 1;
    % Threshold of standard deviation to decide if a block is a ridge 
    % region
    thresh = input('\nEnter Threshold (0.01 - 0.20):\n');
    disp('Normalizing Image and Creating Mask')
    disp(' ')
    [normim, mask] = ridgesegment(im, blksze, thresh);
    show(normim,1);
    show(mask, 2);
    
    % Determine if normalization and mask look good, click on image to
    % accept or press any key to enter new values
    w = input('\nAccept Threshold (yes, no): ','s');
    disp(' ')
    if strcmp(w,'y') || strcmp(w,'Y') || strcmp(w,'yes') || strcmp(w,'Yes')
        disp('Image threshold accepted')
        index = 1;
    else %if w == 'n' || w == 'N' || strcmp(w,'no') || strcmp(w,'No')
        disp('Re-analyze imaging...')
        index = 0;
    end
end

% Create skeleton image from normalized image
norm_bin_skel = maskdata(normim,mask);
show(norm_bin_skel, 3);

% Determine ridge orientations
disp('Calculating Ridge Orientations' )
[orientim, reliability] = ridgeorient(normim, 1, 3, 3);

% Add orientation values by pi*3/2 for direction orthoganol to
% sarcomere in the direction of force generation
orientim_perp = orientim + (pi*3/2);
plotridgeorient(orientim_perp, 10, im, 4)
show(reliability,5)


%% ========================================================================
% Determine ridge frequency values across the image
disp('Calculating Ridge Frequency Values' )
disp(' ')
[~, medfreq] = ridgefreq(normim, mask, orientim, blksze, 3,...
    MinSarcSpacing, MaxSarcSpacing);

freq = medfreq.*mask;

index = 0;

% kx controls the sigma in the x direction which is along the filter, and
% hence its bandwidth.
kx = .4;
% ky controls the sigma across the filter and hence controls the
% orientational selectivity of the filter.
ky = .4;
% A value of 0.4 for both kx and ky is a good starting point, therefore, it
% is the selected default value.

while index < 1;
    newim = ridgefilter(normim, orientim, freq, kx, ky, 1);
    show(newim,6);
    
    binim = newim > 0;
    binim_skel = bwmorph(binim,'skel',Inf);
    binim_skel = bwmorph(binim_skel,'clean',Inf);
    show(binim_skel,7);
    
    b = input('Accept Sarcomere Detection (yes, no): ','s');
    if strcmp(b,'y') || strcmp(b,'Y') || strcmp(b,'yes') || strcmp(b,'Yes')
        disp('Image accepted' )
        index = 1;
    else
        disp('Re-analyze imaging...')
        kx = input('Enter Sigma along filter (default 0.4):\n');
        ky = input('Enter Sigma across filter (default 0.4):\n');
    end
end

% save filtered image of sarcomere skeleton before manual removal of
% non-sarcomeres
type = '_sarcomere_FULL.tif';
filename2 = [filename type];
imwrite(binim_skel,filename2,'Compression','none');


%% ========================================================================
% Remove false sarcomeres at tissue borders by selecting regions to be
% masked using ROI. ROI's selected will be REMOVED!!!
index = 0;
hold on
while index < 1;
    disp('Select ROI to exclude from further analysis');
    disp('(double-click to close the ROI)');
    BW = roipoly(binim_skel);
    BW2 = ~BW;
    binim_skel = binim_skel.*BW2;
    
    % Use mask to remove false sarcomeres from binary skeleton
    show(binim_skel)
    b = input('Select another ROI to exclude? (yes, no): ','s');
    if strcmp(b,'y') || strcmp(b,'Y') || strcmp(b,'yes') || strcmp(b,'Yes')
        disp('Select another ROI...')
        disp(' ')
    else
        disp('Image accepted' )
        disp(' ')
        index = 1;
    end
end
hold off


%% ========================================================================
% Multiply orientation angles by the binary skeleton image to remove
% false sarcomeres
orientim = orientim.*binim_skel;
orientim_perp = orientim_perp.*binim_skel;

% Convert 2D-array to 1D vector
orientation = orientim_perp(:);

% Keep non-zero values only
nonzero_orientation = orientation((orientation)~=0) - 2*pi;

% Find sum of orientation unit vectors in (X) horizontal and (Y) vertical
% directions
x_vector = sum(cos(nonzero_orientation));
disp(x_vector); disp(' ')
y_vector = sum(sin(nonzero_orientation));
disp(y_vector); disp(' ')

% Convert radians to degrees
nonzero_orientation_angles = rad2deg(nonzero_orientation);

% Analyze orientational data
Mean = mean(nonzero_orientation_angles);
disp(Mean); disp(' ')
Std = std(nonzero_orientation_angles);
disp(Std); disp(' ')
Median = median(nonzero_orientation_angles);
disp(Median); disp(' ')

% Create histogram
[n,xout] = hist(nonzero_orientation_angles,180);
dx = xout(2)-xout(1); % calc a single bin width
n = n / sum( n*dx ); % normalize histogram to have area of 1

% Find mode
[C,I] = max(n);
Mode = xout(I);
disp(Mode); disp(' ')

% Plot normalized histogram of raw orientation
figure, bar(xout,n,'hist')

% Making sure that the axis is squeezed to its limits
xlim( [xout(1)-dx/2,xout(end)+dx/2] );
title('Histogram of Sarcomere Orientation Angles')
xlabel('Degrees')
ylabel('Normalized Occurence')


%% ========================================================================
% Total number of sarcomere positive pixels in the skeleton image
% Sarcomere density = total/(image area)
TotalSarc = length(nonzero_orientation_angles);
disp(TotalSarc); disp(' ')


%% ========================================================================
% save filtered image of sarcomere skeleton
type = '_sarcomere_CLEAN.tif';
filename2 = [filename type];
imwrite(binim_skel,filename2,'Compression','none');

% save settings file for sarcDetectMulti
ext_settings = '_sarcDetect_Settings.mat';
filename3 = [filename ext_settings];
save(filename3,'blksze','thresh','MinSarcSpacing','MaxSarcSpacing',...
    'kx','ky','pix2um','width','height','orientim','orientim_perp',...
    'nonzero_orientation');

clear
