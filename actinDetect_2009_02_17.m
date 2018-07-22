% actinDetect
%
% Code for automatic detection of actin filament alignment in cells
% stained for F-actin (such as phalloidin conjugated to FITC)
%
% Adapted from Function to demonstrate use of fingerprint code
%
% Argument:   Load image of actin stained cells, file should
%             TIF format grayscale at least 8-bit depth
%
% Returns:    *.actinDetect.Settings.mat file containing actin alignment
%               data and parameters needed to run actinDetectMulti on a folder
%               full of image files

% Adapted from:
% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
% January 2005
%
% Created by Adam W. Feinberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Havard University, Cambridge, MA 02138

% Last updated Feb 2009 by Adam W. Feinberg

   % ask user to select file
   [file,path]=uigetfile({'*.TIF';'*.tif';'*.bmp';'*.jpg';'*.*'},'Select Image File...','D:/Manuscripts/2D Tissue Manuscripts/');
   filename = [path file];
       disp(filename)
   im = imread(filename);
   info = imfinfo(filename);
   width = info.Width
   height = info.Height

   % Enter pixel to micrometer converesion factor to define appropriate
   % sarcomere spacing
   pix2um = input('Enter calibration factor in pixels per micrometer (Coolsnap 63x=9.8, 40x=6.22, 20x=3.085): ');
   SarcSpacing = 2*pix2um;          % convert 2 um sarcomere spacing into pixels
   MinSarcSpacing = round(0.5*SarcSpacing); % minimum length of sarcomere spacing
   MaxSarcSpacing = round(1.5*SarcSpacing); % maximum length of sarcomere spacing
   blksze = MaxSarcSpacing;
   
   % Identify ridge-like regions and normalise image
   index = 0;
   while index < 1;
%        % Input blocksize over which the standard deviation is determined
%        blksze = input('Enter blocksize (16): ');
       % Threshold of standard deviation to decide if a block is a ridge region
       thresh = input('Enter Threshold (0.1 - 0.2): ');
       disp('Normalizing Image and Creating Mask' )
       [normim, mask] = ridgesegment(im, blksze, thresh);
       show(normim,1);
       show(mask, 2);
       % Determine if normalization and mask look good, click on image to
       % accept or press any key to enter new values
       w = input('Accept Threshold (yes = 0, no = 1): ');
       if w == 0
           disp('Image threshold accepted' )
           index = 1;
       else
           disp('Re-analyze imaging...')
       end
   end

   % Create skeleton image from normalized image
   normim_bin = normim > 0;
   mask_bin = mask > 0;
   normim_mask = normim_bin.*mask_bin;
   norm_bin_skel = bwmorph(normim_mask,'skel',Inf);
   show(norm_bin_skel, 3);

   % Determine ridge orientations
   disp('Calculating Ridge Orientations' )
   [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
   plotridgeorient(orientim, 10, im, 4)
   show(reliability,5)
 
   % Convert 2D-array to 1D vector
   orientation = orientim(:);
      
   % Find sum of orientation unit vectors in (X) horizontal and (Y) vertical directions
   x_vector = sum(cos(orientation))
   y_vector = sum(sin(orientation))
   
   % Convert radians to degrees
   orientation_angles = rad2deg(orientation);

   % hist(nonzero_orientation_angles)
   Mean = mean(orientation_angles)
   Std = std(orientation_angles)
   Median = median(orientation_angles)
   
   % Create histogram
   [n,xout] = hist(orientation_angles,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1
   
   % Find mode
   [C,I] = max(n);
   Mode = xout(I)
      
   % Plot histogram of raw orientation
   figure, bar(xout,n,'hist')              % plot normalized histogram
   xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
   title('Histogram of Actin Orientation Angles')
   xlabel('Degrees')
   ylabel('Normalized Occurance')
   
   % Total number of actin positive pixels in the skeleton image
   % Sarcomere density = total/(image area)
   TotalActin = length(orientation_angles)

   % save settings file for sarcDetectMulti
   ext_settings = '.actinDetect.Settings.mat';
   filename3 = [filename ext_settings];
   save(filename3,'blksze','thresh','MinSarcSpacing','MaxSarcSpacing','pix2um','width','height','orientim','orientation'); 

   clear
   