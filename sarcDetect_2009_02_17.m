%             
%             Binary skeleton image of restored sarcomeres after ROI based
%              removal of non-sarcomers
%             
%             *.sarcDetect.Settings.mat file containing sarcomere alignment
%               data and parameters needed to run sarcDetectMulti on a folder
%               full of image files6gineering and Applied Sciences
% Havard University, Cambridge, MA 02138

% Last updated Feb 2009 by Adam W. Feinberg

   % ask user to select file
   [file,path]=uigetfile({'*.TIF';'*.tif';'*.bmp';'*.jpg';'*.*'},...
       'Select Image File...');%,...
       %'/Users/davileite/Desktop/LSE Data/Processed Images');
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
   % Add1 orientation values by pi*3/2 for direction orthoganol to
   % sarcomere in the direction of force generation
   orientim_perp = orientim + (pi*3/2);
   plotridgeorient(orientim_perp, 10, im, 4)
   show(reliability,5)
 
   % Determine ridge frequency values across the image
   disp('Calculating Ridge Frequency Values' )
   [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 3, MinSarcSpacing, MaxSarcSpacing);

   % Actually I find the median frequency value used across the whole
   % fingerprint gives a more satisfactory result...
   freq = medfreq.*mask;

   hold off;
   index = 0;
   while index < 1;
       % kx controls the sigma in the x direction which is along the
       % filter, and hence controls the bandwidth of the filter.
       kx = input('Enter Sigma along filter (0.4): ');
       % ky controls the sigma across the filter and hence controls the
       % orientational selectivity of the filter. A value of 0.4 for both
       % kx and ky is a good starting point.
       ky = input('Enter Sigma across filter (0.4): ');
       newim = ridgefilter(normim, orientim, freq, kx, ky, 1);
       show(newim,6);
       binim = newim > 0;
       binim_skel = bwmorph(binim,'skel',Inf);
       binim_skel = bwmorph(binim_skel,'clean',Inf);
       show(binim_skel,7);
       b = input('Accept Sarcomere Detection (yes = 0, no = 1): ');
       if b == 0
           disp('Image accepted' )
           index = 1;
       else
           disp('Re-analyze imaging...')
       end
   end

   % save filtered image of sarcomere skeleton before manual removal of
   % non-sarcomeres
   type = '.sarcomere.FULL.tif';
   filename2 = [filename type];
   imwrite(binim_skel,filename2,'Compression','none');

   % Remove false sarcomeres at tissue borders by selecting regions to be
   % masked using ROI.  ROI's selected will be REMOVED!!!
   index = 0;
   hold on
   while index < 1;
       disp('Select ROI to exclude from further analysis (double-click to close the ROI)');
       BW = roipoly(binim_skel);
       BW2 = ~BW;
       binim_skel = binim_skel.*BW2;
       % Use mask to remove false sarcomeres from binary skeleton 
       show(binim_skel)
       b = input('Select another ROI to exclude? (yes = 0, no = 1): ');
       if b == 0
           disp('Image accepted' )
       else
           disp('Select another ROI...')
           index = 1;
       end
   end
   hold off
   
   % Multiply orientation angles by the binary skeleton image to remove
   % false sarcomeres
   orientim = orientim.*binim_skel;
   orientim_perp = orientim_perp.*binim_skel;
   % Convert 2D-array to 1D vector
   orientation = orientim_perp(:);
   % Keep non-zero values only
   nonzero_orientation = orientation(find(orientation))-2*pi;
   
   % Find sum of orientation unit vectors in (X) horizontal and (Y) vertical directions
   x_vector = sum(cos(nonzero_orientation))
   y_vector = sum(sin(nonzero_orientation))
   
   % Convert radians to degrees
   nonzero_orientation_angles = rad2ang(nonzero_orientation);

   % hist(nonzero_orientation_angles)
   Mean = mean(nonzero_orientation_angles)
   Std = std(nonzero_orientation_angles)
   Median = median(nonzero_orientation_angles)
   
   % Create histogram
   [n,xout] = hist(nonzero_orientation_angles,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1
   
   % Find mode
   [C,I] = max(n);
   Mode = xout(I)
      
   % Plot histogram of raw orientation
   figure, bar(xout,n,'hist')              % plot normalized histogram
   xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
   title('Histogram of Sarcomere Orientation Angles')
   xlabel('Degrees')
   ylabel('Normalized Occurence')
   
   % Total number of sarcomere positive pixels in the skeleton image
   % Sarcomere density = total/(image area)
   TotalSarc = length(nonzero_orientation_angles)

   % save filtered image of sarcomere skeleton
   type = '.sarcomere.CLEAN.tif';
   filename2 = [filename type];
   imwrite(binim_skel,filename2,'Compression','none');
   
   % save settings file for sarcDetectMulti
   ext_settings = '.sarcDetect.Settings.mat';
   filename3 = [filename ext_settings];
   save(filename3,'blksze','thresh','MinSarcSpacing','MaxSarcSpacing','kx','ky','pix2um','width','height','orientim','orientim_perp','nonzero_orientation'); 

   clear
   