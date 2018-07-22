% sarcDetectMulti
%
% Code for automatic detection of sarcomere alignment in cardiomyocytes
% stained for sarcomeric alpha-actinin
%
% This code assumes you have already run sarcDetect and generated a
% *.sarcDetect.Settings.mat file for the set of images to be analyzed

% Arguments:  *.sarcDetect.Settings.mat file

%             Images of alpha-actinin stained sarcomeres captured from the 
%               same cover slip under the same exact aquisition conditions, 
%               file should TIF format grayscale at least 8-bit depth
%
% Returns:    Binary skeleton images of restored sarcomeres for each image
%             
%             Binary skeleton images of restored sarcomeres after ROI based
%               removal of non-sarcomeres for each image
%             
%             *.sarcOrientation.mat file containing sarcomere alignment
%
%             *.sarcOrientationALL.mat file containing summarized sarcomer
%               alginment for all images analyzed

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
   [file,path]=uigetfile({'*.mat';'*.*'},'Select Settings File...','D:/Manuscripts/2D Tissue Manuscripts/');
   filename = [path file];
       disp(filename)
   load (filename);

   % ask user to select file
    [file,path]=uigetfile({'*.TIF';'*.*'},'Select Images...',[path],'MultiSelect','on');   

    NumFiles = length(file);

    for i=1:NumFiles
        temp = file{i};
        FileIndex{i} = temp;
    end

    % Set variable to store SUM orientation angles for all images analyzed
    All_angles = 0;
    
    for count=1:NumFiles
       filename = FileIndex{count};
        path_and_filename = [path filename] 
        filename2 = path_and_filename

           im = imread(filename2);
        
           % Identify ridge-like regions and normalise image
           disp('Normalizing Image and Creating Mask' )
           [normim, mask] = ridgesegment(im, blksze, thresh);
           show(normim,1);
           show(mask, 2);
           normim_bin = normim > 0;
           mask_bin = mask > 0;
           normim_mask = normim_bin.*mask_bin;
           norm_bin_skel = bwmorph(normim_mask,'skel',Inf);
           show(norm_bin_skel, 10);

           % Determine ridge orientations
           disp('Calculating Ridge Orientations' )
           [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
           % Multiply orientation values by pi*3/2 for direction orthoganol to
           % sarcomere in the direction of force generation
           orientim_perp = orientim + (pi*3/2);
           plotridgeorient(orientim_perp, 25, im, 4)
           show(reliability,5)

           % Determine ridge frequency values across the image
           disp('Calculating Ridge Frequency Values' )
           [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 3, MinSarcSpacing, MaxSarcSpacing);

           % Actually I find the median frequency value used across the whole
           % fingerprint gives a more satisfactory result...
           freq = medfreq.*mask;

           newim = ridgefilter(normim, orientim, freq, kx, ky, 1);
           show(newim,6);
           binim = newim > 0;
           binim_skel = bwmorph(binim,'skel',Inf);
           binim_skel = bwmorph(binim_skel,'clean',Inf);
           show(binim_skel,11);
           
           % save filtered image of sarcomere skeleton before manual removal of
           % non-sarcomeres
           type = '.sarcomere.FULL.tif';
           filename3 = [filename2 type];
           imwrite(binim_skel,filename3,'Compression','none');

           % Remove false sarcomeres at tissue borders selecting regions to be
           % masked
           index = 0;
           while index < 1;
               disp('Select ROI to exclude from further analysis');
               BW = roipoly(binim_skel);
               BW2 = ~BW;
               binim_skel = binim_skel.*BW2;
               % Use erode mask to remove false sarcomeres from binary skeleton 
               show(binim_skel)
               b = input('Select another ROI to exclude? (yes = 0, no = 1): ');
               if b == 0
                   disp('Image accepted' )
               else
                   disp('Select another ROI...')
                   index = 1;
               end
           end

           % Multiply orientation angles by the binary skeleton image
           orientim_perp = orientim_perp.*binim_skel;
           % Convert 2D-array to 1D vector
           orientation = orientim_perp(:);
           % Keep non-zero values only
           nonzero_orientation = orientation(find(orientation));

           % Convert radians to degrees
           nonzero_orientation_angles = rad2deg(nonzero_orientation);

           % hist(nonzero_orientation_angles)
           Mean = mean(nonzero_orientation_angles)
           Std = std(nonzero_orientation_angles)
           Median = median(nonzero_orientation_angles)

           % Create histogram
           [n,xout] = hist(nonzero_orientation_angles,180);
           title([filename])
           dx = xout(2)-xout(1);                   % calc a single bin width
           n = n / sum( n*dx );                    % normalize histogram to have area of 1

           % Find mode
           [C,I] = max(n);
           Mode = xout(I)

           % Total number of sarcomere positive pixels in the skeleton image
           % Sarcomere density = total/(image area)
           Total = length(nonzero_orientation_angles)

           % save filtered image of sarcomere skeleton
           type = '.sarcomere.CLEAN.tif';
           filename3 = [filename2 type];
           imwrite(binim_skel,filename3,'Compression','none');
           
           % save filtered image of sarcomere skeleton
           type = '.sarcomere.CLEAN.tif';
           filename3 = [filename2 type];
           imwrite(binim_skel,filename3,'Compression','none');

           % save orientation angles
           ext_settings = '.sarcOrientation.mat';
           filename4 = [filename2 ext_settings];
           save(filename4,'orientim_perp','nonzero_orientation','nonzero_orientation_angles','Mode','Total'); 
           
           % Add data to SUM file
           All_angles = cat(1,All_angles,nonzero_orientation_angles);
           
           close all           
           
    end
    
   % Create SUM histogram
   [n,xout] = hist(All_angles,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1

   % Find mode
   [C,I] = max(n);
   Mode = xout(I)

   % Total number of sarcomere positive pixels in the skeleton image
   % Sarcomere density = total/(image area)
   Total = length(All_angles)/NumFiles
   
   % save SUM orientation angles
   ext_settings = '.sarcOrientationALL.mat';
   filename5 = [filename2 ext_settings];
   save(filename5,'All_angles','Mode','Total'); 
    
clear
