% TESTFIN 
%
% Function to demonstrate use of fingerprint code
%
% Usage:  [newim, binim, mask, reliability] =  testfin(im);
%
% Argument:   im -  Fingerprint image to be enhanced.
%
% Returns:    newim - Ridge enhanced image.
%             binim - Binary version of enhanced image.
%             mask  - Ridge-like regions of the image
%             reliability - 'Reliability' of orientation data

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005

    % ask user to select file
    [file,path]=uigetfile({'*.tif';'*.bmp';'*.jpg';'*.*'},'Select Image File...','D:/Parker Lab/');   
    filename = [path file];
	im = imread(filename);
    
    % Identify ridge-like regions and normalise image
    blksze = 16; thresh = 0.2;
    [normim, mask] = ridgesegment(im, blksze, thresh);
    show(normim,1);
    
    % Determine ridge orientations
    [orientim, reliability] = ridgeorient(normim, 1, 5, 5);
    orientim_perp = orientim + (pi/2);
    plotridgeorient(orientim, 20, im, 2)
    plotridgeorient(orientim_perp, 20, im, 3)
 
    show(reliability,6)
    type = '_orientim.txt';
    savefile1 = [filename type];
    save savefile1 orientim -ASCII
    
    % Determine ridge frequency values across the image
    blksze = 36; 
    [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 5, 5, 15);
%    show(freq,3) 
    
    % Actually I find the median frequency value used across the whole
    % fingerprint gives a more satisfactory result...
    freq = medfreq.*mask;
    
    % Now apply filters to enhance the ridge pattern
    newim = ridgefilter(normim, orientim, freq, 0.5, 0.5, 1);
    show(newim,4);
    
    % Binarise, ridge/valley threshold is 0
    binim = newim > 0;
    show(binim,5);

    % Display binary image for where the mask values are one and where
    % the orientation reliability is greater than 0.5
    show(binim.*mask.*(reliability>0.5), 7)
    orientim_perp_mask = orientim_perp.*binim;
    plotridgeorient(orientim_perp_mask, 5, im, 8)
    masked_orientation = orientim_perp.*binim;
    plotridgeorient(masked_orientation, 10, orientim_perp_mask, 9)
    
    orientation = orientim(:);
    nonzero_orientation = orientation(find(orientation));
    nonzero_orientation_angles = rad2deg(nonzero_orientation + pi/2);
    
    % hist(nonzero_orientation_angles)
    Mean = mean(nonzero_orientation_angles)
    Std = std(nonzero_orientation_angles)
    Mean_text = num2str(Mean);
    Std_text = num2str(Std);
    disp(Mean_text)
    disp(Std_text)
   
