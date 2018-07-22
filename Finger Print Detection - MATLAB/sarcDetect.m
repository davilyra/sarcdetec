% sarcDetect
% 
% Code for automatic detection of sarcomer alignment in cardiomyocytes
% stained for sarcomeric alpha-actinin
%
% Adapted from Function to demonstrate use of fingerprint code
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
% January 2005
%
% Last updated May 2007 by Adam W. Feinberg

    % ask user to select file
    [file,path]=uigetfile({'*.TIF';'*.bmp';'*.jpg';'*.*'},'Select Image File...','F:/Immunofluorescent data/');   
    filename = [path file];
	im = imread(filename);
    
    % Identify ridge-like regions and normalise image
    index = 0;
    while index < 1;
        % Input blocksize over which the standard deviation is determined
        blksze = input('Enter blocksize (16): ');
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
    
    % Determine ridge orientations
%     while index < 1;
%         % Sigma of the derivative of Gaussian used to compute image gradients.
%         gradientsigma = input('Enter blocksize (16): ');
%         % Threshold of standard deviation to decide if a block is a ridge region
%         thresh = input('Enter Threshold (0.1 - 0.2): ');
%         [normim, mask] = ridgesegment(im, blksze, thresh);
%         show(normim,1);
%         show(mask, 2);
%         w = waitforbuttonpress;
%         if w == 0
%             disp('Image threshold accepted' )
%             index = 1;
%         else
%             disp('Re-analyze imaging...')
%         end
%     end
    disp('Calculating Ridge Orientations' )
    [orientim, reliability] = ridgeorient(normim, 1, 3, 3);
    orientim_perp = orientim + (pi/2);
    plotridgeorient(orientim_perp, 5, im, 3)

    show(reliability,5)
    type = '_orientim.txt';
    savefile1 = [filename type];
    save savefile1 orientim -ASCII
    
    % Determine ridge frequency values across the image
    blksze = 8;
    disp('Calculating Ridge Frequency Values' )
    [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 3, 2, 9);
%    show(freq,3) 
    
    % Actually I find the median frequency value used across the whole
    % fingerprint gives a more satisfactory result...
    freq = medfreq.*mask;
    
    % Now apply filters to enhance the ridge pattern
    Size = 4;
    Min = 0.3;
    Max = 0.6;
    cell = 1;
    for i=Min:0.1:Max
        for j=Min:0.1:Max
            disp('Applying Filter Step')
            disp(cell)
            newim = ridgefilter(normim, orientim, freq, i, j, 1);
            binim = newim > 0;
            binim_skel = bwmorph(binim,'skel',Inf);
            binim_skel = bwmorph(binim_skel,'clean',Inf);
            kx_i = num2str(i);
            ky_j = num2str(j);
            comma = ', ';
            title_cell = [kx_i comma ky_j];
            imshow(binim_skel)
            subplot(4,4,cell), imshow(binim_skel), title(title_cell)
            cell = cell + 1;
        end
    end

    cell = 1;
    for i=Min:0.1:Max
        for j=Min:0.1:Max
            kx_i = num2str(i);
            ky_j = num2str(j);
            comma = ', ';
            title_cell = [kx_i comma ky_j];
            subplot(4,4,cell), title(title_cell)
            cell = cell + 1;
        end
    end  
    
    index = 0;
    while index < 1;
        % kx controls the sigma in the x direction which is along the
        % filter, and hence controls the bandwidth of the filter.  
        kx = input('Enter Sigma along filter (0.5): ');
        % ky controls the sigma across the filter and hence controls the
        % orientational selectivity of the filter. A value of 0.5 for both
        % kx and ky is a good starting point. 
        ky = input('Enter Sigma across filter (0.5): ');
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
 
    orientation = orientim(:);
    nonzero_orientation = orientation(find(orientation));
    nonzero_orientation_angles = rad2deg(nonzero_orientation + pi/2);
    
    % hist(nonzero_orientation_angles)
    Mean = mean(nonzero_orientation_angles)
    Std = std(nonzero_orientation_angles)
    Median = median(nonzero_orientation_angles)
     
    [n,xout] = hist(nonzero_orientation_angles,180);
    [C,I] = max(n);
    Mode = xout(I)
    figure, hist(nonzero_orientation_angles,180);
    [u,sig,t,iter] = fit_mix_gaussian( n,M )

    Total = length(nonzero_orientation_angles)
    
    type = '_skel.tif';
    filename2 = [filename type];
    imwrite(binim_skel,filename2,'Compression','none');
