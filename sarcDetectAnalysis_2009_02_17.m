% sarcDetectAnalysis
%
% Code for automatic detection of sarcomere alignment in cardiomyocytes
% stained for sarcomeric alpha-actinin
%
% This code assumes you have already run sarcDetect and sarcDetectMulti and
% generated multiple .sarcOrientation.mat files for a set of images

% Arguments:  *.sarcOrientation.mat files
%
% Returns:    *.ALL.mat file containing summarized sarcomere alginment for 
%               all images analyzed
%
% Created by Adam W. Feinberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Havard University, Cambridge, MA 02138

% Last updated Feb 2009 by Adam W. Feinberg

   % ask user to select file
    [file,path]=uigetfile({'*.mat';'*.*'},'Select Orientation files...','D:/Manuscripts/2D Tissue Manuscripts/','MultiSelect','on');   

    NumFiles = length(file);

    for i=1:NumFiles
        temp = file{i};
        FileIndex{i} = temp;
    end

    All_rad=2*pi;
    
    for count=1:NumFiles
       filename = FileIndex{count};
        path_and_filename = [path filename] 
        filename2 = path_and_filename

           load(filename2);
           All_rad = cat(1,All_rad,nonzero_orientation);
           ModeAll(count) = Mode;
           TotalAll(count) = Total;
           figure
           hist(nonzero_orientation)
           title('filename')
           title(filename)           
    end
    
   % Create SUM histogram
   figure
   [n,xout] = hist(All_rad,180);
   bar(xout,n);
   title('Summary Histogram')
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1

   % Find mode
   [C,I] = max(n);
   Mode = xout(I)
   
   % Center data on MODE frequency
   shift = Mode - 2*pi;
   for i=1:1:length(All_rad)
       All_rad(i) = All_rad(i) - shift;
       if All_rad(i) > 5*pi/2
           All_rad(i) = All_rad(i) - pi;
       elseif All_rad(i) < 3*pi/2
           All_rad(i) = All_rad(i) + pi;
       end
   end

   [n,xout] = hist(All_rad,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1

   % Plot histogram of recentered data
   figure, bar(xout,n,'hist')              % plot normalized histogram
   title('Centered on the Mode')
   xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
      
   % Creat Rose histogram
   figure
   rose(All_rad, 360)
   title('Centered on the Mode')
   
   % Center data on MEAN frequency
   All_rad_mean = mean(All_rad);
   shift = All_rad_mean - 2*pi;
   for i=1:1:length(All_rad)
       All_rad(i) = All_rad(i) - shift;
       if All_rad(i) > 5*pi/2
           All_rad(i) = All_rad(i) - pi;
       elseif All_rad(i) < 3*pi/2
           All_rad(i) = All_rad(i) + pi;
       end
   end

   [n,xout] = hist(All_rad,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1

   % Plot histogram of recentered data
   figure, bar(xout,n,'hist')              % plot normalized histogram
   title('Centered on the Mean')
   xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
      
   % Creat Rose histogram
   figure
   rose(All_rad, 360)
   title('Centered on the Mean')
   
   % Find mode
   [C,I] = max(n);
   Mode = xout(I)
   
   sarcX = sum(cos(All_rad));
   sarcY = sum(sin(All_rad));
   maxRad = atan(sarcY/sarcX);
   
   % Center data on max Vector frequency
   shift = maxRad;
   for i=1:1:length(All_rad)
       All_rad(i) = All_rad(i) - shift;
       if All_rad(i) > 5*pi/2
           All_rad(i) = All_rad(i) - pi;
       elseif All_rad(i) < 3*pi/2
           All_rad(i) = All_rad(i) + pi;
       end
   end

   [n,xout] = hist(All_rad,180);
   dx = xout(2)-xout(1);                   % calc a single bin width
   n = n / sum( n*dx );                    % normalize histogram to have area of 1

   % Plot histogram of recentered data
   figure, bar(xout,n,'hist')              % plot normalized histogram
   title('Centered on the Total X Vector')
   xlim( [xout(1)-dx/2,xout(end)+dx/2] );  % make sure that the axis is squeezed to it's limits
      
   % Creat Rose histogram
   figure
   rose(All_rad, 360)
   title('Centered on the Total X Vector')

   sarcX = sum(cos(All_rad));
   sarcY = sum(sin(All_rad));
      
   % Total number of sarcomere positive pixels in the skeleton image
   % Sarcomere density = total/(image area)
   Total = length(All_rad)
   TotalSarcX = sarcX/Total*100
   TotalSarcY = sarcY/Total*100
   
   % save SUM orientation angles
   ext_settings = '.ALL.mat';
   filename5 = [filename2 ext_settings];
   save(filename5,'All_rad','Mode','Total','TotalAll','ModeAll'); 
    
clear
