% sarcDetect Area calculation
%
% This code takes the binary skeleton images of the sarcomeres generated by
% sarcDetect and determines the Area covered by sarcomeres based on the full
% detection image ("FULL") and the manual segmented image to remove
% non-sarcomeres ("CLEAN")
%
% Adam Feinberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University
% feinberg@seas.harvard.edu
% http://www.csse.uwa.edu.au/~pk
% Ferbuary 2009
%
% Last updated Feb 2009 by Adam W. Feinberg

   % ask user to select image files
    [file,path]=uigetfile({'*.TIF';'*.*'},'Select Images...','D:/Manuscripts/2D Tissue Manuscripts/','MultiSelect','on');   

    NumFiles = length(file);

    for i=1:NumFiles
        temp = file{i};
        FileIndex{i} = temp;
    end

   % Enter pixel to micrometer converesion factor to define appropriate
   % sarcomere spacing
   pix2um = input('Enter calibration factor in pixels per micrometer (Coolsnap 63x=9.8, 40x=6.22, 20x=3.085): ');
   SarcSpacing = 2*pix2um;          % convert 2 um sarcomere spacing into pixels
   MinSarcSpacing = round(0.5*SarcSpacing); % minimum length of sarcomere spacing
   MaxSarcSpacing = round(1.5*SarcSpacing); % maximum length of sarcomere spacing
   
   % Number of times to dilate and then erode
   dilate_steps = MaxSarcSpacing/2;
   
   index=1
    
    for count=1:2:NumFiles
       % load "CLEAN" image file
       filename = FileIndex{count};
       path_and_filename = [path filename] 
       filename2 = path_and_filename
       imClean = imread(filename2);
       
       % load "FULL" image file
       filename3 = FileIndex{count+1};
       path_and_filename2 = [path filename3] 
       filename4 = path_and_filename2
       imFull = imread(filename4);
        
       % Merge sarcomers for "FULL" image
       show(imFull)
           for i=1:dilate_steps
               imFull=bwmorph(imFull,'dilate');
           end

           for i=1:dilate_steps
               imFull=bwmorph(imFull,'erode');
           end
       show(imFull)
       
       % Merge sarcomeres for "CLEAN" image
       show(imClean)
           for i=1:dilate_steps
               imClean=bwmorph(imClean,'dilate');
           end

           for i=1:dilate_steps
               imClean=bwmorph(imClean,'erode');
           end
       show(imClean)
       
       % Save images of sarcomere area
       type = '.AREA.tif';
       filename5 = [filename2 type];
       filename6 = [filename4 type];
       imwrite(imClean,filename5,'Compression','none');
       imwrite(imFull,filename6,'Compression','none');
       
       % Clean Area
       AreaTotal = imClean(:);
       % find non zero pixels
       AreaNonZero = AreaTotal(find(AreaTotal));
       % Convert area in pixels to micrometers
       sarcArea(count)=sum(AreaNonZero)/(pix2um^2);
       
       % Full Area
       AreaTotal = imFull(:);
       % find non zero pixels
       AreaNonZero = AreaTotal(find(AreaTotal));
       % Convert area in pixels to micrometers
       sarcArea(count+1)=sum(AreaNonZero)/(pix2um^2);
       
       sarcAreaClean(index)=sarcArea(count);
       sarcAreaFull(index)=sarcArea(count+1);
       sarcAreaCleanPercent(index)=sarcArea(count)/(length(AreaTotal)/(pix2um^2))*100;
       sarcAreaFullPercent(index)=sarcArea(count+1)/(length(AreaTotal)/(pix2um^2))*100;
       index=index+1;
       
       
    end

   AreaTotal_um= sum(AreaTotal)/(pix2um^2);

   MeanAreaClean = mean(sarcAreaClean)
   StdAreaClean = std(sarcAreaClean)
   MeanAreaFull = mean(sarcAreaFull)
   StdAreaFull = std(sarcAreaFull)
   MeanPercentAreaClean = mean(sarcAreaCleanPercent)
   StdPercentAreaClean = std(sarcAreaCleanPercent)
   MeanPercentAreaFull = mean(sarcAreaFullPercent)
   StdPercentAreaFull = std(sarcAreaFullPercent)   

   % save data as .MAT file
   ext_settings = '.AREA_DATA.mat';
   filename3 = [filename2 ext_settings];
   save(filename3,'AreaTotal_um','sarcArea','sarcAreaClean','sarcAreaCleanPercent','sarcAreaFull','sarcAreaFullPercent'); 
   
   clear
