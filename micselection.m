%% Microscope selection routine
% This is script was written to be added as a routine to identify the type
% of microscope being used, the objective and make the process of defining
% the corresponding pixel to micrometer conversion easier.

% By Davi M. Lyra Leite <davi@ieee.org>
% Last modified: July 2014


%% Select the microscope manufacturer
disp('Select your microscope manufacturer: ')
mictype = input('1 - Nikon\n2 - Nikon Confocal\n3 - Zeiss:\n');
disp(' ')

% Obtaining the calibration factor according to the microscope brand and 
% the objective magnification
disp('Select your objective:')
switch mictype    
    % Nikon
    case 1
        magnification = input('1 - 10x;\n2 - 20x;\n3 - 60x:\n');
        if magnification == 1
            pix2um = 1.5385;
        elseif magnification == 2
           pix2um = 3.0303;
        elseif magnification == 3
            pix2um = 9.0909;
        end
        
        % Nikon
    case 2
        disp('1x:')
        disp('1 - 60x: 512x512;')
        disp('2 - 60x: 1024x1024;')
        disp('3 - 60x: 2048x2048;')
        disp('2x:')
        disp('4 - 60x: 512x512;')
        magnification = input('5 - 60x: 1024x1024:\n');
        if magnification == 1
            pix2um = 2.4136;
        elseif magnification == 2
           pix2um = 4.8272;
        elseif magnification == 3
            pix2um = 1/0.1035801;
        elseif magnification == 4
           pix2um = 4.8272;
        elseif magnification == 5
           pix2um = 9.6544;
        end
        
	% Zeiss
    case 3
        magnification = input('1 - 20x;\n2 - 40x;\n3 - 63x:\n');
        if magnification == 1
            pix2um = 3.085;
        elseif magnification == 2
           pix2um = 6.22;
        elseif magnification == 3
            pix2um = 9.8;
        end
end

% convert 2 um sarcomere spacing into pixels
SarcSpacing = 2*pix2um;

% minimum length of sarcomere spacing
MinSarcSpacing = round(0.5*SarcSpacing);

% maximum length of sarcomere spacing
MaxSarcSpacing = round(1.5*SarcSpacing);
blksze = MaxSarcSpacing;

