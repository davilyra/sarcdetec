%{
orrientationalOrder_Matrix

Last updated: 07/25/2014

written for use with sarcDetect
Created by: Anya Grosberg
Disease Biophysics Group
School of Engineering and Applied Sciences
Havard University, Cambridge, MA 02138

Updated by: Davi M. Lyra Leite <davi@ieee.org>
Laboratory for Living Systems Engineering
Viterbi School of Engineering
University of Southern California, Los Angeles, CA 90089


The purpose of this code is to calculate the orientational order parameter 
for the detected sarcomeres.

Input: the code asks for the .mat file produced by the sarcDetect code.
        it will use the nonzero_orientation variable
Output: director - this is the major orientation direction of the sarcomeres
                    and can be used as a check.
        orientational order parameter - a number between 0 (isotropic)
                                        and 1 (perfect alignment).

%}

%% ========================================================================
% Clearing workspace
clear all
close all
% clc

%% ========================================================================
% ask user to select file produced by the sarcDetect code
[file,path]=uigetfile({'*.mat';'*.*'},'Select Settings File...');
filename = [path file];
disp(filename)
disp(' ')
load(filename);


%% ========================================================================
% Take the non-zero angles -- this is important, because if the zero
% angles are not removed the orientational order parameter is very high
% as the zero vectors drown out the information from the other vectors.
% angles = all_images;
angles = deg2rad(all_images); %this is for use with Megan McCains generalized variables.
% angles = nonzero_orientation;


%% ========================================================================
% Tensor Method
% Calculate x and y components of each vector r
r(1,:) = cos(angles);
r(2,:) = sin(angles);

% Calculate the Orientational Order Tensor for each r and the average
% Orientational Order Tensor (OOT_Mean)
OOT_All = zeros(2,2,length(r));
for i=1:2
    for j=1:2
        OOT_All(i,j,:)=r(i,:).*r(j,:);
    end
end
OOT_Mean = mean(OOT_All,3);


%% ========================================================================
% Processing the Orientational Order Tensor (OOT)
% Normalize the orientational Order Tensor, this is necessary to get the
% order paramter in the range from 0 to 1
OOT = 2.*OOT_Mean - eye(2);

% Find the eigenvalues (orientational parameters) and
% eigenvectors (directions) of the OOT
[directions,orient_parameters] = eig(OOT);

% orientational order parameters is the maximum eigenvalue, while the
% director is the corresponding eigenvector
[Orientation_order_parameter,I] = max(max(orient_parameters));
director = directions(:,I);


% calculate the angle corresponding to the director, note that by symmetry
% the director = - director. This implies that any angle with a period of
% 180 degrees will match this director. To help compare these results to
% the plot results we enforce the period to match the period of the
% original data.
directionAngle_defualt = acosd(director(1)/sum(director.^2));
directionAngle = directionAngle_defualt + 180*(floor(min(angles)/pi()));


%% ========================================================================
% Calculate the difference between the director and the mean of the
% angles. Note, that these are not necessarily the same thing because we
% have a finite number of vectors, so there is some inacuracy introduced
% in both methods. We can expect the difference to be very large for
% isotropic and small for well aligned structures. The output of this is
% suppressed unless someone needs it for something.
direction_error = directionAngle - (180/pi())*mean(angles);

% save orientational order parameters, director, and the direction Angle
ext_settings = '.OrientationOrderParameter.mat';
filename2 = [filename ext_settings];
save(filename2,'Orientation_order_parameter','director','directionAngle');

% Display the values
display(['Orientational Order Parameter = ',...
    num2str(Orientation_order_parameter)])
display(['Direction Angle = ',...
    num2str(directionAngle)])
