%{
orrientationalOrder_Matrix

Last updated: 02/24/2009

written for use with sarcDetect
Created by: Anya Grosberg
Disease Biophysics Group
School of Engineering and Applied Sciences
Havard University, Cambridge, MA 02138

The purpose of this code is to calculate the orientational order parameter 
for the detected sarcomeres.

Input: the code asks for the .mat file produced by the sarcDetect code.
        it will use the nonzero_orientation variable
Output: director - this is the major orientation direction of the sarcomeres
                    and can be used as a check.
        orientational order parameter - a number between 0 (isotropic)
                                        and 1 (perfect alignment).

%}
% ask user to select file produced by the sarcDetect code
   [file,path]=uigetfile({'*.mat';'*.*'},...
       'Select Settings File...');
   filename = [path file];
   disp(filename)
   load(filename);
   
   %Take the non-zero angles -- this is important, because if the zero
   %angles are not removed the orientational order parameter is very high
   %as the zero vectors drown out the information from the other vectors.
   angles = ang2rad(all_images); %change this if your angles are called something else
   %angles must be in radians
   %  
   %%%%%%%%%%%%% Tensor Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Calculate x and y components of each vector r
   r(1,:) = cos(angles);
   r(2,:) = sin(angles);
   %Calculate the Orientational Order Tensor for each r and the average
   %Orientational Order Tensor (OOT_Mean)
   for i=1:2
       for j=1:2
           OOT_All(i,j,:)=r(i,:).*r(j,:);
           %The mean of a matrix is the same as the mean of each element so
           %we calculate the OOT for each element separately.
           OOT_Mean(i,j) = mean(OOT_All(i,j,:));
       end
   end
   %Normalize the orientational Order Tensor (OOT), this is necessary to get the
   %order paramter in the range from 0 to 1
   OOT = 2.*OOT_Mean-eye(2);
   %Find the eigenvalues (orientational parameters) and 
   %eigenvectors (directions) of the Orientational Order Tensor
   [directions,orient_parameters]=eig(OOT);
   %orientational order parameters is the maximal eigenvalue, while the
   %direcotor is the corresponding eigenvector
   [Orientation_order_parameter,I] = max(max(orient_parameters));
   director = directions(:,I);
   %calculate the angle corresponding to the director, note that by symmetry
   %the director = - director. This implies that any angle with a period of
   %180 degrees will match this director. To help compare these results to
   %the plot results we enforce the period to match the period of the
   %original data.
   directionAngle_defualt = acosd(director(1)/sum(director.^2));
   directionAngle = directionAngle_defualt+180*(floor(min(angles)/pi()));
   %
   %Calculate the difference between the director and the mean of the
   %angles. Note, that these are not necessarily the same thing because we
   %have a finite number of vectors, so there is some inacuracy introduced
   %in both methods. We can expect the difference to be very large for
   %isotropic and small for well aligned structures. The output of this is
   %suppressed unless someone needs it for something.
   direction_error = directionAngle-(180/pi())*mean(angles);
    
   %save orientational order parameters, director, and the direction Angle
   ext_settings = '.OrientationOrderParameter.mat';
   filename2 = [filename ext_settings];
   save(filename2,'Orientation_order_parameter','director','directionAngle');
   
   %Display the values
   display([Orientation_order_parameter,directionAngle])
   
   clear