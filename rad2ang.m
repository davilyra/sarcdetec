function degrees = rad2ang(radians)
% function angles = rad2ang(radians)
% Converts radians to degrees.
% Inputs:
%   - radians: a vector, matrix, or tensor containing the angles in radians
% Outputs:
%   - degrees: the input angles converted in degrees
%
% July 2014
% By Davi M. Lyra-Leite <davi@ieee.org>

degrees = (180/pi)*radians;