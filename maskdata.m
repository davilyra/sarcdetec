function masked = maskdata(normim,mask)
% function masked = maskdata(normim,mask)
% Masks the data.
% Inputs:
%   - normim: normalized image
%   - mask: mask
% Outputs:
%   - masked: the masked image
%
% July 2014
% By Davi M. Lyra-Leite <davi@ieee.org>

normim_bin = normim > 0;
mask_bin = mask > 0;
normim_mask = normim_bin.*mask_bin;
masked = bwmorph(normim_mask,'skel',Inf);