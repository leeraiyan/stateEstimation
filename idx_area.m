function [AREA_I, PRICE_REF_BUS] = idx_area
%IDX_AREA   Defines constants for named column indices to areas matrix.
%   [AREA_I, PRICE_REF_BUS] = idx_area
%
%   The index, name and meaning of each column of the areas matrix is given
%   below:
% 
%   columns 1-2
%    1  AREA_I       	area number
%    2  PRICE_REF_BUS	price reference bus for this area


%% define the indices
AREA_I          = 1;    %% area number
PRICE_REF_BUS   = 2;    %% price reference bus for this area

return;
