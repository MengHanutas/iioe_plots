function [var_q] = ra_projvar(var, x, xq)

% This function projects both physical and biogeochemical variables on to
% isopycnal surfaces
% ALGORITHM: Projecting all the profiles to given density range and then
% averaging all the profiles along the isopycnal surfaces
% NOTE: Nutrient data must be in dimension as levels x stations
%           rho data must be in same dimension as nutrients
%

% var = ctdm; xq = neutralgrid; x = gm_ndm;

stn = size(var, 2);

var_q = NaN(length(xq), stn);

for prof = 1:stn
    z = x(:, prof);
    para = var(:, prof);
    xx = find(~isnan(z .* para));
    % Reference nutrient to nutrientr(sigmarho)
    if ~isempty(xx)
        varz = interp1(z(xx), para(xx), xq);
        var_q(:, prof) = varz;
    else
        continue
    end%endif
end%endfor
