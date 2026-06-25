function point_xyz = ray_intersect(F, origin, direction, t_min, t_max,varargin)
% Finds the first point of intersection between a ray emanating from the
% specified "origin" point, propagating along "direction", and a surface
% characterized by a scalar field "F" set equal to zero. Searches for a
% point of intersection between "t_min" and "t_max" for this point of
% intersection.
% 
% Inputs:
% - F should be provided as a scalar field equal to zero. For a circle with
% radius 1, F should be provided as the following function handle: F =
% @(x) x(1).^2 + x(2).^2 + x(3).^2 - 1.
% - Origin and direction are three-component vectors specifying positions
% and directions, respectively.
% - t_min and t_max are scalars for the search distance
% - varargin is for later use.

n_t_steps = 1e5;

% Parametrize the ray
origin = origin(:)'; direction = direction(:)';
ray = @(t) origin + t.*direction;

% Find the values of F along the ray direction
F_path = @(t) F(ray(t));
t = linspace(t_min,t_max,n_t_steps);
F_values = arrayfun(F_path,t); % Now F is a series of function values

% Find intersection points
idxs = find(diff(sign(F_values)) ~= 0); % indices of intersection in t-space
if isempty(idxs)
    error('No intersection found. Consider a different t range, finer t-stepping, and the ray direction.')
end
t_guess = t(idxs(1)); % We know that the first (and closest) intersection point is between t(idxs(1)) and t(idxs(2))
F_opt = @(t) (F_path(t)).^2;
point_t = fminsearch(F_opt,t_guess);
point_xyz = ray(point_t);
end