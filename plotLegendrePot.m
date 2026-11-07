
function plotLegendrePot(l_in,A_in,B_in, l_out,A_out,B_out, r_plane, R, ngrid)
% plotLegendrePot  Plot Legendre-series electric potential in x-z plane.
% Arguments:
%  l_in, A_in, B_in   : vectors of same length (may be empty) for terms valid inside r<=R
%  l_out,A_out,B_out  : vectors of same length (may be empty) for terms valid outside r>=R
%  r_plane            : half-width of plotting window (plot from -r_plane to +r_plane)
%  R                  : sphere radius (matching continuity region)
%  ngrid              : number of grid points per axis
%
% The potential is assumed axisymmetric (m = 0):
%  Phi(r,theta) = sum_l [ A_l * r^l + B_l * r^(-l-1) ] * P_l(cos(theta))
% with the specified terms applied only inside or outside the sphere.

% --- Input checks
if nargin<9, ngrid = 300; end
if nargin<8, R = 1; end
if nargin<7, r_plane = 2*R; end

if isempty(l_in), A_in = []; B_in = []; end
if isempty(l_out), A_out = []; B_out = []; end

% --- Build grid in x-z plane (y=0)
x = linspace(-r_plane,r_plane,ngrid);
z = linspace(-r_plane,r_plane,ngrid);
[X,Z] = meshgrid(x,z);
Y = zeros(size(X));
r = sqrt(X.^2 + Y.^2 + Z.^2);
% avoid singularity at r=0 for negative powers
r_safe = r;
r_safe(r_safe==0) = eps;

% colatitude theta (0 at +z axis)
cos_theta = zeros(size(r));
inside_mask = (r <= R);
outside_mask = (r >= R);
% define cos(theta) where r>0; at r=0 set theta=0 (cos_theta = 1)
cos_theta(r>0) = Z(r>0) ./ r(r>0);
cos_theta(r==0) = 1;

Phi = zeros(size(r));

% --- helper to evaluate P_l(cos_theta) (m=0)
for k=1:numel(l_in)
    l = l_in(k);
    Pl = legendre(l, cos_theta);      % size (m+1) x size(cos_theta)
    Plm0 = squeeze(Pl(1,:,:));        % m=0 row
    a = 0; b = 0;
    if ~isempty(A_in), a = A_in(k); end
    if ~isempty(B_in), b = B_in(k); end
    term = ( a .* (r_safe.^l) + b .* (r_safe.^(-l-1)) ) .* Plm0;
    term(~inside_mask) = 0;            % apply only inside
    Phi = Phi + term;
end

for k=1:numel(l_out)
    l = l_out(k);
    Pl = legendre(l, cos_theta);
    Plm0 = squeeze(Pl(1,:,:));
    a = 0; b = 0;
    if ~isempty(A_out), a = A_out(k); end
    if ~isempty(B_out), b = B_out(k); end
    term = ( a .* (r_safe.^l) + b .* (r_safe.^(-l-1)) ) .* Plm0;
    term(~outside_mask) = 0;           % apply only outside
    Phi = Phi + term;
end

% --- Plot
figure;
h = pcolor(X,Z,Phi);
set(h,'EdgeColor','none');
axis equal tight;
xlabel('x'); ylabel('z');
title('Electric potential (x-z plane)');

% Custom diverging colormap: red for negative, blue for positive with white center
ncol = 256;
neg = [linspace(1,1,ncol/2)' linspace(0,1,ncol/2)' linspace(0,1,ncol/2)']; % red->white
pos = [linspace(1,0,ncol/2)' linspace(1,0,ncol/2)' linspace(1,1,ncol/2)']; % white->blue
cmap = [neg; pos];
colormap(cmap);
colorbar;

% draw sphere cross-section (circle of radius R)
hold on;
theta_c = linspace(0,2*pi,360);
plot(R*cos(theta_c), R*sin(theta_c), 'k-', 'LineWidth',1.2);
hold off;

% center color at zero
cmax = max(abs(Phi(:)));
if cmax==0, cmax = 1; end
caxis([-cmax cmax]);

end