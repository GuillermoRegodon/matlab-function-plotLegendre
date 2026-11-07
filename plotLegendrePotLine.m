
function plotLegendrePotLine(l_in,A_in,B_in, l_out,A_out,B_out, r_plane, R, ngrid)
% plotLegendrePot  Plot Legendre-series potential only on the circle r = R
% in the x-z plane.

if nargin < 9 || isempty(ngrid), ngrid = 300; end
if nargin < 8 || isempty(R), R = 1; end
if nargin < 7 || isempty(r_plane), r_plane = 2*R; end %#ok<NASGU>

% Angle around the circle
theta = linspace(0, 2*pi, ngrid);
x = R * sin(theta);
z = R * cos(theta);
cos_theta = cos(theta);

Phi = zeros(size(theta));

% Inside-series contribution evaluated on r = R
for k = 1:numel(l_in)
    l = l_in(k);
    Pl = legendre(l, cos_theta);
    Pl0 = squeeze(Pl(1,:,:));

    a = 0;
    b = 0;
    if ~isempty(A_in), a = A_in(k); end
    if ~isempty(B_in), b = B_in(k); end

    Phi = Phi + (a * R^l + b * R^(-l-1)) .* Pl0;
end

% Outside-series contribution evaluated on r = R
for k = 1:numel(l_out)
    l = l_out(k);
    Pl = legendre(l, cos_theta);
    Pl0 = squeeze(Pl(1,:,:));

    a = 0;
    b = 0;
    if ~isempty(A_out), a = A_out(k); end
    if ~isempty(B_out), b = B_out(k); end

    Phi = Phi + (a * R^l + b * R^(-l-1)) .* Pl0;
end

% Plot only the circle
%figure
scatter(x, z, 20, Phi, 'filled')
hold on
plot(x, z, 'k-', 'LineWidth', 1)
hold off
axis equal tight
xlabel('x')
ylabel('z')
title('Electric potential on r = R in the x-z plane')

% Custom diverging colormap: red for negative, blue for positive with white center
ncol = 256;
neg = [linspace(1,1,ncol/2)' linspace(0,1,ncol/2)' linspace(0,1,ncol/2)']; % red->white
pos = [linspace(1,0,ncol/2)' linspace(1,0,ncol/2)' linspace(1,1,ncol/2)']; % white->blue
cmap = [neg; pos];
colormap(cmap);
colorbar;

cmax = max(abs(Phi));
if cmax == 0
    cmax = 1;
end
caxis([-cmax cmax])

end