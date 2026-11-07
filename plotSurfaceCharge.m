
function plotSurfaceCharge(l_in,A_in,B_in, l_out,A_out,B_out, R, nsphere)
% plotSurfaceCharge  Plot surface charge density on sphere r=R
% l_in,A_in,B_in : vectors for interior multipoles (same length)
% l_out,A_out,B_out : vectors for exterior multipoles
% R : sphere radius (scalar)
% nsphere : sphere resolution (e.g. 80)

if nargin<8, nsphere = 80; end
if nargin<7, R = 1; end
eps0 = 8.854187817e-12;

% spherical grid
[phi,theta] = meshgrid(linspace(0,2*pi,nsphere), linspace(0,pi,nsphere));
mu = cos(theta);            % size nsphere x nsphere
sinth = sin(theta);

% helper to compute P_l(mu) for many l: returns (nsphere x nsphere x L)
function P = compute_Pl_array(lvec, mu)
    L = numel(lvec);
    sz = size(mu);
    P = zeros([sz, L]);
    mu_row = reshape(mu,1,[]);
    for k=1:L
        Pall = legendre(lvec(k), mu_row);      % (m+1) x Nmu, m=0 row is P_l
        P(:,:,k) = reshape(Pall(1,:), sz);
    end
end

% combine l vectors and compute Pl for all unique ls needed
l_all = unique([l_in(:); l_out(:)]);
Pl_all = compute_Pl_array(l_all, mu);

% function to get P_l for a given l using lookup in l_all
function P = getP(l)
    idx = find(l_all==l,1);
    P = Pl_all(:,:,idx);
end

% compute radial derivative of Vr at r=R for interior and exterior terms
% Vr = A r^l + B r^{-l-1}  => dVr/dr = l A r^{l-1} - (l+1) B r^{-l-2}
dVr_in = zeros(size(mu)); dVr_out = zeros(size(mu));
for k=1:numel(l_in)
    lcur = l_in(k);
    P_l = getP(lcur);
    A = A_in(k); B = B_in(k);
    dVr_in = dVr_in + ( lcur*A*R^(lcur-1) - (lcur+1)*B*R^(-lcur-2) ) .* P_l;
end
for k=1:numel(l_out)
    lcur = l_out(k);
    P_l = getP(lcur);
    A = A_out(k); B = B_out(k);
    dVr_out = dVr_out + ( lcur*A*R^(lcur-1) - (lcur+1)*B*R^(-lcur-2) ) .* P_l;
end

% radial E is -dVr/dr * P_l contribution already included above:
% note we summed dVr*P_l so E_r = - (that sum)
E_r_in_surf  = - dVr_in;
E_r_out_surf = - dVr_out;

% surface charge density
sigma = eps0 * (E_r_out_surf - E_r_in_surf);  % nsphere x nsphere

% build sphere coordinates
x = R .* sinth .* cos(phi);
y = R .* sinth .* sin(phi);
z = R .* cos(theta);

% plot colored sphere: blue positive, red negative (bwr)
figure('Color','w');
s = surf(x,y,z, sigma, 'FaceColor','interp', 'EdgeColor','none', 'FaceAlpha',1);
axis equal; xlabel('x'); ylabel('y'); zlabel('z');
camlight head;
lighting gouraud;


% create symmetric colormap red->white->blue
n = 256;
half = floor(n/2);
reds  = [linspace(1,1,half)' linspace(0,1,half)' linspace(0,1,half)'];   % red -> white
blues = [linspace(1,0,n-half)' linspace(1,1,n-half)' linspace(1,1,n-half)']; % white -> blue
cmap = [reds; blues];

colormap(cmap);
% center color scale at zero
cmax = max(abs(sigma(:)));
clim([-cmax cmax]);
colorbar;

title('Surface charge density \sigma (blue>0, red<0)');
end