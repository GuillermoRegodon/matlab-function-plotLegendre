
% function signature (add arrowScale as last optional arg)
function plotLegendreField(l_in,A_in,B_in, l_out,A_out,B_out, r_plane, R, ngrid, arrowScale)
% plotLegendreField  Plot E-field arrows in x-z plane from Legendre-term potentials.
%
% Signature:
%   plotLegendreField(l_in,A_in,B_in, l_out,A_out,B_out, r_plane, R, ngrid)
% Inputs:
%   l_in,A_in,B_in   - vectors (or empty []) for interior terms (valid for r<R)
%   l_out,A_out,B_out- vectors (or empty []) for exterior terms (valid for r>R)
%   r_plane          - scalar half-width for x and z axes (plot region = [-r_plane,r_plane])
%   R                - sphere radius (for reference)
%   ngrid            - resolution (default 60)
%   arrowscale       - number to scale size of arrows, [] for the same size
%   
% Each term: V_l = (A*r^l + B*r^(-l-1)) * P_l(cosθ). Superposition applied.

if nargin<10, arrowScale = 0; end
if nargin<9, ngrid = 60; end
if nargin<8, R = 1; end
if nargin<7, r_plane = 2*R; end

% ensure column vectors
l_in = l_in(:); A_in=A_in(:); B_in=B_in(:);
l_out=l_out(:); A_out=A_out(:); B_out=B_out(:);

% grid in x-z plane (y=0)
xv = linspace(-r_plane,r_plane,ngrid);
zv = linspace(-r_plane,r_plane,ngrid);
[X,Z] = meshgrid(xv,zv);
Y = zeros(size(X));
r = sqrt(X.^2 + Y.^2 + Z.^2);
theta = acos( min(max( Z./(r + (r==0)), -1), 1) ); % θ in [0,pi], avoid div0
mu = cos(theta);
sinth = sin(theta);


r_safe = max(r, eps);
mask_in  = (r < R);    % inside sphere (strictly inside)
mask_out = (r >= R);   % outside sphere (including surface)

% initialize fields
Ex = zeros(size(X)); Ez = zeros(size(X));

% helper to compute Pl(mu) and P'l(mu)
    function [Pl, dPl] = Pl_and_dPl(lvec, mu_mat)
        % returns Pl summed over all l in lvec weighted later; but this helper returns arrays per l
        L = numel(lvec);
        Pl = zeros([size(mu_mat) L]);
        dPl = zeros([size(mu_mat) L]);
        for k=1:L
            lcur = lvec(k);
            % use legendre to get P_l and P_{l-1}
            Pall = legendre(lcur, mu_mat(:).'); % (m+1) x N
            Pl_k = reshape(Pall(1,:), size(mu_mat)); % m=0
            if lcur==0
                dPk = zeros(size(mu_mat));
            else
                Pallm1 = legendre(lcur-1, mu_mat(:).');
                P_lm1 = reshape(Pallm1(1,:), size(mu_mat));
                % derivative w.r.t mu
                dP_dmu = lcur./(mu_mat.^2 - 1) .* (mu_mat .* Pl_k - P_lm1);
                dPk = dP_dmu;
            end
            Pl(:,:,k) = Pl_k;
            dPl(:,:,k) = dPk;
        end
    end


% before loops compute:
r_safe = max(r, eps);
r_xy = sqrt(X.^2 + Y.^2);
cosphi = zeros(size(r)); sinphi = zeros(size(r));
nz = r_xy > eps;
cosphi(nz) = X(nz) ./ r_xy(nz);
sinphi(nz) = Y(nz) ./ r_xy(nz);
sinth = sqrt(max(0, 1 - mu.^2));
cos_theta = mu;


% process interior terms (valid for r>0 and typically r<R)
if ~isempty(l_in)
    [Pl_in, dPl_in] = Pl_and_dPl(l_in, mu);
    for k=1:numel(l_in)
        lcur = l_in(k); A = A_in(k); B = B_in(k);

        % inside-loop for each l_in,k:
        Vr_in = A_in(k).*r.^lcur + B_in(k).*r.^(-lcur-1);          % full grid
        dVdr_in = lcur.*A_in(k).*r.^(lcur-1) - (lcur+1).*B_in(k).*r.^(-lcur-2);
        
        E_r_in  = - dVdr_in .* Pl_in(:,:,k);
        E_th_in =  (1./r_safe) .* Vr_in .* sinth .* dPl_in(:,:,k);
        
        % apply mask so only inside points contribute
        E_r_in(~mask_in)  = 0;
        E_th_in(~mask_in) = 0;
        
        % convert and accumulate
        Ex = Ex + (E_r_in .* (sinth .* cosphi) + E_th_in .* (cos_theta .* cosphi));
        Ez = Ez + (E_r_in .* cos_theta - E_th_in .* sinth);

    end
end

% process exterior terms
if ~isempty(l_out)
    [Pl_out, dPl_out] = Pl_and_dPl(l_out, mu);
    for k=1:numel(l_out)
        lcur = l_out(k); A = A_out(k); B = B_out(k);

        Vr_out = A_out(k).*r.^lcur + B_out(k).*r.^(-lcur-1);
        dVdr_out = lcur.*A_out(k).*r.^(lcur-1) - (lcur+1).*B_out(k).*r.^(-lcur-2);
        
        E_r_out  = - dVdr_out .* Pl_out(:,:,k);
        E_th_out =  (1./r_safe) .* Vr_out .* sinth .* dPl_out(:,:,k);
        
        E_r_out(~mask_out)  = 0;
        E_th_out(~mask_out) = 0;
        
        Ex = Ex + (E_r_out .* (sinth .* cosphi) + E_th_out .* (cos_theta .* cosphi));
        Ez = Ez + (E_r_out .* cos_theta - E_th_out .* sinth);

    end
end

% Mask points inside singularity if r==0
Ex(r==0)=0; Ez(r==0)=0;

% Optionally mask inside/outside regions to show only valid region
mask = true(size(r));
% interior terms plotted for r<=R if any interior present, exterior for r>=R if any exterior present
if ~isempty(l_in) && isempty(l_out)
    mask = r<=1.01*R;
elseif isempty(l_in) && ~isempty(l_out)
    mask = r>=0.99*R;
elseif ~isempty(l_in) && ~isempty(l_out)
    % show both; keep whole plane
    mask = true(size(r));
end


% downsample for plotting
skip = max(1, round(ngrid/25));
Xs = X(1:skip:end,1:skip:end);
Zs = Z(1:skip:end,1:skip:end);
Exs = Ex(1:skip:end,1:skip:end).*mask(1:skip:end,1:skip:end);
Ezs = Ez(1:skip:end,1:skip:end).*mask(1:skip:end,1:skip:end);

% magnitudes
mag = sqrt(Exs.^2 + Ezs.^2);
maxm = max(mag(:));
minm = min(mag(mag>0));

if ~isempty(arrowScale) && arrowScale == 0
    % Equal-length arrows mode
    arrowLength = 3 * r_plane / ngrid;           % default fixed arrow length (tunable)
    % unit directions
    ux = zeros(size(mag)); uz = zeros(size(mag));
    nz = mag > eps;
    ux(nz) = Exs(nz) ./ mag(nz);
    uz(nz) = Ezs(nz) ./ mag(nz);
    Uplot = arrowLength * ux;
    Vplot = arrowLength * uz;
    autoscale_flag = 0;                      % quiver uses absolute values
else
    % Normal variable-length mode: compute suggested scale if arrowScale empty
    if isempty(arrowScale) || ~isnumeric(arrowScale) || ~isfinite(arrowScale)
        if maxm > 0
            base = 0.18 * r_plane / maxm;
            if ~isempty(minm) && minm>0 && maxm/minm > 1e2
                dyn = sqrt(minm*maxm);
                base = 0.18 * r_plane / dyn;
            end
        else
            base = 1;
        end
        scale = base;
    else
        scale = arrowScale;
    end
    Uplot = Exs * scale;
    Vplot = Ezs * scale;
    autoscale_flag = 0;                      % set to 0 to use values directly
end

figure('Color','w')
quiver(Xs, Zs, Uplot, Vplot, autoscale_flag, 'k');
axis equal; hold on
th = linspace(0,2*pi,200);
plot(R*cos(th), R*sin(th), 'b--', 'LineWidth',1);
xlabel('x'); ylabel('z');
if exist('arrowLength','var')
    title('Electric field (equal-length arrows)');
else
    title('Electric field (variable-length arrows)');
end
xlim([-r_plane r_plane]); ylim([-r_plane r_plane])
hold off

end