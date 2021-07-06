classdef ViscProf
    % generates an eddy viscosity-based velocity profile
    % ViscProf and LamProf use y as wall-normal component, every other
    % code uses z as wall-normal component
    
    properties
        %% Channel parameters
        y = []; % grid to interpolate onto
        Nnu = 100 % dimensions of profile
        h = 1; % wall-height
        Re = []; % Reynolds number
        Aplus = 25.4; % viscosity constant
        kappa = 0.426; % viscosity constant
        utau = 1; % velocity constant
        z01; % interpolation variable
        z; % interpolation variable
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Public, dependent
        nu_T; % Eddy viscosity profile and derivatives
        dnu_T;
        ddnu_T;
        U; % mean velocity profile and derivatives
        dU
        ddU;
        nu;
    end
    properties (Dependent, Access = protected)
        eta; % grid points for velocity (half channel)
        geom; % geometry class
        nu_T_half;
        U_half;
        dU_half;
        yind; % indicies
    end
    %%
    methods
        %% Setup function - default overwrite
        function chanViscProf = ViscProf(y,Nnu,h,Re)
            if nargin == 0 || nargin > 4
            else
                chanViscProf.y = y;
                if nargin > 1
                    chanViscProf.Nnu = Nnu;
                end
                if nargin > 2
                    chanViscProf.h = h;
                end
                if nargin > 3
                    chanViscProf.Re = Re;
                end
            end
        end
        %% Get geom: geometry of the channel
        function geom = get.geom(ViscProf)
            geom = channelOSS.chebyshev(ViscProf.Nnu,ViscProf.h);
        end
        %% Get indicies for full channel height
        function yind = get.yind(ViscProf)
            yind.n0 = find(ViscProf.y <= 0); % case 1: y \in [0,1]
            yind.n =sort(find(ViscProf.y < 0),'descend'); % indices for y \in [-1,0)
            yind.p =sort(find(ViscProf.y > 0),'ascend');  % indices for y \in (0,1]
        end
        %% Get the grid from 0 to 1
        function z01 = get.z01(ViscProf)
            z01 = (ViscProf.geom.x +  1) / 2;
        end
        %% Get eta, the grid that we generate our profiles on
        function eta = get.eta(ViscProf)
            eta = 1 - ViscProf.z01;
        end
        %% z the ture grid from -1 to 0
        function z = get.z(ViscProf)
            z = ViscProf.z01 - 1;
        end
        %% Get half of the nu_T profile
        function nu_T_half = get.nu_T_half(ViscProf)
            nu_Tplus = ViscProf.kappa .* ViscProf.Re ./ 3 .* (1 - ViscProf.eta.^2 ) .* ( 1 + 2 .* ViscProf.eta.^2) .* ( 1 - exp( ( abs(ViscProf.eta) - 1 ) .* ViscProf.Re / ViscProf.Aplus ));
            nu_tplus = 0.5 .* ( 1 + nu_Tplus.^2 ).^0.5 + 0.5;
            nu_T_half = nu_tplus * ViscProf.nu;
        end
        %% Get nu - for non-dimensionalisation
        function nu = get.nu(ViscProf)
            nu = ViscProf.utau * ViscProf.h / ViscProf.Re;
        end
        %% Get half of dU
        function dU_half = get.dU_half(ViscProf)
            dUplus = (ViscProf.Re * (ViscProf.eta ./ ViscProf.nu_T_half) * ViscProf.nu );
            dU_half = ViscProf.utau * dUplus; % uplus = u / utau
        end
        %% Get half of U
        function U_half = get.U_half(ViscProf)
            U_half = zeros(ViscProf.Nnu + 1, 1);
            for i = 1:(ViscProf.Nnu + 1)
                z0 = ViscProf.z01(i); % we want to integrate between [1,eta0]
                zlim = ViscProf.z01 * z0; % grid points on which to integrate: we need to flip to get  [1,eta0]
                f = bary_interp(ViscProf.z01, zlim, ViscProf.geom.b, ViscProf.dU_half); % interpolate dU
                U_half(i) = ( (ViscProf.geom.w)' * f * (z0) / 2 ); % integrate using Clenshaw-Curtis
            end
        end
        %% Get U
        function U = get.U(ViscProf)
            U = FullInterpolate(ViscProf,ViscProf.U_half,'sym'); % multiply by 2 because z \in [0,1] rather than [-1,1]
        end
        %% Get dU
        function dU = get.dU(ViscProf)
            dU = FullInterpolate(ViscProf,ViscProf.dU_half,'asym'); % multiply by 2 because z \in [0,1] rather than [-1,1]
        end
        %% Get ddU
        function ddU = get.ddU(ViscProf)
            ddU = FullInterpolate(ViscProf,ViscProf.geom.D1 * ViscProf.dU_half * 2,'sym'); % multiply by 2 because z \in [0,1] rather than [-1,1]
        end
        
        %% Get nu_T
        function nu_T = get.nu_T(ViscProf)
            nu_T = FullInterpolate(ViscProf,ViscProf.nu_T_half,'sym'); % multiply by 2 because z \in [0,1] rather than [-1,1]
        end
        %% Get dnu_T
        function dnu_T = get.dnu_T(ViscProf)
            dnu_T = FullInterpolate(ViscProf,ViscProf.geom.D1 * ViscProf.nu_T_half * 2,'asym'); % multiply by 2 because z \in [0,1] rather than [-1,1]
        end
        %% Get ddnu_T
        function ddnu_T = get.ddnu_T(ViscProf)
            ddnu_T = FullInterpolate(ViscProf,ViscProf.geom.D2 * ViscProf.nu_T_half * 4,'sym'); % multiply by 4 because z \in [0,1] rather than [-1,1]
        end
    end
end

function ff = bary_interp(x , xx , c , f)

N = length(x);

numer = zeros(size(xx));
denom = zeros(size(xx));

exactx  = zeros(size(x));
exactxx = zeros(size(xx));

for j = 1:N % for each Chebyshev grid point
    xdiff = xx-x(j);
    temp = c(j)./xdiff;
    numer = numer + temp*f(j);
    denom = denom + temp;
    if min(abs(xdiff))==0, exactx(j) = 1; end
    exactxx(xdiff==0) = 1;
end
ff = numer./denom;

indx = find(exactx); indxx = find(exactxx);
ff(indxx) = f(indx);
end

function Out = FullInterpolate(ViscProf,Inp,arg)
Out(ViscProf.yind.n0) = bary_interp(ViscProf.z,ViscProf.y(ViscProf.yind.n0),ViscProf.geom.b,  Inp);
if strcmp(arg,'sym')
    Out(ViscProf.yind.p) =+  Out(ViscProf.yind.n); % symmetric
elseif strcmp(arg,'asym')
    Out(ViscProf.yind.p) =-  Out(ViscProf.yind.n); % asymmetric
end
end