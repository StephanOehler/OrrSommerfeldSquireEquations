classdef OrrSommerfeldSquire
    % This class generates the terms necessary for the Orr Sommerfeld
    % Squire Equations  
    %%  Private Properties
    properties
        %% Flow Properties - default
        N = 200; % resolution
        type = 'turbulent'
        %'turbulent' 'laminar-channel' 'laminar-cuette'\ 'turbulent-noeddy'
        h = 1; % wall half height
        kx = 0.5; % streamwise wavenumber
        ky = 1; % spanwise wavenumber
        Re = 2000; % Re_\tau friction Reynolds number
        Nnu; % resolution of the viscosity profile
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Outputs
        geom; % channel geometry class
        profile; % class that generates the mean velocity profile
        Los; % Orr sommerfeld operator
        Lsq; % Squire operator
        Lac; % cross terms
        
        Delta; % delta operator
        I; % Identity
        Z; % zeros
        k; % sqrt(kx^2 + ky^2)
    end
    %% Protected terms
    properties (Dependent, Access = protected)
        Delta2; % delta squared
    end
    %% Methods
    methods
        %% Setup function - default overwrite
        function thisOSS = OrrSommerfeldSquire(kx,ky,Re,N,type,h,Nnu)
            if nargin == 0 || nargin > 7
            else
                thisOSS.kx = kx;
                if nargin > 1
                    thisOSS.ky = ky;
                end
                if nargin > 2
                    thisOSS.Re = Re;
                end
                if nargin > 3
                    thisOSS.N = N;
                end
                if nargin > 4
                    thisOSS.type = type;
                end
                if nargin > 5
                    thisOSS.h = h;
                end
                if nargin > 6
                    thisOSS.Nnu = Nnu;
                end
            end
            if nargin ~= 7
                thisOSS.Nnu = thisOSS.N;
            end
        end
        %% Identify
        function identify(OrrSommerfeldSquire)
            disp('I am the Orr Sommerfeld Squire equation!')
        end
        %% Get geom
        function geom = get.geom(OSS)
            geom =  channelOSS.Geom(OSS.N,OSS.h);
        end
        %% Get profile
        function profile = get.profile(OSS)
            if any(strcmp(OSS.type,{'laminar-channel','laminar-cuette'}))
                if strcmp(OSS.type,'laminar-channel')
                    profile =  channelOSS.LamProf(OSS.geom.x,'channel',OSS.h);
                    
                elseif strcmp(OSS.type,'laminar-cuette')
                    profile =  channelOSS.LamProf(OSS.geom.x,'cuette',OSS.h);
                end
            elseif any(strcmp(OSS.type,{'turbulent','turbulent-noeddy'}))
                profile =  channelOSS.ViscProf(OSS.geom.x,OSS.Nnu,OSS.h,OSS.Re);
            end
        end
        %% Get I (identity)
        function I = get.I(OSS)
            I =  (eye(OSS.N + 1));
        end
        %% Get Z (zeros)
        function Z = get.Z(OSS)
            Z =  (zeros(OSS.N + 1));
        end
        %% Get k2 operator
        function k = get.k(OSS)
            k =  sqrt(OSS.kx^2 + OSS.ky^2);
        end
        %% Get Delta operator
        function Delta = get.Delta(OSS)
            Delta = OSS.geom.D2 - OSS.I * OSS.k^2;
        end
        %% Get Delta2 operator
        function Delta2 = get.Delta2(OSS)
            Delta2 = OSS.geom.D4 + OSS.I * OSS.k^4 - 2 * OSS.k^2 * OSS.geom.D2;
        end
        %% Get Orr-Sommerfeld Operator
        function Los = get.Los(OSS)
            if any(strcmp(OSS.type,{'laminar-channel','laminar-cuette','turbulent-noeddy'}))
                Los = - 1i * OSS.kx * diag(OSS.profile.U) * OSS.Delta + ...
                    1i * OSS.kx * diag(OSS.profile.ddU) + OSS.Delta2 / OSS.Re;
            elseif strcmp(OSS.type,'turbulent')
                Los = - 1i * OSS.kx * diag(OSS.profile.U) * OSS.Delta + ...
                    1i * OSS.kx * diag(OSS.profile.ddU) + ...
                    diag(OSS.profile.nu_T) * OSS.Delta2 + ...
                    2 * diag(OSS.profile.dnu_T) * ( OSS.geom.D3 - OSS.k^2 * OSS.geom.D1 ) + ...
                    diag(OSS.profile.ddnu_T) * (OSS.geom.D2 + OSS.k^2 * OSS.I);
            end
        end
        %% Squire
        function Lsq = get.Lsq(OSS)
            if any(strcmp(OSS.type,{'laminar-channel','laminar-cuette','turbulent-noeddy'}))
                Lsq = - 1i * OSS.kx *  diag(OSS.profile.U) + ...
                    OSS.Delta / OSS.Re;
            elseif strcmp(OSS.type,'turbulent')
                Lsq = - 1i * OSS.kx *  diag(OSS.profile.U) + ...
                    diag(OSS.profile.nu_T) * OSS.Delta + ...
                    diag(OSS.profile.dnu_T) * OSS.geom.D1;
            end
        end
        %% Cross term
        function Lac = get.Lac(OSS)
            Lac = -1i * OSS.ky * diag(OSS.profile.dU); % cross term
        end
        %% Restrict modes
        function OSS = set.type(OSS,type)
            if any(strcmp(type,{'laminar-channel','laminar-cuette','turbulent','turbulent-noeddy'}))
                OSS.type = type;
            else
                error('type must be either: ''laminar-channel'', ''laminar-cuette'', ''turbulent'', ''turbulent-noeddy''')
            end
        end
    end
end
