classdef StateSpace < channelOSS.OrrSommerfeldSquire
    % This class turns the Orr Sommerfeld Squire equation into a state
    % space model.
    properties (Dependent, Access = public)
        %% Discretised System terms
        A % System matrix
        B % Standard Input matrix
        Bw % Weighted Input matrix (random forcing everywhere)
        Bz % Input at specific wall height.
        Cx % Recovers the boundary conditions
        C % Standard output matrix
        Cw % Weighted Output matrix
        Cz % Output at specific wall height
        Ctaup % shear stress
        %% only needed if Bbs is used (use in conjunction with Cx, C or Cw).
        Bbs % Bbs = B matix for blowing and suction  
        Dx 
        D; 
        Dw;
        Dz;
        Dtaup;
    end
    properties
        %% Flow Properties
        half = 'both'; % top bot both
        N_out = 200; % Output resolution
        Cb % Uses baryentric interpolation to map to output resolution. 
        Bb % Uses baryentric interpolation to map from output resolution. 
        zBz = 0; % actuator location
        zCz = 0; %sensor location
    end
    properties(Dependent, Access = public)
        z_out % Output planes (chebyshev points in z - wall normal)
        w_out % Corresponding output clenshaw curtis weights
    end
    
    properties (Dependent, Access = private)
        Losq; % OSS operator
        E; % E matrix (LHS)
        CC; % temporary C matrix
        CCtaup; % temporary Ctaup matrix
        CCz; % temporary Cz matrix
    end
    
    methods
        %% Identify
        function identify(OSS)
            disp(['I am an Orr-Sommerfeld-Squire State space model with N = ',num2str(OSS.N),' chebyshev states'])
        end
        %% Setup function
        function OSS = StateSpace(kx,ky,Re,N,type,half,N_out,h,zBz,zCz,Nnu)
            if nargin == 0 || nargin > 9
            else
                OSS.kx = kx;
                if nargin > 1
                    OSS.ky = ky;
                end
                if nargin > 2
                    OSS.Re = Re;
                end
                if nargin > 3
                    OSS.N = N;
                end
                if nargin > 4
                    OSS.type = type;
                end
                if nargin > 5
                    OSS.half = half;
                end
                if nargin > 6
                    OSS.N_out = N_out;
                end
                if nargin < 7
                    OSS.N_out = N;
                end
                if nargin > 7
                    OSS.h = h;
                end
                if nargin > 8
                    OSS.zBz = zBz;
                end
                if nargin > 9
                    OSS.zCz = zCz;
                end
                if nargin > 10
                    OSS.Nnu = Nnu;
                end
                if nargin < 11
                    OSS.Nnu = OSS.N;
                end
            end
        end
        %% Get ss operator Losq
        function Losq = get.Losq(OSS)
            if OSS.k == 0
                if any(strcmp(OSS.type,{'laminar-channel','laminar-cuette','turbulent-noeddy'}))
                    Losq = OSS.geom.D2 / OSS.Re;
                elseif any(strcmp(OSS.type,{'turbulent'}))
                    Losq = diag(OSS.profile.nu_T) * OSS.geom.D2 + diag(OSS.profile.dnu_T) * OSS.geom.D1;
                end             
                Losq = [Losq OSS.Z; OSS.Z Losq];
            else
                Losq =  [OSS.Los , OSS.Z ; ...
                    OSS.Lac, OSS.Lsq];
            end
            
        end
        %% Get ss operator A
        function A = get.A(OSS) 
            LLosq = OSS.Losq(OSS.geom.ik,OSS.geom.ik) + OSS.Losq(OSS.geom.ik,OSS.geom.ir) * OSS.geom.cG; %implement boundary conditions
            A = OSS.E \ LLosq;
        end
        %% Left hand side
        function E = get.E(OSS)
            if OSS.k == 0
                E = [OSS.I, OSS.Z; OSS.Z OSS.I ];
            else
                E = [OSS.Delta, OSS.Z; OSS.Z OSS.I ];
            end
            E = E(OSS.geom.ik,OSS.geom.ik) + E(OSS.geom.ik,OSS.geom.ir) * OSS.geom.cG;
        end
        %% Get ss operator F
        function B = get.B(OSS)
            if OSS.k == 0
                B = [OSS.I OSS.Z OSS.Z; OSS.Z OSS.I OSS.Z];
            else
                B =  [-1i * OSS.kx * OSS.geom.D1, -1i * OSS.ky * OSS.geom.D1, -OSS.k^2 * OSS.I; ...
                    1i * OSS.ky * OSS.I,       -1i * OSS.kx * OSS.I,        OSS.Z];
            end      
            B = OSS.E \ B(OSS.geom.ik,:);
        end
        %% Get Cx the recovery matrix
        function Cx = get.Cx(OSS)
            I_kk = eye(length(OSS.geom.ik));
            Cx = [I_kk; OSS.geom.cG];
            Cx = Cx(OSS.geom.ii,:);
        end
         %% Get ss operator D, which is only needed if BBs is used (use in conjunction with C).
        function Dx = get.Dx(OSS)
            Dx =  OSS.geom.cC \ eye(length(OSS.geom.ir));
        end
        %% Get ss operator CC 
        function CC = get.CC(OSS)
            if OSS.k == 0
                CC = [OSS.I OSS.Z; OSS.Z OSS.I; OSS.Z OSS.Z];
            else
                D1Dirichelet = OSS.geom.D1;
                D1Dirichelet([1,OSS.geom.N + 1],:) = 0;
                CC =  [ 1i * OSS.kx * D1Dirichelet, -1i * OSS.ky * OSS.I; ...
                    1i * OSS.ky * D1Dirichelet,  1i * OSS.kx * OSS.I; ...
                    OSS.k^2 * OSS.I, OSS.Z] / (OSS.k^2);
            end
        end
        %% Get ss operator C
        function C = get.C(OSS)
            C = OSS.Cb * OSS.CC * OSS.Cx;
        end
        %% Get ss operator D
        function D = get.D(OSS)
            D = OSS.Cb * OSS.CC * OSS.Dx;
        end
        %% Get energy weighted forcing shapes Bw
        function Bw = get.Bw(OSS)
            iw = 1 ./ sqrt(OSS.geom.w);
            Bw = OSS.B * diag([iw; iw; iw]);
        end
        %% Get energy weighted forcing velocities Cw
        function Cw = get.Cw(OSS)
            w = sqrt(OSS.w_out);
            Cw = diag([w; w; w]) * OSS.C;
        end
        %% Get energy weighted forcing velocities Dw
        function Dw = get.Dw(OSS)
            w = sqrt(OSS.w_out);
            Dw = diag([w; w; w]) * OSS.D;
        end
        %% blowing and suction
        function Bbs = get.Bbs(OSS) 
           
            Bbs = OSS.E \ OSS.Losq(OSS.geom.ik,OSS.geom.ir) * OSS.geom.cH;

            if OSS.k == 0
                E_temp = [OSS.I, OSS.Z; OSS.Z OSS.I ];
            else
                E_temp = [OSS.Delta, OSS.Z; OSS.Z OSS.I ];
            end
            Bbs2 = - OSS.E \ E_temp(OSS.geom.ik,OSS.geom.ir) * OSS.geom.cH;
            Bbs = [Bbs,Bbs2];
        end
        %% wall measurements
        function CCtaup = get.CCtaup(OSS) 
            % dx, dy, dz^3h
            C0{1} = OSS.geom.D1([1,end],:) / OSS.Re;
            C0{2} = C0{1};
            D1tilde = OSS.geom.D1;
            D1tilde([1,end],:) =  0;
            D3tilde = OSS.geom.D2 * D1tilde;
            if OSS.k == 0
                C0{3} = zeros(size(D3tilde([1,end],:)));
            else
                C0{3} = D3tilde([1,end],:) / OSS.Re / OSS.k^2;
            end
            % build overall C matrix in terms of C just found
            q = 2;
            p = OSS.geom.N + 1;
            Z = zeros(q,p);
            % indices corresponding to...
            ind{1} = 1:p;    % streamwise velocity components
            ind{2} = p + ind{1}; % wall-normal
            ind{3} = p + ind{2}; % spanwise
            % for each velocity component in vel, construct relevant rows of CC
            CCtaup = []; % initialize CC to an empty matrix
            for i = 1:3
                temp = [Z Z Z];
                temp(:,ind{i}) = C0{i};
                CCtaup = [CCtaup;temp]; % append to Cuy the rows found in temp
            end
        end
        %% wall measurements
        function Ctaup = get.Ctaup(OSS)
            Ctaup = OSS.CCtaup * OSS.CC * OSS.Cx;
        end
        %% wall measurements
        function Dtaup = get.Dtaup(OSS)
            Dtaup = OSS.CCtaup * OSS.CC * OSS.Dx;
        end        
        %% Cz point sensor (can handle multiple points)
        function CCz = get.CCz(OSS)
            q = length(OSS.zCz);
            p = OSS.N + 1; 
            Cint = zeros(q, p);
            for j = 1:q % for each interpolant grid point (each row of B)
                % Barycentric formula as matrix
                diff = OSS.zCz(j) - OSS.geom.x;
                temp = OSS.geom.b./ diff;
                Cint(j,:) = temp / sum(temp);
                % deal with case where z(j)=(one of Chebyshev points)
                if min(abs(diff)) == 0
                    Cint(j,:) = 0;
                    Cint(j, diff == 0) = 1;
                end
            end
            Z = zeros(q,p);
            % indices corresponding to...
            ind{1} = 1:p;      % streamwise velocity components
            ind{2} = p + ind{1}; % wall-normal
            ind{3} = p + ind{2}; % spanwise
            % for each velocity component in vel, construct relevant rows
            % of BB
            CCz = []; % initialize BB to an empty matrix
            for i = 1:3
                temp = [Z Z Z];
                temp(:,ind{i}) = Cint;
                CCz = [CCz; temp]; % append to Bu the rows found in temp
            end
        end
         %% Cz point sensor (can handle multiple points)
        function Cz = get.Cz(OSS)
           Cz = OSS.CCz * OSS.CC * OSS.Cx; % should I scale the three velocity components with diag(OSS.geom.b)?
        end 
        %% Dz point sensor (can handle multiple points)
        function Dz = get.Dz(OSS)
           Dz = OSS.CCz * OSS.CC * OSS.Dx; % should I scale the three velocity components with diag(OSS.geom.b)?
        end          
        %% Bz point actuator (can handle multiple points)
        function Bz = get.Bz(OSS)
            q = length(OSS.zBz);
            p = OSS.N + 1; 
            Bint = zeros(q, p);
            for j = 1:q % for each interpolant grid point (each row of B)
                % Barycentric formula as matrix
                diff = OSS.zBz(j) - OSS.geom.x;
                temp = OSS.geom.b./ diff;
                Bint(j,:) = temp / sum(temp);
                % deal with case where z(j)=(one of Chebyshev points)
                if min(abs(diff)) == 0
                    Bint(j,:) = 0;
                    Bint(j, diff == 0) = 1;
                end
            end
            Z = zeros(q,p);
            % indices corresponding to...
            ind{1} = 1:p;      % streamwise velocity components
            ind{2} = p + ind{1}; % wall-normal
            ind{3} = p + ind{2}; % spanwise
            % for each velocity component in vel, construct relevant rows
            % of BB
            Bz = []; % initialize Bz to an empty matrix
            for i = 1:3
                temp = [Z Z Z];
                temp(:,ind{i}) = Bint;
                Bz = [Bz; temp]; % append to Bu the rows found in temp
            end
            Bz = OSS.B * Bz.';
        end
        %% set w_out
        function w_out = get.w_out(OSS)
            w_out = gl_chebw(OSS.N_out);
            switch OSS.half
                case 'bot'; w_out = w_out / 2;
                case 'top'; w_out = w_out / 2;
                case 'both'
            end
        end
        %% set z_out
        function z_out = get.z_out(OSS)
            z_out = gl_chebx(OSS.N_out);
            switch OSS.half
                case 'bot'; z_out = (z_out - 1) / 2;
                case 'top'; z_out = (z_out + 1) / 2;
                case 'both'
            end
        end
        %% set Cb
        function Cb = get.Cb(OSS)    
            b = [1/2; ones(OSS.N-1,1); 1/2] .* (-1).^((0:OSS.N)');
            % get the interpolating matrix C
            Cb_temp = zeros(OSS.N_out+1,OSS.N+1);
            for j = 1:(OSS.N_out+1) % for each interpolant grid point (each row of C)
                % Barycentric formula as matrix
                diff = OSS.z_out(j)-OSS.geom.x;
                temp = b./diff;
                Cb_temp(j,:) = temp/sum(temp);
                % deal with case where z(j)=(one of Chebyshev points)
                if min(abs(diff))==0
                    Cb_temp(j,:) = 0;
                    Cb_temp(j,diff==0) = 1;
                end
            end
            Z = zeros(OSS.N_out+1,OSS.N+1);
            Cb = [Cb_temp,Z,Z;Z,Cb_temp,Z;Z,Z,Cb_temp];
        end%% set Cb
        %% set Bb
        function Bb = get.Bb(OSS)    
            b = [1/2; ones(OSS.N_out-1,1); 1/2] .* (-1).^((0:OSS.N_out)');
            % get the interpolating matrix C
            Bb_temp = zeros(OSS.N+1,OSS.N_out+1);
            for j = 1:(OSS.N+1) % for each interpolant grid point (each row of C)                
                % Barycentric formula as matrix
                diff = OSS.geom.x(j)-OSS.z_out;
                temp = b./diff;
                Bb_temp(j,:) = temp/sum(temp);
                % deal with case where z(j)=(one of Chebyshev points)
                if min(abs(diff))==0
                    Bb_temp(j,:) = 0;
                    Bb_temp(j,diff==0) = 1;
                end
            end
            switch OSS.half
                case 'bot'; Bb_temp(1:(floor(OSS.N/2)+1),:) = 0;
                case 'top'; Bb_temp((ceil(OSS.N/2) + 1):end,:) = 0;
                case 'both'
            end
            Z = zeros(OSS.N+1,OSS.N_out+1);
            Bb = [Bb_temp,Z,Z;Z,Bb_temp,Z;Z,Z,Bb_temp];
        end
        %% Restrict half
        function cheb = set.half(cheb,half)
            if any(strcmp(half,{'both','top','bot'}))
                cheb.half = half;
            else
                error('half must be either ''both'',''top'', or ''bot''')
            end
        end   
        %% Restrict half N
        function cheb = set.N_out(cheb,N_out)
            if N_out > 0
                cheb.N_out = N_out;
            else
                error('N must be > 0')
            end
        end
        %% Get energy norm
        function Output = GetH2(OSS,norm)
            if all(ismember(norm,[1,2,3,4])) ~= 1
                error('Give norm type as two digit vector: first entry: x = 1, y = 2, z = 3, all = 4, second entry: u = 1, v = 2, w = 3, all = 4. Example: [2,3] gives H2yw ')
            end
            %% Get correct B
            if norm(1) == 4
                Bin = OSS.B;
            else
                Bin = OSS.B(:,( (norm(1) - 1) * OSS.N + 1):( norm(1) * OSS.N));
            end
            if norm(2) == 4
                Cout = OSS.C ;
            else
                Cout = OSS.C(( (norm(2) - 1) * OSS.N + 1):(norm(2) * OSS.N),:);
            end
            %% Get norm
            X = lyap(OSS.A,Bin*Bin');
            Output = trace(X * (Cout' * Cout));
        end
        %% Get eigenvalues
        function [eigOS, eigSq] = Geteig(OSS)
            % Get Orr-Sommerfeld eigenvalues
            eigOS = eig(OSS.Los,OSS.Delta);
            [~,ind] = sort(real(eigOS),'descend');
            eigOS = eigOS(ind);
            % Get Squire eigenvalues
            eigSq = eig(OSS.Lsq);
            [~,ind] = sort(real(eigSq),'descend');
            eigSq = eigSq(ind);
        end
        %% G - give initial time, final time, and number of steps
        function [G, t] = GetG(OSS,t0,tf,nt)
            %do an eigenvalue decomposition
            [X, D] = eig(OSS.A); d = diag(D);
            [~, ind] = sort(real(d), 'descend'); d = d(ind); X = X(:,ind);
            
            % form F matrix in F*F=M. This comes from the energy-weighted C matrix
            CW = OSS.Cw; M=CW'*CW; [F,~]=sqrth(M);
            
            % form transient growth curve
            t = linspace(t0,tf,nt); dt = t(2)-t(1);
            H = eye(2 * OSS.N-4);
            % note two ways to get evolution matrix--both give same result
            expmat = X * diag(exp(d * dt)) / X; %expmat=expm(A*dt); this is slower
            G = zeros(1,nt);
            for tind = 1:nt
                G(tind)=norm(F * H / F)^2;
                H = H * expmat;
            end
        end
    end
end

%% compute decomposition A=F'*F and inv(F)
% using singular value decomposition
% (source forgotten but not mine)
function varargout=sqrth(A)
% 
%
% If you want to decompose an hermitian matrix
% as A=F*F' (covariance stuff)
% you can do F=sqrth(A)' ; and check with norm(A-F*F')

[Uexp,Sexp,~]=svd(A);
s=sqrt( diag(Sexp) );
F=diag(s)*Uexp';
Finv=Uexp*diag(ones(size(s))./s);

switch nargout
    case 1
        varargout{1}=F;
    case 2
        varargout{1}=F;
        varargout{2}=Finv;
    otherwise
end
end

%% Generate chebyshev polynomials
% See Trefethen, Spectral Methods in Matlab
function x = gl_chebx(N)
x = cos( pi * (0:N) / N).';
x( abs (x) < 1e-15) = 0; % ensures that the centre point is
end

%% Generate chebyshev integration vector
% See Trefethen, Spectral Methods in Matlab
function w = gl_chebw(N)
theta = pi*(0:N)'/N;
w = zeros(1,N+1);
ii = 2:N;
v = ones(N-1,1);
if mod(N,2)==0
    w(1) = 1/(N^2-1);
    w(N+1) = w(1);
    for k=1:N/2-1
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
end
w(ii) = 2*v/N;
w = w.';
end
