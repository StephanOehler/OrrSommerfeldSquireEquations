classdef chebyshev
    % Generates all you need if you want to work with chebyshev polynomials
    % Supports up to fourth order differentiation.
    % The entire code is based on Trefethen's Spectral Methods in Matlab
    %
    % Inputs:
    % N: chebyshev dimension
    %        Default = 100;
    % h:     linear scaling factor
    %        Default = 1
    % integration: integration method:
    %        chose {'clenshawcurtis'.'trapezoidal'}
    %        Default = 'clenshawcurtis';
    %
    % Outputs:
    % x  : grid_points (N chebyshev grid points)
    % b  : barycentric interpolation weights
    % D1 : first order differentiation matrix
    % D2 : second order differentiation matrix
    % D3 : third order differentiation matrix
    % D4 : fourth order differentiation matrix
    % w : the integration weights
    %
    % Call methods:
    % Example1 = chebyshev()
    % Example2 = chebyshev(nx)
    % Example3 = chebyshev(nx,h)
    % Example5 = chebyshev(N,h,integration)
    %%  Private Properties
    properties
        %% Discretisation
        N = 100;
        h = 1;
        integration = 'clenshawcurtis';
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Discretised System term
        x % unscaled chebyshev grid
        b % barycentric interpolation weights
        w % clenshaw curtis integration weights
        D1 % Derivatives
        D2
        D3
        D4
    end
    %% Unscaled Chebyshev terms
    properties (Dependent, Access = protected)
        xh % pure, unstacled chebyshev points
        DC % 1st order differentiation, unscaled
        DC2 % 2nd order conditionless, unscaled
        DC3 % 3rd order conditionless, unscaled
        DC4 % 4th order conditionless, unscaled
    end
    %% Methods
    methods
        %% Setup function - default overwrite
        function thischebyshev = chebyshev(N,h,integration)
            if nargin == 0 || nargin > 3
            else
                thischebyshev.N = N;
                if nargin > 1
                    thischebyshev.h = h;
                end
                if nargin > 2
                    thischebyshev.integration = integration;
                end
            end
        end
        %% Identify
        function identify(cheb)
            disp(['I generate everything you need for chebyshev integration and up to fourth order differentiation with N = ',num2str(cheb.N),' chebyshev points'])
        end
        %% Get DC - unscaled differentiation matrix - c stands for conditionless
        function DC = get.DC(cheb)
            DC = gl_chebD(cheb.N);
        end
        %% Get second order unstacled differentiation matrix
        function DC2 = get.DC2(cheb)
            DC2 =cheb.DC^2;
        end
        %% Get third order unstacled differentiation matrix
        function DC3 = get.DC3(cheb)
            DC3 =cheb.DC2 * cheb.DC;
        end
        %% Get fourth order unstacled differentiation matrix
        function DC4 = get.DC4(cheb)
            DC4 =cheb.DC3 * cheb.DC;
        end
        %% Get xh - unscaled chebyshev points with boundary conditions applied
        function xh = get.xh(cheb)
            xh = gl_chebx(cheb.N);
            xh( abs (xh) < 1e-15) = 0; % ensures that the centre point is
            %zero, might need this later for barycentric interpolation
        end
        %% Get xgrid - scaled chebyshev points with boundary conditions applied
        function x = get.x(cheb)
            x = cheb.xh * cheb.h;
        end
        %% Get D1 - scaled first order differentiation matrix
        function D1 = get.D1(cheb)
            D1 = cheb.DC / cheb.h;
        end
        %% Get D2 - scaled second order differentiation matrix
        function D2 = get.D2(cheb)
            D2 = cheb.DC2 / cheb.h^2;
        end
        %% Get D3 - scaled second order differentiation matrix
        function D3 = get.D3(cheb)          
            D3 = cheb.DC3 / cheb.h ^ 3;
        end
        %% Get D4 - scaled second order differentiation matrix
        function D4 = get.D4(cheb)
                D4 = cheb.DC4 / cheb.h ^ 4;
        end
        %% Get b, the barycentric interpolation weights
        function b = get.b(cheb)
            b = [1/2; ones(cheb.N-1,1); 1/2] .* (-1).^((0:cheb.N)');
        end
        %% Get w (diag(Q)) - scaled integration matrix with boundary conditions applied
        function w = get.w(cheb)
            if strcmp(cheb.integration,'clenshawcurtis')
                w = gl_chebw(cheb.N);
            elseif strcmp(cheb.integration,'trapizoidal')
                w = ( ([ diff(cheb.x); 0 ] + [ 0; diff(cheb.x) ]) / 2 ); % Trapezoidal rule
            end
            w = cheb.h * w';
        end
        %% Restrict integration
        function cheb = set.integration(cheb,integration)
            if  any(strcmp(integration,{'clenshawcurtis','trapizoidal'}))
                cheb.integration = integration;
            else
                error('integration must be either: ''clenshawcurtis'',''trapizoidal''')
            end
        end
        %% Restrict N
        function cheb = set.N(cheb,N)
            if N > 0
                cheb.N = N;
            else
                error('N must be > 0')
            end
        end
    end
end
%% Generate chebyshev polynomials
% See Trefethen, Spectral Methods in Matlab
function x = gl_chebx(N)
x = cos( pi * (0:N) / N)';
end
%% Generate chebyshev differentiation matrix
% See Trefethen, Spectral Methods in Matlab
function D = gl_chebD(N)
x = gl_chebx(N);
c = [2; ones(N - 1 , 1); 2] .* (-1).^(0:N)';
X = repmat(x , 1 , N + 1 );
dX = X - X';
D = ( c * (1./ c)') ./ (dX + (eye (N + 1)));
D = D - diag( sum(D'));
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
end