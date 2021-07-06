classdef Geom < channelOSS.chebyshev
    % Generates the geometry for the channel flow by applying suitable
    % boundary conditions.
    % Created using 'Implementation of boundary condition' 
    % by Jerome Hoepffner (June 2007)    
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Discretised System terms
        ik % DOF that are kept
        ir % DOF that are removed
        ii % sorted DOFS
        cC % C matrix
        cG % Give back matrix
        cH % non-homogenou boundary conditions
    end
    
    properties (Dependent, Access = private)
        %% Discretised System terms    
    end  
    %% Methods
    methods
        %% Setup function - default overwrite
        function thisgeom = Geom(N,h,integration)
            if nargin == 0 || nargin > 3
            else
                thisgeom.N = N;
                if nargin > 1
                    thisgeom.h = h;
                end
                if nargin > 2
                    thisgeom.integration = integration;
                end
            end
            %thisgeom = thisgeom.updateN;
        end
        %% Get ik: kept degrees of freedom
        function ik = get.ik(geom)
            ik = [(3:geom.N-1)';(2:geom.N-0)'+geom.N+1];
        end
        %% Get ir: % removed degrees of freedom
        function ir = get.ir(geom)
            ir = [[1:2,geom.N,geom.N + 1], [1,geom.N + 1] + geom.N + 1 ]';
        end
        %% Get ii: % sort degrees of freedom
        function ii = get.ii(geom)
            [~,ii] = sort([geom.ik; geom.ir]);
        end
        %% Get cCtemp
        function cC = get.cC(geom)
            c1 = zeros(2,geom.N+1);
            c1([1,2*geom.N+2]) = 1; % dirichlet at first and last mesh point
            c2 = geom.D1([1,geom.N+1],:); % Neuman at first and last mesh point
            cC2 = [ c1; c2]; % the constraint matrix
            cC = [cC2, zeros(4,geom.N+1); zeros(2,geom.N+1), c1];
        end
        %% Get cG
        function cG = get.cG(geom)
            %cG = -( diag([1,1,1/2,1/2,1,1]) * geom.cC(:,geom.ir)) \ geom.cC(:,geom.ik); % give-back matrix
            cG = -( geom.cC(:,geom.ir)) \ geom.cC(:,geom.ik); % give-back matrix
            %cG( abs (cG) < 1e-16) = 0;
        end
        %% Get cH
        function cH = get.cH(geom)
            cH =  inv(diag([1,1,1/2,1/2,1,1]) * geom.cC(:,geom.ir) ); % give-back matrix   
            %cH( abs (cH) < 1e-16) = 0;
        end
    end
    
end

