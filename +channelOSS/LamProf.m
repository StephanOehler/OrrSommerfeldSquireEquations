classdef LamProf
    % Generates a laminar mean profile (y for wall normal location)
    % ViscProf and LamProf use y as wall-normal component, every other
    % code uses z as wall-normal component
    
    properties
        %% Channel parameters
        y = [];% y-grid points
        type = 'channel'
        % 'channel','couette'
        h = 1; % wall height
    end
    %% Dependent properties
    properties (Dependent, Access = public)
        %% Derivatives
        U % velocity
        dU
        ddU
    end
    properties (Dependent, Access = protected)
        %% vector of ones
        onevec
    end
    %%
    methods
        %% Setup function - default overwrite
        function thisLamProf = LamProf(y,type,h)
            if nargin < 1 || nargin > 3
            else
                thisLamProf.y = y;
                if nargin > 1
                    thisLamProf.type = type;
                end
                if nargin > 2
                    thisLamProf.h = h;
                end
            end
        end
        %% Restrict type
        function thisLamProf = set.type(thisLamProf,type)
            if  any(strcmp(type,{'channel','couette'}))
                thisLamProf.type = type;
            else
                error('type must be either: ''channel'',''couette''')
            end
        end
        %% Get U: Velocity Profile
        function U = get.U(LamProf)
            if strcmp(LamProf.type,'channel')
                U = LamProf.h - LamProf.y.^2;
            else
                U = LamProf.y;
            end
        end
        %% Get dU: derrivative of Velocity Profile
        function dU = get.dU(LamProf)
            if strcmp(LamProf.type,'channel')
                dU  = - 2 * LamProf.y;
            else
                dU  = LamProf.onevec;
            end
        end
        %% Get ddU: 2nd derrivative of Velocity Profile
        function ddU = get.ddU(LamProf)
            if strcmp(LamProf.type,'channel')
                ddU  = -2 * LamProf.onevec;
            else
                ddU  = 0 * LamProf.onevec;
            end
        end
        %% Get nx: number of grid points
        function onevec = get.onevec(LamProf)
            onevec = ones(length(LamProf.y ), 1);
        end
    end
end

