function [K1Top, K1Base, wf, A] = PollardTractionIntegrationVerticalCrack(tnFunc, c, mu, nu, z, zmid, Precomp)
% PollardTractionIntegrationVerticalCrack Calculates complete set of crack parameters 
% using Townsend and Pollard's (2018) method
%
% This function calculates mode I stress intensity factors at both crack tips,
% crack opening displacement profile, and crack area for a vertical crack under
% specified traction conditions. The calculation uses a discretized approach where
% the crack is divided into segments and contributions from each segment are summed.
%
% The calculation follows the methodology presented in:
% "Fluid-filled fractures in Earth's lithosphere: Gravitational loading, 
% interpenetration, and stable height of dikes and veins"
% by Townsend and Pollard (2018), Journal of Structural Geology.
%
% Arguments: (input)
%   tnFunc  - Function handle for normal traction distribution
%             Must accept a position argument z and return traction value
%   c       - Crack half-length [m]
%   mu      - Shear modulus [Pa]
%   nu      - Poisson's ratio [-]
%   z       - Array of positions along crack [-c to c] [m]
%   zmid    - Midpoints between z positions [m]
%   Precomp - Structure containing precomputed parameters from 
%             TownsendPollardPrecompute2018, must include fields:
%             .KIb1   - Precomputed stress intensity factor term 1
%             .KIb2   - Precomputed stress intensity factor term 2
%             .part2  - Precomputed displacement term 2
%             .part3  - Precomputed displacement term 3
%             .part4  - Precomputed displacement term 4
%             .part5  - Precomputed displacement term 5
%             .part6  - Precomputed displacement term 6
%
% Arguments: (output)
%   K1Top   - Mode I stress intensity factor at top of crack [Pa*m^0.5]
%   K1Base  - Mode I stress intensity factor at base of crack [Pa*m^0.5]
%   wf      - Crack opening displacement profile [m]
%             (Only computed if requested in output arguments)
%   A       - Crack area [m^2]
%             (Only computed if requested in output arguments)
%
% Notes:
%   - The function uses optimized versions of displacement and K calculations
%     (TownsendPollardPartialPressure_KSPEED and TownsendPollardPartialPressure_dispSPEED)
%   - Displacement and area calculations are only performed if requested in the
%     output arguments to optimize performance
%   - The crack opening displacement is interpolated to match the input z positions
%   - Area is calculated using trapezoidal integration of the opening profile
%
% Example usage:
%   c = 1;              % Crack half-length
%   mu = 1e9;           % Shear modulus
%   nu = 0.25;          % Poisson's ratio
%   tnFunc = @(z) 1e6;  % Constant traction
%   n = 50;             % Number of discretization points
%   Precomp = TownsendPollardPrecompute2018(c, n, nu);
%   [K1Top, K1Base, wf, A] = PollardTractionIntegrationVerticalCrack(...
%           tnFunc, c, mu, nu, Precomp.z, Precomp.zmid, Precomp);
%
% References:
%   Townsend, M.R., and Pollard, D.D., 2018, Fluid-filled fractures in
%   Earth's lithosphere: Gravitational loading, interpenetration, and
%   stable height of dikes and veins: Journal of Structural Geology,
%   v. 109, p. 38-54, doi:10.1016/j.jsg.2017.11.015
%
% See also: TownsendPollardPrecompute2018, TownsendPollardPartialPressure_KSPEED,
%          TownsendPollardPartialPressure_dispSPEED

    % Initialize stress intensity factors
    K1Top = 0;
    K1Base = 0;
    
    % Calculate stress intensity factors for each segment
    for j = 1:numel(zmid)
        % Calculate segment contribution using precomputed parameters
        [K1TopSeg, K1BaseSeg] = TownsendPollardPartialPressure_KSPEED(...
            z(j+1), z(j), c, tnFunc((z(j)+z(j+1))/2), ...
            Precomp.KIb1(j), Precomp.KIb2(j));
        
        % Add segment contributions to totals
        K1Top = K1Top + K1TopSeg;
        K1Base = K1Base + K1BaseSeg;
    end
    
    % Only calculate displacement if requested
    if nargout > 2
        % Initialize displacement array
        wfPoll = zeros(size(zmid));
        
        % Calculate displacement for each segment
        for j = 1:numel(zmid)
            % Calculate segment contribution using precomputed parameters
            [~, Dn] = TownsendPollardPartialPressure_dispSPEED(...
                z(j+1), z(j), c, tnFunc((z(j)+z(j+1))/2), zmid, nu, mu, ...
                Precomp.part2(j,:), Precomp.part3(j,:), Precomp.part4(j,:), ...
                Precomp.part5(j,:), Precomp.part6(j,:));
            
            % Add segment contribution to total displacement
            wfPoll = wfPoll + Dn;
        end
        
        % Add zero displacement at crack tips
        wfPoll = [0 wfPoll 0];
        zmid = [-c zmid c];
        
        % Interpolate displacement to match input z positions
        wf = interp1(zmid, wfPoll, z);
        
        % Calculate area if requested
        if nargout == 4
            A = trapz(zmid, wfPoll);
        end
    end
end