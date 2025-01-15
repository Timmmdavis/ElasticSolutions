function [K1Top, K1Base, wf, A] = PollardTractionIntegrationVerticalCrackWithLinGrads(...
    tnFunc, c, mu, nu, z, zmid, Precomp)
% PollardTractionIntegrationVerticalCrackWithLinGrads Calculates crack parameters 
% incorporating linear stress gradients using Townsend and Pollard's (2018) method
%
% This function extends PollardTractionIntegrationVerticalCrack to handle linear
% gradients in the stress field. It calculates mode I stress intensity factors,
% crack opening displacement, and area by decomposing the problem into three parts:
% 1. Base pressure field
% 2. Linear gradient term
% 3. Pressure pullback correction
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
%             .KIb3   - Precomputed stress intensity factor term 3
%             .part2  - Precomputed displacement term 2
%             .part3  - Precomputed displacement term 3
%             .part4  - Precomputed displacement term 4
%             .part5  - Precomputed displacement term 5
%             .part6  - Precomputed displacement term 6
%             .part2G - Precomputed gradient displacement term 2
%             .part3G - Precomputed gradient displacement term 3
%             .part4G - Precomputed gradient displacement term 4
%             .part5G - Precomputed gradient displacement term 5
%             .part6G - Precomputed gradient displacement term 6
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
%   - The function decomposes the problem into three components:
%     1. Base pressure (Z terms)
%     2. Linear gradient (Y terms)
%     3. Pressure pullback correction (X terms)
%   - Uses optimized versions of displacement and K calculations with 
%     precomputed terms for efficiency
%   - Gradient calculations require additional precomputed terms compared
%     to the basic version
%   - Displacement and area calculations are only performed if requested
%
% Example usage:
%   c = 1;              % Crack half-length
%   mu = 1e9;           % Shear modulus
%   nu = 0.25;          % Poisson's ratio
%   tnFunc = @(z) 1e6;  % Constant traction
%   n = 50;             % Number of discretization points
%   Precomp = TownsendPollardPrecompute2018(c, n, nu);
%   [K1Top, K1Base, wf, A] = PollardTractionIntegrationVerticalCrackWithLinGrads(...
%           tnFunc, c, mu, nu, Precomp.z, Precomp.zmid, Precomp);
%
% References:
%   Townsend, M.R., and Pollard, D.D., 2018, Fluid-filled fractures in
%   Earth's lithosphere: Gravitational loading, interpenetration, and
%   stable height of dikes and veins: Journal of Structural Geology,
%   v. 109, p. 38-54, doi:10.1016/j.jsg.2017.11.015
%
% See also: TownsendPollardPrecompute2018, PollardTractionIntegrationVerticalCrack,
%          TownsendPollardPartialPressure_KSPEED, TownsendPollardPartialLinGrad_KSPEED

    % Evaluate traction and calculate its gradient
    tn = tnFunc(zmid);
    dzmid = gradient(zmid);
    dtn = gradient(tn)./dzmid;
    
    % Initialize stress intensity factors
    K1Top = 0;
    K1Base = 0;
    
    % Calculate stress intensity factors for each segment
    for j = 1:numel(zmid)
        % 1. Pressure pullback correction
        [K1TopX, K1BaseX] = TownsendPollardPartialPressure_KSPEED(...
            z(j+1), z(j), c, -zmid(j)*dtn(j), ...
            Precomp.KIb1(j), Precomp.KIb2(j));
        
        % 2. Linear gradient term
        [K1TopY, K1BaseY] = TownsendPollardPartialLinGrad_KSPEED(...
            z(j+1), z(j), c, dtn(j), ...
            Precomp.KIb1(j)/2, Precomp.KIb2(j), Precomp.KIb3(j));
        
        % 3. Base pressure term
        [K1TopZ, K1BaseZ] = TownsendPollardPartialPressure_KSPEED(...
            z(j+1), z(j), c, tn(j), ...
            Precomp.KIb1(j), Precomp.KIb2(j));
        
        % Sum all contributions
        K1Top = K1Top + K1TopX + K1TopY + K1TopZ;
        K1Base = K1Base + K1BaseX + K1BaseY + K1BaseZ;
    end
    
    % Calculate displacement if requested
    if nargout > 2
        % Initialize displacement array
        wfPoll = zeros(size(zmid));
        
        % Calculate displacement for each segment
        for j = 1:numel(zmid)
            % 1. Pressure pullback correction
            [~, DnX] = TownsendPollardPartialPressure_dispSPEED(...
                z(j+1), z(j), c, -zmid(j)*dtn(j), zmid, nu, mu, ...
                Precomp.part2(j,:), Precomp.part3(j,:), ...
                Precomp.part4(j,:), Precomp.part5(j,:), ...
                Precomp.part6(j,:));
            
            % 2. Linear gradient term
            [~, DnY] = TownsendPollardPartialLinGrad_dispSPEED(...
                z(j+1), z(j), c, dtn(j), zmid, nu, mu, ...
                Precomp.part2G(j,:), Precomp.part3G(j,:), ...
                Precomp.part4G(j,:), Precomp.part5G(j,:), ...
                Precomp.part6G(j,:));
            
            % 3. Base pressure term
            [~, DnZ] = TownsendPollardPartialPressure_dispSPEED(...
                z(j+1), z(j), c, tn(j), zmid, nu, mu, ...
                Precomp.part2(j,:), Precomp.part3(j,:), ...
                Precomp.part4(j,:), Precomp.part5(j,:), ...
                Precomp.part6(j,:));
            
            % Sum all contributions
            wfPoll = sum([DnX + DnY + DnZ; wfPoll]);
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