function Precomp = TownsendPollardPrecompute2018(c, n, nu)
% TownsendPollardPrecompute2018 Precomputes parameters needed for fluid-filled
% fracture analysis following Townsend and Pollard (2018)
%
% This function precomputes various parameters used in the analysis of 
% fluid-filled fractures based on the methodology presented in:
% "Fluid-filled fractures in Earth's lithosphere: Gravitational loading, 
% interpenetration, and stable height of dikes and veins"
% by Townsend and Pollard (2018), Journal of Structural Geology.
%
% Arguments: (input)
%   c   - Crack half-length [m]
%   n   - Number of points for discretization
%   nu  - Poisson's ratio of the material [-]
%
% Arguments: (output)
%   Precomp - Structure containing precomputed parameters:
%     .z      - Array of positions along crack [-c to c]
%     .zmid   - Midpoints between z positions
%     .dz     - Spacing between z points
%     .KIb1   - Precomputed stress intensity factor term 1
%     .KIb2   - Precomputed stress intensity factor term 2
%     .KIb3   - Precomputed stress intensity factor term 3
%     .part2  - Precomputed displacement term 2
%     .part3  - Precomputed displacement term 3
%     .part4  - Precomputed displacement term 4
%     .part5  - Precomputed displacement term 5
%     .part6  - Precomputed displacement term 6
%     .part2G - Precomputed gradient term 2
%     .part3G - Precomputed gradient term 3
%     .part4G - Precomputed gradient term 4
%     .part5G - Precomputed gradient term 5
%     .part6G - Precomputed gradient term 6
%
% Example usage:
%   c = 1;              % Crack half-length
%   n = 50;             % Number of discretization points
%   nu = 0.25;          % Poisson's ratio
%   Precomp = TownsendPollardPrecompute2018(c, n, nu);
%
% References:
%   Townsend, M.R., and Pollard, D.D., 2018, Fluid-filled fractures in
%   Earth's lithosphere: Gravitational loading, interpenetration, and
%   stable height of dikes and veins: Journal of Structural Geology,
%   v. 109, p. 38-54, doi:10.1016/j.jsg.2017.11.015
%
% See also: PollardTractionIntegrationVerticalCrack, 
%          TownsendPollardPartialLinGrad_disp

    % Generate position arrays
    z = linspace(-c, c, n);
    dz = diff(z);
    zmid = z(1:end-1) + dz/2;
    
    % Convert to complex circle coordinates
    VarTheta = acos(zmid./c);
    Z = exp(1i.*VarTheta);
    
    % Initialize arrays
    [part2, part3, part4, part5, part6] = deal(zeros(numel(Z)));
    [part2GRAD, part3GRAD, part4GRAD, part5GRAD, part6GRAD] = deal(zeros(numel(Z)));
    [KIb1, KIb2, KIb3] = deal(zeros(size(Z)));
    
    % Compute parameters for each position
    for j = 1:numel(zmid)
        b1 = z(j+1);
        b2 = z(j);
        
        % Boundary condition coordinates
        Theta1 = acos(b1./c);
        Z1 = exp(1i.*Theta1);
        Theta2 = acos(b2./c);
        Z2 = exp(1i.*Theta2);
        
        % Compute stress intensity factors
        KIb1(j) = 2*log(Z2/Z1);
        KIb2(j) = (Z1-(Z1^-1)-Z2+(Z2^-1));
        KIb3(j) = 0.25*((Z1^2)-(Z1^-2)-(Z2^2)+(Z2^-2));
        
        X = 3-4.*nu;
        
        % Compute displacement terms
        part2(j,:) = 2.*(X.*(Z.^-1)-Z).*log(Z2./Z1);
        part3(j,:) = Z+(Z.^-1)-Z1-(Z1.^-1);
        part4(j,:) = X.*log((Z1-Z)./((Z1.^-1)-Z))-log((Z1-(Z.^-1))./((Z1.^-1)-(Z.^-1)));
        part5(j,:) = Z+(Z.^-1)-Z2-(Z2.^-1);
        part6(j,:) = X.*log((Z2-Z)./((Z2.^-1)-Z))-log((Z2-(Z.^-1))./((Z2.^-1)-(Z.^-1)));
        
        % Compute gradient terms
        part2GRAD(j,:) = 2*(X*(Z.^-2)-(Z.^2)).*log(Z2/Z1)+(X+1)*(Z-(Z.^-1))*(Z1-(Z1.^-1)-Z2+(Z2.^-1));
        part3GRAD(j,:) = (Z-(Z.^-1)).^2-(Z1-(Z1.^-1)).^2;
        part4GRAD(j,:) = X.*log((Z1-Z)./((Z1.^-1)-Z))-log((Z1-(Z.^-1))./((Z1.^-1)-(Z.^-1)));
        part5GRAD(j,:) = (Z-(Z.^-1)).^2-(Z2-(Z2.^-1)).^2;
        part6GRAD(j,:) = X.*log((Z2-Z)./((Z2.^-1)-Z))-log((Z2-(Z.^-1))./((Z2.^-1)-(Z.^-1)));
    end
    
    % Store in output structure
    Precomp = struct('z', z, 'zmid', zmid, 'dz', dz(1), ...
                    'KIb1', KIb1, 'KIb2', KIb2, 'KIb3', KIb3, ...
                    'part2', part2, 'part3', part3, 'part4', part4, ...
                    'part5', part5, 'part6', part6, ...
                    'part2G', part2GRAD, 'part3G', part3GRAD, ...
                    'part4G', part4GRAD, 'part5G', part5GRAD, ...
                    'part6G', part6GRAD);
end