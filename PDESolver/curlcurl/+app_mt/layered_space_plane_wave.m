function E = layered_space_plane_wave(f, mu, rho, d, z)
    % Compute electric field for plane-wave excitation of 1D layered space
    %
    % INPUT PARAMETERS
    %
    %   f   ... frequency [Hz]
    %   mu  ... vacuum permeability [H/m]
    %   rho ... vector [num_layers] of electrical resistivities [Ohm*m]
    %   d   ... vector [num_layers+1] depths s(z-coordinate) of the layers [m]
    %   z   ... z-coordinates where to evaluate electrical field [m]
    %
    % OUTPUT PARAMETERS
    %
    %   E ... electrical field [V/m] at points z
    %
    % NOTES
    %
    %   rho is a vector of length matching the number of layers;
    %   d are depths (z-coordinate) of layer interfaces including
    %   -inf and +inf; hence d is a vector of length number of
    %   layers plus one.
    %
    %   z is vector z-coordinates where to evaluate electrical field;
    %   the function is vectorized and returns E of the same shape.
    %
    %   E is normalized such that magnetic field H = 1 A/m at z = d(2).

    % Check input
    assert(isscalar(f) && isreal(f) && f > 0, 'Expecting positive ''f''');
    assert(isscalar(mu) && isreal(mu) && mu > 0, 'Expecting positive ''mu''');
    assert(numel(rho)+1 == numel(d), ...
           ['Expecting matching number of layers in ''rho'' and ''d'', ', ...
            'i.e., numel(d) == numel(rho) + 1']);
    assert(d(1) == -inf, 'First depth should be -inf');
    assert(d(end) == +inf, 'Last depth should be +inf');
    assert(all(d(:) == sort(d(:))), 'Depths should be in increasing order');
    assert(numel(d) >= 3, 'Need at least two layers');

    % Convert arguments
    iwm = 1i * 2 * pi * f * mu;  % FIXME: Would getE1dMT() work with the eps*omega^2*mu part?
    sigair = 1/rho(1);
    layer_resistivities = rho(2:end);
    layer_thicknesses = d(3:end-1) - d(2:end-2);
    z = z - d(2);

    % Actual computation
    E = getE1dMT(iwm, sigair, layer_resistivities, layer_thicknesses, z);

    % Return E in array of the same shape as z
    E = reshape(E, size(z));
end


function E = getE1dMT(iwm,sigair,rho,d,zi)
% E = getE1DMT(F,RHO,D,Z)
%
% Compute horizontal electric field E within and above a layered half-space.
% Plane wave excitation is assumed with magnetic field H = 1 A/m at z >= 0.
%
% Input:
%
% IWM: Constant 1i*2*pi*f*mu
% RHO: NL vector of layer resistivites within conducting half-space
% D:   vector of layer thicknesses, length NL - 1
% Z:   NZ Vector of z coordinates where electric fields have to be computed.
%      Note that z axis is pointing downward.
%
% Output:
%
% E:   NZ vector of complex electric fields at depths Z.
% Note that E is normalized such that magnetic field H = 1 A/m at z = 0.
%
% Ralph-Uwe Boerner (2009)

assert(nargin == 5, 'getE1DMT requires four input arguments.');

nl = length(rho);
rho = rho(:);
d = d(:);
% coordinates of interfaces
h = [0; cumsum(d)];
% Constants
alpha = complex(zeros(nl, 1));
b = alpha;
aa = complex(zeros(nl - 1, 1));
nz = length(zi);
E = complex(zeros(nz, 1));

alpha = sqrt(iwm ./ rho);
if nl == 1
    c1 = iwm / alpha;
else
    alphad = alpha(1:(nl - 1)) .* d;
    talphad = tanh(alphad);
    % Compute admittance at surface of layer recursively
    b(nl) = alpha(nl);
    for nn = nl - 1:-1:1
        b(nn) = alpha(nn) * (b(nn + 1) + alpha(nn) * talphad(nn)) ./...
            (alpha(nn) + b(nn + 1) * talphad(nn));
    end
    % Impedance
    c1 = iwm ./ b(1);
    % Continuation from boundary to boundary
    for nn = 1:nl - 1
        aa(nn) = (b(nn) + alpha(nn)) / (b(nn + 1) + alpha(nn)) * ...
            exp( - alpha(nn) * d(nn));
    end
end

for ii = 1:nz
    z = zi(ii);
    if z >= 0
        if nl == 1
            a = exp(-alpha(nl) * z);
        else
            ind = find(z >= h, 1, 'last');
            if ind < nl
                a = prod(aa(1:ind - 1)) * 0.5 * (1 + b(ind) / alpha(ind)) * ...
                    (exp( - alpha(ind) * (z - h(ind))) - ...
                    (b(ind + 1) - alpha(ind)) / (b(ind + 1) + alpha(ind)) * ...
                    exp( - alpha(ind) * (d(ind) + h(ind + 1) - z)));
            else
                a = prod(aa) * exp( - alpha(ind) .* (z - h(ind)));
            end
        end
    else
        k0 = sqrt(iwm * sigair);
        pr = (c1 - iwm/k0) / (c1 + iwm/k0);
        ar = k0 * z;
        a = exp(-ar) * (1 + pr * exp(2 * ar)) * 1 ./ (1 + pr);
    end
    E(ii) = a * c1;
end

end
