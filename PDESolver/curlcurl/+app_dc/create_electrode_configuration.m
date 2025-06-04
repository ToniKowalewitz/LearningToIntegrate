function [Mtx, Mrx] = create_electrode_configuration(type, widths, num_electrodes)
    % Create 1D DC survey
    %
    % INPUT PARAMETERS
    %   type           ... char array, currenty 'wenner' or 'pdp'
    %   widths         ... vector of config widths, e.g., [2, 4, 8]
    %   num_electrodes ... number of electrodes
    %
    % OUTPUT PARAMETERS
    %   Mtx, Mrx ... matrix [num_electrodes, num_observations], the two matrices
    %                define observations (voltage per unit current) in terms of
    %                potential point values at tx electrodes for point sources
    %                at rx electrodes

    assert(ischar(type));
    assert(isvector(widths) && ispositiveint(widths));
    assert(isscalar(num_electrodes) && ispositiveint(num_electrodes));

    switch lower(type)
    case 'wenner'
        [Mtx, Mrx] = create_electrode_configuration_wenner(widths, num_electrodes);
    case 'pdp'
        [Mtx, Mrx] = create_electrode_configuration_pdp(widths, num_electrodes);
    otherwise
        error('Electrode configuration type "%s" unknown!', type);
    end

    assert(issparse(Mtx) && issparse(Mrx));

end


function [Mtx, Mrx] = create_electrode_configuration_wenner(widths, num_electrodes)
    Mtx = spalloc(num_electrodes, 0, 0);
    Mrx = spalloc(num_electrodes, 0, 0);
    for width = widths(:).'
        [Mtx_, Mrx_] = create_electrode_configuration_wenner_(width, num_electrodes);
        Mtx = [Mtx, Mtx_];  %#ok<AGROW>
        Mrx = [Mrx, Mrx_];  %#ok<AGROW>
    end
end


function [Mtx, Mrx] = create_electrode_configuration_wenner_(width, num_electrodes)
    num_configs = max(num_electrodes - 3*width, 0);
    Mtx = spalloc(num_electrodes, num_configs, 2*num_configs);
    Mrx = spalloc(num_electrodes, num_configs, 2*num_configs);
    ia = 0*width;
    im = 1*width;
    in = 2*width;
    ib = 3*width;
    for k = 1:num_configs
        Mtx(ia+k, k) = +1;  %#ok<SPRIX>
        Mtx(ib+k, k) = -1;  %#ok<SPRIX>
        Mrx(im+k, k) = +1;  %#ok<SPRIX>
        Mrx(in+k, k) = -1;  %#ok<SPRIX>
    end
end


function [Mtx, Mrx] = create_electrode_configuration_pdp(widths, num_electrodes)
    Mtx = spalloc(num_electrodes, 0, 0);
    Mrx = spalloc(num_electrodes, 0, 0);
    for width = widths(:).'
        [Mtx_, Mrx_] = create_electrode_configuration_pdp_(width, num_electrodes);
        Mtx = [Mtx, Mtx_];  %#ok<AGROW>
        Mrx = [Mrx, Mrx_];  %#ok<AGROW>
    end
end


function [Mtx, Mrx] = create_electrode_configuration_pdp_(width, num_electrodes)
    num_configs = max(num_electrodes - 2*width, 0);
    Mtx = spalloc(num_electrodes, 2*num_configs, 1*num_configs);
    Mrx = spalloc(num_electrodes, 2*num_configs, 2*num_configs);
    ia = 0*width;
    im = 1*width;
    in = 2*width;
    for k = 1:num_configs
        Mtx(ia+k, k) = +1;  %#ok<SPRIX>
        Mrx(im+k, k) = +1;  %#ok<SPRIX>
        Mrx(in+k, k) = -1;  %#ok<SPRIX>
    end
    for k = 1:num_configs
        Mtx(num_electrodes+1-ia-k, num_configs+k) = +1;  %#ok<SPRIX>
        Mrx(num_electrodes+1-im-k, num_configs+k) = +1;  %#ok<SPRIX>
        Mrx(num_electrodes+1-in-k, num_configs+k) = -1;  %#ok<SPRIX>
    end
end


function flag = ispositiveint(mat)
    % Is matrix of positive integers (possibly of floating point type)
    flag = isint(mat) && isreal(mat) && all(mat > 0, 'all');
end


function flag = isint(mat)
    % Is matrix of integers (possibly of floating point type)
    flag = all(floor(mat) == ceil(mat), 'all');
end
