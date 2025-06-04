%#ok<*INUSD>

function tests = test_demos
    % Run all function tests in this file
    tests = functiontests(localfunctions);
end

function test_cavity(t)
    evalc('demo.cavity');
end

function test_darcy(t)
    evalc('demo.darcy');
end

function test_darcy2(t)
    evalc('demo.darcy2');
end

function test_eit(t)
    evalc('demo.eit');
end

function test_harmonic(t)
    evalc('demo.harmonic');
end

function test_helmholtz(t)
    evalc('demo.helmholtz');
end

function test_mt(t)
    evalc('demo.mt');
end

function test_plot_rba_error(t)
    evalc('demo.plot_rba_error');
end

function test_point(t)
    evalc('demo.point');
end

function test_robin(t)
    evalc('demo.robin');
end

function test_tem(t)
    evalc('demo.tem');
end

function test_mt_inversion(t)
    evalc('app_mt.demo_inversion');
end

function test_mt_layered_space(t)
    evalc('app_mt.demo_layered_space');
end

function test_dc_checkerboard(t)
    evalc('dc_checkerboard_example');
    evalc('app_dc.checkerboard_example.postprocess');
    evalc('app_dc.checkerboard_example.postprocess2');
end
