% This script runs simplified version of
% app_dc.checkerboard_example.full_benchmark
% which is possible to fit on CI

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dfw-xx-ref2-n0008', 'krylov-full-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dfw-xx-ref2-n0016', 'krylov-full-woodbury');
save('checkerboard-timings-2dfw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dnw-xx-ref2-n0008', 'krylov-no-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dnw-xx-ref2-n0016', 'krylov-no-woodbury');
save('checkerboard-timings-2dnw.mat', 'timings');

clear('timings');
%timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 1, 0, '3dfw-xx-2x1', 'krylov-full-woodbury');  % Too slow...
timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 0, 0, '3dfw-xx-2x1', 'krylov-full-woodbury');  % Too crude...
save('checkerboard-timings-3dfw.mat', 'timings');

clear('timings');
%timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 1, 0, '3dnw-xx-2x1', 'krylov-no-woodbury');  % Too slow...
timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 0, 0, '3dnw-xx-2x1', 'krylov-no-woodbury');  % Too crude...
save('checkerboard-timings-3dnw.mat', 'timings');
