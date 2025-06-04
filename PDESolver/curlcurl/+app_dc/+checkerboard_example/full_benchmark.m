% This script runs full DC checkerboards benchmark
% to stress the Woodbury-MINRES-AMG solver. This takes
% too much time and memory and such is not tested on CI.
%
% Postprocess by other script in this directory as
% appropriate.

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dfw-ref2-n0004', 'krylov-full-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dfw-ref2-n0008', 'krylov-full-woodbury');
timings(3) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dfw-ref2-n0016', 'krylov-full-woodbury');
timings(4) = app_dc.checkerboard_example.solve(2,   32, 1, 2, '2dfw-ref2-n0032', 'krylov-full-woodbury');
timings(5) = app_dc.checkerboard_example.solve(2,   64, 1, 2, '2dfw-ref2-n0064', 'krylov-full-woodbury');
timings(6) = app_dc.checkerboard_example.solve(2,  128, 1, 2, '2dfw-ref2-n0128', 'krylov-full-woodbury');
timings(7) = app_dc.checkerboard_example.solve(2,  256, 1, 2, '2dfw-ref2-n0256', 'krylov-full-woodbury');
save('checkerboard-timings-2dfw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(2,    4, 1, 2, '2dnw-ref2-n0004', 'krylov-no-woodbury');
timings(2) = app_dc.checkerboard_example.solve(2,    8, 1, 2, '2dnw-ref2-n0008', 'krylov-no-woodbury');
timings(3) = app_dc.checkerboard_example.solve(2,   16, 1, 2, '2dnw-ref2-n0016', 'krylov-no-woodbury');
timings(4) = app_dc.checkerboard_example.solve(2,   32, 1, 2, '2dnw-ref2-n0032', 'krylov-no-woodbury');
timings(5) = app_dc.checkerboard_example.solve(2,   64, 1, 2, '2dnw-ref2-n0064', 'krylov-no-woodbury');
save('checkerboard-timings-2dnw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 2, 1, '3dfw-2x1', 'krylov-full-woodbury');
timings(2) = app_dc.checkerboard_example.solve(3, [2, 2], 2, 1, '3dfw-2x2', 'krylov-full-woodbury');
timings(3) = app_dc.checkerboard_example.solve(3, [2, 3], 2, 1, '3dfw-2x3', 'krylov-full-woodbury');
timings(4) = app_dc.checkerboard_example.solve(3, [2, 4], 2, 1, '3dfw-2x4', 'krylov-full-woodbury');
save('checkerboard-timings-3dfw.mat', 'timings');

clear('timings');
timings(1) = app_dc.checkerboard_example.solve(3, [2, 1], 2, 1, '3dnw-2x1', 'krylov-no-woodbury');
timings(2) = app_dc.checkerboard_example.solve(3, [2, 2], 2, 1, '3dnw-2x2', 'krylov-no-woodbury');
timings(3) = app_dc.checkerboard_example.solve(3, [2, 3], 2, 1, '3dnw-2x3', 'krylov-no-woodbury');
save('checkerboard-timings-3dnw.mat', 'timings');
