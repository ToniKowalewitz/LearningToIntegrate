%% Build documenation
% Running this script builds all documenation for the package.
% Currently it works by enumerating below all files we want
% to build documentation from/for.
%
% Note that the following setting strips from all files
% package namespaces, e.g., |demo.common.check_convergence|
% results in |check_convergence.html| in a flat build dir.
% Link as |<check_convergence.html DISPLAY_NAME>|.

%%
% Place the build in |doc| subdirectory
docpath = strcat(pwd, '/doc');

%%
% Build pages for functions, classes and other infrastructure

opts = struct('evalCode', false, 'outputDir', docpath);

publish('demo.common.check_convergence', opts);
publish('doc.index', opts);
publish('doc.build', opts);

%%
% Generate API doc for functions and classes; m2html version
% from https://github.com/blechta/m2html is required
opts = struct;
opts.htmlDir = strcat(docpath, '/api');
opts.recursive = 'on';
opts.globalHypertextLinks = 'on';
opts.todo = 'on';
opts.graph = 'on';
%opts.search = 'on';  % NB: This is broken
opts.ignoredDir = {'.git', 'm2html'};
m2html(opts);

%%
% Build pages for demos by running them. Note that this can
% take a while, depending on what demos do. To keep this step
% feasible, demos should only feature small problems which
% can be run in several seconds.

opts = struct;
opts.evalCode = true;
opts.outputDir = docpath;
opts.createThumbnail = false;
opts.figureSnapMethod = 'print';
opts.useNewFigure = false;

publish('demo.cavity', opts);
publish('demo.mt', opts);
publish('demo.helmholtz', opts);
publish('demo.harmonic', opts);
publish('demo.point', opts);
