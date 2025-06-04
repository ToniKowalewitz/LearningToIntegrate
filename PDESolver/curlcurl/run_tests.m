import matlab.unittest.TestRunner
import matlab.unittest.plugins.XMLPlugin

% Create runner with both text and XML output
runner = TestRunner.withTextOutput;
plugin = XMLPlugin.producingJUnitFormat('test_result.xml');
runner.addPlugin(plugin);

% Select tests
suite = testsuite('tests', 'IncludeSubPackages', true);

% Run tests
result = runner.run(suite);

% Print results
disp(result.table);

% Exit with return code 0 or 1
exit(any([result.Failed]));
