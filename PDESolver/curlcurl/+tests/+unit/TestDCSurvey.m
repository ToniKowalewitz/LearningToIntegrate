classdef TestDCSurvey < matlab.unittest.TestCase

    properties (TestParameter)
        survey_conf = {
            % First column: survey type,
            % Second column: widths,
            % Third column: range of number of electrodes,
            {'wenner', [1, 2, 4, 8, 16], 4:30};
            {'pdp',    [1, 2, 4, 8, 16], 3:30};
        };
    end

    methods (Test)

        function test_linear_independence(t, survey_conf)

            % Extract parameters
            [survey_type, widths, num_electrodes_range] = survey_conf{:};

            % Loop over number of electrodes
            for num_electrodes = num_electrodes_range
                test_linear_independence_(t, survey_type, widths, num_electrodes);
            end

        end

    end

end


function test_linear_independence_(t, survey_type, widths, num_electrodes)

    % Create survey
    [Mtx, Mrx] = app_dc.create_electrode_configuration(survey_type, widths, num_electrodes);
    t.assertEqual(size(Mtx), size(Mrx));
    t.assertTrue(issparse(Mtx) && issparse(Mrx));

    % Find linearly independent measurements
    cols = app_dc.select_linearly_independent_measurements_impl(Mtx, Mrx);

    % Should have full rank
    t.verifyEqual(cols, 1:size(Mtx, 2));

end
