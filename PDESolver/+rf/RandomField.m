classdef (Abstract) RandomField < handle
	% Class containing model data of random fields
	% and allowing to create and evaluate random samples.

	properties (GetAccess = public, SetAccess = protected)
		rf_type					% string containing name of random field type
		number_realizations		% number counting the number of realizations generated
		evaluation_type			% string containging either 'KL' or 'function', sets the way how RF is evaluated
	end

	methods (Access = public)

		function obj = RandomField(rf_type)
			% Constructor

			obj.rf_type = rf_type;
			obj.number_realizations = 0;
			obj.realization_ = struct;
		end
		
		set_realizations(obj, realizations)
	end

	methods (Access = public, Sealed = true)

		function clear_realizations(obj)
			%clear all saved realizations of RF and reset counter
			obj.realization_ = struct;
			obj.number_realizations = 0;
		end

		function add_realizations(obj, number_new_realizations)
			% append new realizations to already saved ones
			new_realizations = obj.generate_realizations(number_new_realizations);
			if(obj.number_realizations)
				obj.realization_ = [obj.realization_, new_realizations];
			else
				obj.realization_ = new_realizations;
			end
			obj.number_realizations = obj.number_realizations + number_new_realizations;
		end

		function result = evaluate_realization(obj, x, y, c, realization_number, varargin)
			% This function evaluates the saved realization at the given points (x,y) or on some given cell of the triangulation.
			% The parameter realization_number gives the number of the realization that needs to be evaluated.
			if(nargin == 6) %override default settings?
				eval_type = varargin{1};
			elseif(nargin > 6)
				error('Too many input arguments.')
			else
				eval_type = obj.evaluation_type;
			end

			% evaluate RF at given points, choose selected method first
			switch eval_type
				case 'KL'
					result = obj.eval_KL(x, y, realization_number);
				case 'function'
					result = obj.eval_function(x, y, realization_number);
				case 'byCell'
					result = obj.eval_by_cell(c, realization_number);
				otherwise
					error('Evaluation type not known.')
			end
		end

		function plot(obj, realization_number)
			% Plot the random field number realization_number in the unit square
			assert(realization_number <= obj.number_realizations, 'This realization does not exist.')
			x = linspace(0,1,101);
			surf(x,x,obj.evaluate_realization(x',x,1,realization_number)');
			xlabel('x')
			ylabel('y')
		end

	end

	methods (Access = protected)

		realizations = generate_realizations(obj, number_new_realizations)
		result = eval_KL(obj, x, y, realization_number)
		result = eval_function(obj, x, y, realization_number)
		result = eval_by_cell(obj, c, realization_number)

	end
	
	properties (Access = protected)
		realization_
	end
end
