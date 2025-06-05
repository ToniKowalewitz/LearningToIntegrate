classdef EvaluatePredefinedSamples  < handle

	properties (GetAccess = public, SetAccess = protected)
		solution
	end

	properties (Access = protected)
		solver_
		random_field_
		approx_deg_
		law_
		sample_size_
		parallel_
		qoi_type_
	end

	methods (Access = public)
		function obj = EvaluatePredefinedSamples(solver, random_field,parallel,grid_type,varargin)
			% Constructor
			if(nargin > 7)
				error('too many input arguments')
			end
			obj.solver_ = solver;
			obj.parallel_ = parallel;
			obj.qoi_type_ = 'perm_eff';
			obj.random_field_ = random_field;

			switch grid_type
			case 'gmesh'
				path_to_gmesh = varargin{1};
				obj.solver_.import_gmesh(path_to_gmesh);
			case 'element_table'
				path_to_vertex_coords = varargin{1};
				path_to_element_table = varargin{2};
				obj.solver_.import_mesh_by_table(path_to_vertex_coords,path_to_element_table);
			otherwise
				error('wrong grid type')
			end
			obj.solver_.init_PDE(-1,1);

			% Supress warning (stored connectivity matrix for mesh)
			warning('off', 'Mesh:ExtraConnStored')
		end
		
		function import_RF(obj,approx_deg, law, sample_size)
			obj.approx_deg_ = approx_deg;
			obj.law_ = law;
			obj.sample_size_ = sample_size;
			obj.random_field_.import_realizations(approx_deg,law,sample_size);
			switch law
			case 'norm'
				obj.random_field_.apply_affine(1e3,0);
				obj.random_field_.apply_exponential();
			case 'pois'
				obj.random_field_.apply_affine(1e3,1e-2);
				obj.random_field_.apply_exponential();
			case 'gamma'
				obj.random_field_.apply_affine(1e3,1e-2);
				obj.random_field_.apply_exponential();
			case 'bigamma'
				obj.random_field_.apply_affine(1e3,0);
				obj.random_field_.apply_exponential();
			otherwise
				error('law not found')
			end
		end

		function solve(obj)
			switch obj.approx_deg_
			case 2		% sparse grids - 9 modes
				nsample = 1378;
%				nsample = 32658;
			case 3		% sparse grids - 25 modes
%				nsample = 1353;
				nsample = 23554;
			end
%			nsample = obj.sample_size_  % sample size for MC
			solution = zeros(nsample,1);
			if(obj.parallel_)
				parfor sample = 1:nsample
					solution(sample) = obj.solver_.solve_and_eval_PDE(sample);
				end
			else
				for sample = 1:nsample
					solution(sample) = obj.solver_.solve_and_eval_PDE(sample);
				end
			end
			obj.solution = solution;
		end

		function result = get_results(obj)
			result = obj.solution;	
		end

		function export_results(obj)
			export_string = strcat('data/output/mesh_C_A',int2str(obj.approx_deg_),'_L_',obj.law_,'_N_',int2str(obj.sample_size_),'.txt');
			file = fopen(export_string,'w');
			fprintf(file,'%16.15e\n',obj.solution);
			fclose(file);
		end
	end

end
