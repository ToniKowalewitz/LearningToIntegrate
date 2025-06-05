classdef MaternByCell < rf.RandomField

	methods (Access = public)

		function obj = MaternByCell()
			% Constructor

			% Call Superclass constructor first
			obj = obj@rf.RandomField('MaternByCell');

			% Set model data
			obj.evaluation_type = 'byCell';
		end


		function import_realizations(obj,approx_deg,law,total_sample_size)
			obj.approx_deg_ = approx_deg;
			obj.law_ = law;
			obj.total_sample_size_ = total_sample_size;
			import_string = strcat('data/input/sparse_grid_fields_mesh_C_A',int2str(obj.approx_deg_),'_L_',obj.law_,'_N_',int2str(obj.total_sample_size_),'.csv');
			obj.realization_ = readmatrix(import_string);
			obj.realization_ = obj.realization_(2:end,2:end-2);
			obj.number_realizations = size(obj.realization_,1);	
		end

		function apply_exponential(obj)
			obj.realization_ = exp(obj.realization_);
		end

		function apply_affine(obj, m,n)
			obj.realization_ = m*obj.realization_ + n;
		end

		function plot_on_mesh_by_gmesh(obj, path_to_gmesh_file,realization)
			[mesh,~,~,~,~] = io.read_msh(path_to_gmesh_file, [-1 0:2]);
			mesh.compute_connectivity(mesh.dim, mesh.dim-1);
			meshing.plot_mesh(mesh, 'cell_markers',obj.realization_(realization,:))
		end
		
		function plot_on_mesh_by_table(obj,path_to_vertex_coords, path_to_element_table,realization)
			vertex_coords = readmatrix(path_to_vertex_coords)';
			elements = uint32(readmatrix(path_to_element_table)');
			mesh = meshing.Mesh(2,vertex_coords, elements);
			mesh.compute_connectivity(mesh.dim, mesh.dim-1);
			meshing.plot_mesh(mesh, 'cell_markers',obj.realization_(realization,:))
		end

	end

	methods (Access = protected)

		function result = generate_realizations(~,~)
			result = NaN; %#ok<NASGU>
			error('Can only import realizations from external data.')		
		end

		function result = eval_KL(~, ~, ~, ~)
			% KL expansion is not implemented yet.
			result = NaN; %#ok<NASGU>
			error('KL-expansion not implemented. Please use the other evaluation type')
		end

		function result = eval_function(obj, x, y, realization_number)
			result = NaN; %#ok<NASGU>
			error('Point-wise evaluation not implemented. Please use the other evaluation type')
		end
		
		function result = eval_by_cell(obj,c, realization_number)
			result = obj.realization_(realization_number,c);
		end
	end

	properties (Access = private)
		approx_deg_
		law_
		total_sample_size_
	end
end
