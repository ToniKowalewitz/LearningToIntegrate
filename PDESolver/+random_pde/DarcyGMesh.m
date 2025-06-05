classdef DarcyGMesh < handle
% This class contains data needed to solve a random Darcy flow equation on the square [a,b]^2, where a and b are defined in init_PDE.
% We empose homogenous Neumann BC on the top and bottom boundary. The conductiviy is given by a random field.
% We use the mixed formulation and discritize with an uniform mesh and RT0-P0 finite elements.

	methods (Access = public)

		function obj = DarcyGMesh(random_field, boundary_left, boundary_right)
		% Constructor
			obj.random_field_ = random_field;
			obj.boundary_conditions_left_ = boundary_left;
			obj.boundary_conditions_right_ = boundary_right;
		end

		function import_gmesh(obj, path_to_gmesh_file)
			[obj.mesh_, ~, ~, ~, ~] = io.read_msh(path_to_gmesh_file, [-1 0:2]);
		end

		function import_mesh_by_table(obj, path_to_vertex_coords, path_to_element_table)
			vertex_coords = readmatrix(path_to_vertex_coords)';
			elements = uint32(readmatrix(path_to_element_table)');
			obj.mesh_ = meshing.Mesh(2, vertex_coords, elements);
		end

		function init_PDE(obj, lower_boundary, upper_boundary)
		% Import mesh from a gmesh file, set up dofmaps and BC.
			obj.mesh_.compute_connectivity(obj.mesh_.dim, obj.mesh_.dim-1);
			
			element_l2 = fe.create_p0_element(obj.mesh_.dim);
			element_hdiv = fe.create_raviart_thomas_element(obj.mesh_.dim,1);

			obj.dofmap_l2_ = assembling.build_dofmap(obj.mesh_, element_l2);
			obj.dofmap_hdiv_ = assembling.build_dofmap(obj.mesh_, element_hdiv);

			% Ensure homogeneous Neumann BC on top and bottom boundary,
			% this equals Dirchlet BC for the (normal component) of the flux.
			tol = 1e-8;
			where_neumann = @(x,b) b && (x(2) < lower_boundary + tol || x(2) > upper_boundary - tol);
			obj.mesh_.compute_boundary_facets();
			obj.marker_ = meshing.mark_entities(obj.mesh_, obj.mesh_.dim-1, [], where_neumann, 1);
			[bc_dofs, bc_values] = assembling.build_dirichlet_dofs(obj.dofmap_hdiv_, obj.marker_, {1}, {@(x) [0, 0]});
			obj.boundary_conditions_dirichlet_ = struct('dofs', bc_dofs, 'values', bc_values);

			% Also mark facets for Dirichlet BC on left and right boundary
			where_dirichlet_left = @(x,b) b && x(1) < lower_boundary + tol;
			obj.marker_ = meshing.mark_entities(obj.mesh_, obj.mesh_.dim-1, obj.marker_, where_dirichlet_left, 2);
			where_dirichlet_right = @(x,b) b && x(1) > upper_boundary - tol;
			obj.marker_ = meshing.mark_entities(obj.mesh_, obj.mesh_.dim-1, obj.marker_, where_dirichlet_right, 3);
		
			% Allocate space for solutions
			obj.reset_solution();
		end

		function solve_PDE(obj, realization)
		% Solve PDE where the coefficient is given by realization of the random field
			% Assemble matrices and rhs vectors
			coefficient = @(x,c) obj.random_field_.evaluate_realization(x(1),x(2),c,realization);
			[M, D] = assembling.assemble_hdiv_operators(obj.dofmap_hdiv_, obj.dofmap_l2_, coefficient, 0);
			rhs_l2 = zeros(obj.dofmap_l2_.dim, 1);

			% Apply BC to system
			rhs_hdiv = assembling.assemble_neumann_source(obj.dofmap_hdiv_, obj.marker_, 2, @(x,n,c) obj.boundary_conditions_left_, 1);
			rhs_hdiv = rhs_hdiv + assembling.assemble_neumann_source(obj.dofmap_hdiv_, obj.marker_, 3, @(x,n,c) obj.boundary_conditions_right_, 1);
			[M, rhs_hdiv] = assembling.apply_dirichlet_bc(M, rhs_hdiv, obj.boundary_conditions_dirichlet_.dofs, obj.boundary_conditions_dirichlet_.values);
			rhs_l2 = rhs_l2 - D(:,obj.boundary_conditions_dirichlet_.dofs)*obj.boundary_conditions_dirichlet_.values; %TODO: this is not dependent on RF -> put to init_PDE + seperate assembly of M, D?
			D(:, obj.boundary_conditions_dirichlet_.dofs) = 0;

			% Build system matrix and rhs
			Z = sparse([], [], [], size(D, 1), size(D, 1));
			A = [M, D.'; D, Z];
			b = [rhs_hdiv; rhs_l2];

			% Solve system
			x = A\b;

			% Split solution into flux and pressure
			obj.solution_(realization).flux = x(1:size(M, 1), 1);
			obj.solution_(realization).pressure = x(size(M, 1)+1:end, 1);
		end

		function plot_solution(obj, realization)
			% Plot solution (pressure) belonging to given realization of random field
			meshing.plot_mesh(obj.mesh_, 'cell_markers', obj.solution_(realization).pressure)
		end
		
		function plot_random_field(obj, realization)
			% Plot realization of the random field
			cells = 1:size(obj.mesh_.cells,2)
			x_coords = (obj.mesh_.vertex_coords(1,obj.mesh_.cells(1,:)) + obj.mesh_.vertex_coords(1,obj.mesh_.cells(2,:)) + obj.mesh_.vertex_coords(1,obj.mesh_.cells(3,:)))/3
			y_coords = (obj.mesh_.vertex_coords(2,obj.mesh_.cells(1,:)) + obj.mesh_.vertex_coords(2,obj.mesh_.cells(2,:)) + obj.mesh_.vertex_coords(2,obj.mesh_.cells(3,:)))/3
			rf_values = obj.random_field_.evaluate_realization(x_coords, y_coords, cells, realization)
			meshing.plot_mesh(obj.mesh_, 'cell_markers', rf_values)
		end

		function qoi = compute_quantity(obj, realization, type, varargin)
		% Compute quantity of interest for a given random solution of the PDE.
		% varargin can be set to true to enable some plots or details useful for
		% debugging or verification. It should be set to false for e.g. MC simulation
		% to avoid spam or reduction of performance.
			if(nargin > 3)
				show_details = varargin{1};
				if(nargin > 4)
					warning('Too many input arguments, ignoring additional parameters.')
				end
			else
				show_details = false;
			end
			% Trajectory number valid?
			assert(realization <= obj.random_field_.number_realizations, 'This realization does not exist!')
			assert(~isempty(obj.solution_(realization).flux) && ~isempty(obj.solution_(realization).pressure), 'Solution has not been computed yet!')
			% Compute quantity of interest of the solution given by certain random field
			switch type
				case 'perm_eff'
					% Comupte effective permeability: This equals the flow through the right boundary.
					% As the normal component of the flux equals the value of the dof of the boundary edges.
					% Those have been marked above with label 3.
					qoi = sum(abs(obj.solution_(realization).flux(obj.marker_ == 3)));
					obj.solution_(realization).flux(obj.marker_ == 3)
				case 'travel_time'
					% Compute travel time from a given starting point to the boundary of the domain.
					% Initialize some quantities first.
					starting_point = [.05; .5];
					travel_time = 0;
					path = zeros(2,2*sqrt(size(obj.mesh_.cells,2)/2));
					path(:,1) = starting_point;
					current_point = starting_point;
					triangle = obj.point_to_cell(starting_point);
					iteration=1;
					% Main loop: we are still inside our computational domain.
					while(triangle)
						% Compute distance of entry point to triangle edges
						% via projection onto normal vector.
						point1 = obj.mesh_.vertex_coords(:,obj.mesh_.cells(1,triangle));
						point2 = obj.mesh_.vertex_coords(:,obj.mesh_.cells(2,triangle));
						point3 = obj.mesh_.vertex_coords(:,obj.mesh_.cells(3,triangle));
						edge1 = point3-point2; normal1 = [-edge1(2);edge1(1)]; normal1 = (-1)^triangle*normal1/norm(normal1);
						edge2 = point1-point3; normal2 = [-edge2(2);edge2(1)]; normal2 = (-1)^triangle*normal2/norm(normal2);
						edge3 = point2-point1; normal3 = [-edge3(2);edge3(1)]; normal3 = (-1)^triangle*normal3/norm(normal3);
						distance1 = (point3 - current_point)'*normal1;
						distance2 = (point1 - current_point)'*normal2;
						distance3 = (point2 - current_point)'*normal3;

						% Compute velocity of flux towards each boundary
						flux = obj.get_contant_flux(realization, triangle);
						velocity1 = flux'*normal1;
						velocity2 = flux'*normal2;
						velocity3 = flux'*normal3;

						% Compute time needed to get to each edge.
						time1 = obj.travel_time_cases(velocity1, distance1);
						time2 = obj.travel_time_cases(velocity2, distance2);
						time3 = obj.travel_time_cases(velocity3, distance3);

						% Find first edge we cross (=leaving edge).
						[time, index] = min([time1 time2 time3]);
						% Update time and path.
						travel_time = travel_time + time;
						current_point = current_point + time*flux;
						path(:,iteration) = current_point;
						iteration = iteration+1;
						% Check if we are leaving our domain.
						if (current_point(1) <= eps || current_point(1) >= 1 - eps || current_point(2) <= eps || current_point(2) >= 1 - eps)
							% We are done: violate condition in while loop.
							triangle = 0;
						else
							% Find triangle we are entering next.
							% Check if we are leaving through a vertex.
							[leaving_through_vertex, vertex_index] = ismember(current_point', [point1, point2, point3]', 'rows');
							if (leaving_through_vertex)
								% Get local index of vertex.
								vertex_index = obj.mesh_.cells(vertex_index, triangle);
								% Find all triangles containing leaving point (except current triangle).
								adjacent_triangles = find(obj.mesh_.cells(1,:) == vertex_index | obj.mesh_.cells(2,:) == vertex_index | obj.mesh_.cells(3,:) == vertex_index);
								adjacent_triangles = setdiff(adjacent_triangles, triangle);

								% Find triangle the flux is pointing into. This
								% can be done by checking if the flux is contained
								% in the cone spanned by the two adjacent edges.
								% This is equivalent to the flux being a linear
								% combination with non-negative coefficients.
								for triangle_index = adjacent_triangles
									other_vertices = setdiff(obj.mesh_.cells(:,triangle_index), vertex_index);
									cone_matrix = obj.mesh_.vertex_coords(:,other_vertices) - current_point;
									coefficients = cone_matrix\flux;
									if(~sum(coefficients<0))
										triangle = triangle_index;
										break;
									end
								end
							else
								% We are leaving through an edge.
								switch index
								case 1
									if (mod(triangle,2))
										% Go from bottom-left to upper-right
										triangle = triangle + 1;
									else
										% Go from bottom to upper
										triangle = triangle + 2*sqrt(size(obj.mesh_.cells,2)/2) - 1;
									end
								case 2
									if (mod(triangle,2))
										% Go from right to left
										triangle = triangle - 1;
									else
										% Go from left to right
										triangle = triangle + 1;
									end
								case 3
									if (mod(triangle,2))
										% Go from upper to bottom
										triangle = triangle - 2*sqrt(size(obj.mesh_.cells,2)/2) + 1;
									else
										% Go from upper-right to bottom-left
										triangle = triangle - 1;
									end
								end
							end
						end
					end
					path(:,iteration:end) = [];
					if (show_details)
						figure(2); clf;
						obj.plot_solution(realization)
						hold on
						plot(path(1,:), path(2,:),'k-')
					end
					qoi = travel_time;
				otherwise
					error('Unkown quantity of interest %s.', type)
			end
		end

		function qoi = solve_and_eval_PDE(obj, realization)
		% Solve PDE where the coefficient is given by realization of the random field
			% Assemble matrices and rhs vectors
			coefficient = @(x,c) obj.random_field_.evaluate_realization(x(1),x(2),c,realization);
			[M, D] = assembling.assemble_hdiv_operators(obj.dofmap_hdiv_, obj.dofmap_l2_, coefficient, 0);
			rhs_l2 = zeros(obj.dofmap_l2_.dim, 1);

			% Apply BC to system
			rhs_hdiv = assembling.assemble_neumann_source(obj.dofmap_hdiv_, obj.marker_, 2, @(x,n,c) obj.boundary_conditions_left_, 1);
			rhs_hdiv = rhs_hdiv + assembling.assemble_neumann_source(obj.dofmap_hdiv_, obj.marker_, 3, @(x,n,c) obj.boundary_conditions_right_, 1);
			[M, rhs_hdiv] = assembling.apply_dirichlet_bc(M, rhs_hdiv, obj.boundary_conditions_dirichlet_.dofs, obj.boundary_conditions_dirichlet_.values);
			rhs_l2 = rhs_l2 - D(:,obj.boundary_conditions_dirichlet_.dofs)*obj.boundary_conditions_dirichlet_.values; %TODO: this is not dependent on RF -> put to init_PDE + seperate assembly of M, D?
			D(:, obj.boundary_conditions_dirichlet_.dofs) = 0;

			% Build system matrix and rhs
			Z = sparse([], [], [], size(D, 1), size(D, 1));
			A = [M, D.'; D, Z];
			b = [rhs_hdiv; rhs_l2];

			% Solve system
			x = A\b;

			% Split solution into flux and pressure
			flux = x(1:size(M, 1), 1);
			pressure = x(size(M, 1)+1:end, 1);
			qoi = sum(abs(flux(obj.marker_ == 3)));
		end

		function reset_solution(obj)
			% Clear all solutions and create empty struct array for solution vectors, 
			% depending of number of existing realizations of the random field
			obj.solution_ = struct;
			if(obj.random_field_.number_realizations)
				obj.solution_(obj.random_field_.number_realizations).flux = [];
				obj.solution_(obj.random_field_.number_realizations).pressure = [];
			end
		end
	end

	properties (Access = private)
		mesh_
		dofmap_hdiv_
		dofmap_l2_
		boundary_conditions_dirichlet_
		random_field_
		boundary_conditions_left_
		boundary_conditions_right_
		solution_
		marker_
	end

	methods (Access = private)

		function time = travel_time_cases(~, velocity, distance)
		% Return time needed to leave triangle or return Inf
		% if flow is in opposite direction.
			% Take care of round-off errors
			distance = distance * (abs(distance) > eps);
			if (velocity > eps) % Is eps better than 0? Change due to misbehavior with flux parallel to edge.
				time = distance/velocity;
			else
				time = Inf;
			end
		end

		function flux = get_contant_flux(obj, realization, triangle)
		% Compute and return the (for homogenous equations!) contant
		% flow of a random solution on the given triangle.
			% Compute basis on reference element.
			evaluation_point_reference = [1/3 1/3];
			basis_hdiv = obj.dofmap_hdiv_.element.tabulate_basis(evaluation_point_reference)';
			% Get indizes and coordinates of the triangle
			nodes = obj.mesh_.cells(:,triangle);
			edges = obj.dofmap_hdiv_.cell_dofs(:,triangle);
			coords = obj.mesh_.vertex_coords(:,nodes);
			% Get solution on reference element
			sol_ref = basis_hdiv * obj.solution_(realization).flux(edges);
			% Apply Piola transform
			BK = coords(:,1:2) - coords(:,3);
			detBK = det(BK);
			flux = -BK*sol_ref / detBK;
		end
	end
end
