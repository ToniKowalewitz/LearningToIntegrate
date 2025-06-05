path_to_vertex_coords = 'meshes/mesh-48-nodes.txt';
path_to_element_table = 'meshes/mesh-48-elements.txt';
%path_to_vertex_coords = 'meshes/mesh-798-nodes.txt';
%path_to_element_table = 'meshes/mesh-798-elements.txt';
%path_to_vertex_coords = 'meshes/mesh-8510-nodes.txt';
%path_to_element_table = 'meshes/mesh-8510-elements.txt';
grid_type='element_table';
parallel=true;
randomField = rf.MaternByCell();
dirichlet_left = 1;
dirichlet_right = 0;
solver = random_pde.DarcyGMesh(randomField, dirichlet_left, dirichlet_right);
evaluater = uq.EvaluatePredefinedSamples(solver,randomField,parallel,grid_type,path_to_vertex_coords,path_to_element_table);

for iter_deg = 2:2
	switch iter_deg
	case 1
		approx_deg = 2
	case 2
		approx_deg = 3
	end
	for iter_law = 1:4
		switch iter_law
		case 1
			law = 'norm'
		case 2
			law = 'pois'
		case 3
			law = 'gamma'
		case 4
			law = 'bigamma'
		end
		for iter_sample = 1:4
			switch iter_sample  %sample size used for training
				case 0
					sample_size=10
				case 1
					sample_size=100
				case 2
					sample_size=1000
				case 3
					sample_size=10000
				case 4
					sample_size=100000
			end
			evaluater.import_RF(approx_deg, law, sample_size);
			evaluater.solve();
			evaluater.export_results();
		end
	end
end	
