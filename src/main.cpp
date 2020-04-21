//Tests
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <thread>
#include <ctime>
#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#endif

#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include "element.hpp"
#endif

#ifndef CONCENTRATION_MODEL_HPP
#define CONCENTRATION_MODEL_HPP
#include "concentration_model.hpp"
#endif

#ifndef PEAR_FUNCTORS_HPP
#define PEAR_FUNCTORS_HPP
#include "pear_functors.hpp"
#endif

double TEMP = 25 + 273.15;
size_t NUM_NODES_G = 529;
double CON_O2 = 20.8/100.0;
double CON_CO2 = 0.4/100.0;
std::string FILEPATH = "../Input/new_pear_2.msh";

// Element values configuation
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::SIGMA_UR = 2.8e-10;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::SIGMA_UZ = 1.10e-9;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::SIGMA_VR = 2.32e-9;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::SIGMA_VZ = 6.97e-9;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::RESP_Q = 0.97;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::MAX_FERM_CO2 = 1.61e-4*exp(
			(56700/8.314)*(1/293.15 - 1/TEMP));
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::K_MFU = 0.1149;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::K_MU = 0.4103;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::K_MV = 27.2438;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::V_MU = 2.39e-4*exp(
			(80200/8.314)*(1.0/293.15 - 1.0/TEMP));
template <typename P>
typename FEM_module::Element<P>::node_t 
	FEM_module::ElementTriangular<P>::NUM_ELM = 0;
template <typename P>
typename FEM_module::Element<P>::node_t 
	FEM_module::ElementTriangular<P>::NUM_NODES = NUM_NODES_G;

template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::C_U_AMB = 101300.0*CON_O2/(8.314*TEMP);
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementTriangular<P>::C_V_AMB = 101300.0*CON_CO2/(8.314*TEMP);

// boundary values configuration
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementBoundary<P>::C_U_AMB = 101300.0*CON_O2/(8.314*TEMP);
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementBoundary<P>::C_V_AMB = 101300.0*CON_CO2/(8.314*TEMP);
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementBoundary<P>::RHO_U = 7e-7;
template <typename P>
typename FEM_module::Element<P>::precision_t 
	FEM_module::ElementBoundary<P>::RHO_V = 7.5e-7;
template <typename P>
typename FEM_module::Element<P>::node_t 
	FEM_module::ElementBoundary<P>::NUM_NODES = NUM_NODES_G;

int test_importer(){

	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	const std::vector<std::vector<double>> coords = 
		mesh_importer.node_matrix();
	std::vector<size_t> nodes_from_element;
	std::cout<<"Node Num - x - y"<<std::endl;
	size_t count = 1;
    for (auto elem : coords) {
		std::cout<<"Node "<<count<<" - "<<elem[0]<<" - "<<elem[1]<<std::endl;
		count++;
    }
	std::cout<<std::endl;
	for (auto elem : mesh_importer.element_matrix()){
		FEM_module::ElementTriangular<double> temp_element(coords, elem);
		nodes_from_element = temp_element.nodes();
		for (int i = 0; i < 3; i++){
			std::cout<<nodes_from_element[i]<<" ";
		}
		std::cout<<std::endl;
	}
	return EXIT_SUCCESS;
}

int test_pure_diffusion(){
	std::vector<double> interior_point{0.01, 0.06};
	gsl_vector* sol_c;
	gsl_vector* f_vec_mult;
	double c_u_a;
	double c_v_a;
	double condition = 0;
	double f_norm = 0;
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	sol_c = gsl_vector_alloc((size_t)(2*model.number_nodes()));
	f_vec_mult = gsl_vector_calloc((size_t)(2*model.number_nodes()));
	while(condition < 1e-1){
		c_u_a = (double)rand()/(double)RAND_MAX;
		c_v_a = (double)rand()/(double)RAND_MAX;
		FEM_module::ElementBoundary<double>::C_U_AMB = c_u_a;
		FEM_module::ElementBoundary<double>::C_V_AMB = c_v_a;
		model.generate_f_vector();
		for (size_t idx = 0; idx < model.number_nodes(); idx++){
			gsl_vector_set(sol_c, idx, c_u_a);
			gsl_vector_set(sol_c, idx + model.number_nodes(), c_v_a);
		}
		gsl_spblas_dgemv(CblasNoTrans, 1.0, model.stiffness_matrix(), sol_c, 
				0.0, f_vec_mult);
		gsl_vector_add(f_vec_mult, model.f_vector());
		condition = 0;
		f_norm = 0;
		for (size_t idx = 0; idx < model.number_nodes(); idx++){
			f_norm += pow(gsl_vector_get(model.f_vector(), idx), 2) +
				pow(gsl_vector_get(model.f_vector(), 
							idx + model.number_nodes()),2);
			condition += pow(gsl_vector_get(f_vec_mult, idx), 2) + 
				pow(gsl_vector_get(f_vec_mult, idx + 
							model.number_nodes()), 2);
		}
		condition = sqrt(condition/f_norm);
		std::cout<<condition<<std::endl;
	}
	return EXIT_SUCCESS;
};

int test_constant_resp(){
	gsl_vector* helper;
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	model.generate_f_vector();
	helper = gsl_vector_alloc((size_t)(2*model.number_nodes()));
	gsl_vector_set_all(helper, 3e-18);
	gsl_vector_add(helper, model.f_vector());
	gsl_vector_scale(helper, -1.0);
	model.set_f_vector(helper);
	FEM_module::write_matrix_to_file(model.stiffness_matrix(), "../output/stiff");
	FEM_module::write_vector_to_file(model.f_vector(), "../output/f_vector");
	model.solve_linear_model_LU();
	FEM_module::write_vector_to_file(model.coefficients(), "../output/initial_coeff");
	return EXIT_SUCCESS;
}

int test_initial_cond(){
	std::vector<double> interior_point{0.01, 0.06};
	std::string filename;
	filename = "../output/matlab_mesh_6";
	NUM_NODES_G = 8034;
	FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
	FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
	FEM_module::ImporterText<double> mesh_importer(filename);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	
	model.generate_f_vector();
	model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
				FEM_module::ElementTriangular<double>::C_V_AMB);
	model.add_linear_approx_to_stiffness();
	model.add_linear_approx_to_f_vector();
	model.solve_linear_model();
	FEM_module::write_vector_to_file(model.coefficients(), "../output/initial_coeff");
	return EXIT_SUCCESS;
}

int test_timing_linear_model(){
	std::ofstream myfile;
	myfile.open("../output/lin_mod_runtime", std::ios::out);
	std::vector<double> interior_point{0.01, 0.06};
	std::vector<size_t> nodes_v = {15, 40, 193, 285, 663, 3988, 8034};
	std::string filename;
	time_t time1;
	time_t time2;
	double final_time;
	for (size_t node : nodes_v){
		filename = "../output/matlab_mesh_" + std::to_string(node);
		NUM_NODES_G = node;
		FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementTriangular<double>::NUM_ELM = 0;
		FEM_module::ImporterText<double> mesh_importer(filename);
		mesh_importer.process_file();
		FEM_module::ConcentrationModel<double> model(mesh_importer, 
				interior_point);
		
		model.generate_f_vector();
		model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
					FEM_module::ElementTriangular<double>::C_V_AMB);
		model.add_linear_respiration();
		final_time = 0;
		for (int it = 0; it < 25; it++){
			time(&time1);
			model.solve_linear_model();
			time(&time2);
			if (it > 4){
				final_time += (double)(time2 - time1);
			}
		}
		myfile<<node<<" "<<final_time/20<<std::endl;
	}
	myfile.close();
	return EXIT_SUCCESS;
}

int test_concentration_model_sparse_nonlinear_solver(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_sparse_nonlinear_model();
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}

int test_concentration_model_nonlinear(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.generate_f_vector();
	model.generate_stiffness_matrix();
	model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
				FEM_module::ElementTriangular<double>::C_V_AMB);
	model.generate_initial_codition();
	model.solve_nonlinear_model();
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}

int test_concentration_model_stepped_nonlinear(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_stepped_nonlinear_model(14);
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}


int test_linear_resp(){
	std::vector<double> interior_point{0.01, 0.06};
	std::vector<size_t> nodes_v = {15, 40, 193, 285, 663, 3988, 8034};
	std::string filename;
	int count = 1;
	for (size_t node : nodes_v){
		filename = "../output/matlab_mesh_" + std::to_string(node);
		NUM_NODES_G = node;
		FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementTriangular<double>::NUM_ELM = 0;
		FEM_module::ImporterText<double> mesh_importer(filename);
		mesh_importer.process_file();
		FEM_module::ConcentrationModel<double> model(mesh_importer, 
				interior_point);
		
		model.generate_f_vector();
		model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
					FEM_module::ElementTriangular<double>::C_V_AMB);
		model.add_linear_respiration();
		model.solve_linear_model();
		model.write_elements_to_file("../output/elements_" + std::to_string(count));
		model.write_boundaries_to_file("../output/boundaries_" + std::to_string(count));
		model.write_coordinates_to_file("../output/coords_" + std::to_string(count));
		FEM_module::write_vector_to_file(model.coefficients(), 
				"../output/linear_resp_coeff_" + std::to_string(count));
		count++;
	}
	return EXIT_SUCCESS;
}

int test_new_importer(){;
	NUM_NODES_G = 285;
	FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
	FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterText<double> mesh_importer("../output/matlab_mesh_285");
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.generate_f_vector();
	model.generate_stiffness_matrix();
	model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
				FEM_module::ElementTriangular<double>::C_V_AMB);
	model.generate_initial_codition();
	FEM_module::write_vector_to_file(model.coefficients(), 
		"../output/initial_coeff");
	
	model.solve_sparse_nonlinear_model();
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}

int test_linear_solver_analytical(){
	std::ofstream myfile;
	myfile.open("mesh_convergence", std::ios::out);
	FEM_module::ElementTriangular<double>::SIGMA_UZ = 
		FEM_module::ElementTriangular<double>::SIGMA_UR;
	FEM_module::ElementTriangular<double>::SIGMA_VZ = 
		FEM_module::ElementTriangular<double>::SIGMA_VR;
	std::vector<size_t> nodes_v = {13, 21, 64, 146, 367, 802, 4729};
	std::vector<double> interior_point{0.01, 0.06};
	double rel_err_u = 0;
	double rel_err_v = 0;
	double norm_u = 0;
	double norm_v = 0;
	double radius = 0.06;
	double c0_u = 2.0;
	double c2_u = (c0_u - FEM_module::ElementTriangular<double>::C_U_AMB)/
		(-2*FEM_module::ElementTriangular<double>::SIGMA_UR*radius/
		 FEM_module::ElementBoundary<double>::RHO_U - pow(radius, 2));
	double c0_v = 25.0;
	double c2_v = (c0_v - FEM_module::ElementTriangular<double>::C_V_AMB)/
		(-2*FEM_module::ElementTriangular<double>::SIGMA_VR*radius/
		 FEM_module::ElementBoundary<double>::RHO_V - pow(radius, 2));
	std::string filename;
	double r, z, val_u, val_v;
	int count = 1;
	gsl_vector* coeff_diff;
	for (size_t node : nodes_v){
		filename = "../Input/circle_" + std::to_string(count) + ".msh";
		NUM_NODES_G = node;
		FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
		FEM_module::ElementTriangular<double>::NUM_ELM = 0;
		FEM_module::ImporterMsh<double> mesh_importer(filename);
		mesh_importer.process_file();
		FEM_module::ConcentrationModel<double> model(mesh_importer, 
				interior_point);
		if (count == 7){
			coeff_diff = gsl_vector_alloc(2*model.number_nodes());
		}
		model.generate_f_vector();
		model.generate_stiffness_matrix();
		model.add_analytical_resp(c2_u, c2_v);
		model.solve_linear_model();
		rel_err_u = 0;
		rel_err_v = 0;
		norm_u = 0;
		norm_v = 0;
		for (size_t idx = 0; idx < NUM_NODES_G; idx++){
			r = model.coords()[idx][0];
			z = model.coords()[idx][1];
			val_u = c0_u + c2_u*(pow(r, 2) + pow(z - 0.07, 2));
			val_v = c0_v + c2_v*(pow(r, 2) + pow(z - 0.07, 2));
			if (count == 7){
				gsl_vector_set(coeff_diff, idx, gsl_vector_get(
							model.coefficients(), idx) - val_u);
				gsl_vector_set(coeff_diff, idx + NUM_NODES_G, gsl_vector_get(
							model.coefficients(), idx + NUM_NODES_G) - val_v);
			}
			rel_err_u += pow(gsl_vector_get(model.coefficients(), idx) - val_u, 
					2);
			rel_err_v += pow(gsl_vector_get(model.coefficients(), idx + 
						NUM_NODES_G) - val_v, 2);
			norm_u +=  pow(val_u, 2);
			norm_v +=  pow(val_v, 2);
		}
		myfile<<FEM_module::ElementTriangular<double>::NUM_ELM<<" "<<rel_err_u/norm_u<<
			" "<<rel_err_v/norm_v<<std::endl;
		if (count == 7){
			model.write_elements_to_file("../output/elements");
			model.write_boundaries_to_file("../output/boundaries");
			model.write_coordinates_to_file("../output/coords");
			FEM_module::write_vector_to_file(coeff_diff, "../output/initial_coeff");
			gsl_vector_free(coeff_diff);
		}
		count += 1;
	}
	myfile.close();	
	return EXIT_SUCCESS;
	
}

int test_no_o2(){
	FEM_module::ElementBoundary<double>::C_U_AMB = 0.0;
	FEM_module::ElementTriangular<double>::C_U_AMB = 0.0;
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.solve_nonlinear_model();
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}

int test_jacobian(){
	std::vector<double> interior_point{0.01, 0.06};
	gsl_vector* x;
	gsl_vector* f_val;
	gsl_permutation* permut;
	gsl_matrix* jac;
	gsl_matrix* jac_fd;
	gsl_matrix* stiff;
	gsl_matrix* inverse;
	int signum = 0;
	double norm_jac = 0;
	double norm_delta = 0;
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	x = gsl_vector_alloc(2*model.number_nodes());
	f_val = gsl_vector_alloc(2*model.number_nodes());
	permut = gsl_permutation_calloc(2*model.number_nodes());
	for (size_t idx_0 = 0; idx_0 < (size_t)model.number_nodes(); idx_0 ++){
		gsl_vector_set(x, idx_0, 101300*0.208/(8.314*298.15));
		gsl_vector_set(x, idx_0 + model.number_nodes(), 
				101300*0.0004/(8.314*298.15));
	}
	jac = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());
	jac_fd = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());
	stiff = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());
	inverse = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());

	gsl_spmatrix_sp2d(stiff, model.stiffness_matrix());
	FEM_module::NonLinearSystemFunctor<double> 
		nls_functor(model);
	FEM_module::JacobianFunctor<double> 
		jac_functor(model);
	FEM_module::FiniteDifferenceFunctor<double> 
		fd_functor(model, nls_functor, 1e-8);
	for (int iter = 0; iter < 10; iter++){
		jac_functor(x, NULL, jac);
		fd_functor(x, NULL, jac_fd);
		gsl_matrix_sub(jac, stiff);
		gsl_matrix_sub(jac_fd, stiff);
		for (size_t idx_1 = 0; idx_1 < jac->size1; idx_1++){
			for (size_t idx_2 = 0; idx_2 < jac->size1; idx_2++){
				norm_jac += std::pow(gsl_matrix_get(jac, idx_1, idx_2), 2);
				norm_delta += std::pow(gsl_matrix_get(jac, idx_1, idx_2) -
						gsl_matrix_get(jac_fd, idx_1, idx_2), 2);
			}
		}
		std::cout<<"Norm Jacobian: "<<std::sqrt(norm_jac)<<std::endl<<
			"Norm Difference: "<<std::sqrt(norm_delta)<<std::endl;
		nls_functor(x, NULL, f_val);
		gsl_linalg_LU_decomp(jac, permut, &signum);
		gsl_linalg_LU_invert(jac, permut, inverse);
		gsl_blas_dgemv(CblasNoTrans, -1.0, inverse, f_val, 1.0, x);
	}
	gsl_vector_free(x);
	gsl_matrix_free(jac);
	gsl_matrix_free(jac_fd);
	return EXIT_SUCCESS;
}

int test_jacobian_convergence(){
	std::ofstream myfile;
	myfile.open("fd_convergence", std::ios::out);
	std::vector<double> interior_point{0.01, 0.06};
	gsl_vector* x;
	gsl_matrix* jac;
	gsl_matrix* jac_fd;
	gsl_matrix* stiff;
	double norm_jac = 0;
	double norm_delta = 0;
	std::vector<double> eps_vect = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	x = gsl_vector_alloc(2*model.number_nodes());
	for (size_t idx_0 = 0; idx_0 < (size_t)model.number_nodes(); idx_0 ++){
		gsl_vector_set(x, idx_0, 101300*0.208/(8.314*298.15));
		gsl_vector_set(x, idx_0 + model.number_nodes(), 
				101300*0.0004/(8.314*298.15));
	}
	jac = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());
	jac_fd = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());
	stiff = gsl_matrix_alloc(2*model.number_nodes(), 2*model.number_nodes());

	gsl_spmatrix_sp2d(stiff, model.stiffness_matrix());
	FEM_module::NonLinearSystemFunctor<double> 
		nls_functor(model);
	FEM_module::JacobianFunctor<double> 
		jac_functor(model);
	FEM_module::FiniteDifferenceFunctor<double> 
		fd_functor(model, nls_functor, eps_vect[0]);
	for (size_t iter = 0; iter < eps_vect.size(); iter++){
		fd_functor.change_epsilon(eps_vect[iter]);
		jac_functor(x, NULL, jac);
		fd_functor(x, NULL, jac_fd);
		gsl_matrix_sub(jac, stiff);
		gsl_matrix_sub(jac_fd, stiff);
		for (size_t idx_1 = 0; idx_1 < jac->size1; idx_1++){
			for (size_t idx_2 = 1; idx_2 < jac->size1; idx_2++){
				norm_jac += std::pow(gsl_matrix_get(jac, idx_1, idx_2), 2);
				norm_delta += std::pow(gsl_matrix_get(jac, idx_1, idx_2) -
						gsl_matrix_get(jac_fd, idx_1, idx_2), 2);
			}
		}
		myfile<<eps_vect[iter]<<" "<<norm_delta/norm_jac<<std::endl;
	}
	gsl_vector_free(x);
	gsl_matrix_free(jac);
	gsl_matrix_free(jac_fd);
	myfile.close();
	return EXIT_SUCCESS;
}

int configure_case(std::ifstream& input_stream, std::string& line){
	double help_val;
	size_t help_int;
	std::vector<std::string> line_data;
	std::stringstream str_to_num;

	std::getline(input_stream, line);
	line_data = FEM_module::split(line, ' ');
	str_to_num << line_data[1];
	str_to_num >> help_val;
	str_to_num.clear();
	TEMP = 273.15 + help_val;

	std::getline(input_stream, line);
	line_data = FEM_module::split(line, ' ');
	str_to_num << line_data[1];
	str_to_num >> help_int;
	str_to_num.clear();
	NUM_NODES_G = help_int;

	std::getline(input_stream, line);
	line_data = FEM_module::split(line, ' ');
	str_to_num << line_data[1];
	str_to_num >> help_val;
	str_to_num.clear();
	CON_O2 = help_val/100.0;

	std::getline(input_stream, line);
	line_data = FEM_module::split(line, ' ');
	str_to_num << line_data[1];
	str_to_num >> help_val;
	str_to_num.clear();
	CON_CO2 = help_val/100.0;

	// Element values configuation
	FEM_module::ElementTriangular<double>::MAX_FERM_CO2 = 1.61e-4*exp(
			(56700/8.314)*(1/293.15 - 1/TEMP));
	FEM_module::ElementTriangular<double>::V_MU = 2.39e-4*exp(
			(80200/8.314)*(1.0/293.15 - 1.0/TEMP));
	FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
	
	FEM_module::ElementTriangular<double>::C_U_AMB = 101300.0*CON_O2/
		(8.314*TEMP);
	FEM_module::ElementTriangular<double>::C_V_AMB = 101300.0*CON_CO2/
		(8.314*TEMP);
	
	FEM_module::ElementBoundary<double>::C_U_AMB = 101300.0*CON_O2/
		(8.314*TEMP);
	FEM_module::ElementBoundary<double>::C_V_AMB = 101300.0*CON_CO2/
		(8.314*TEMP);
	FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
	return EXIT_SUCCESS;	
}

int program(std::string input_path){
	std::string mesh_file;
	std::string line;
	std::ifstream input_stream;
	input_stream.open(input_path);

	std::getline(input_stream, line);
	mesh_file = line;
	configure_case(input_stream, line);
	
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterText<double> mesh_importer(mesh_file);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.generate_stiffness_matrix();
	model.generate_f_vector();
	model.set_coefficients_to_amb(FEM_module::ElementTriangular<double>::C_U_AMB, 
				FEM_module::ElementTriangular<double>::C_V_AMB);
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/initial_coeff");
	model.generate_initial_codition();
	model.solve_nonlinear_model();
	FEM_module::write_vector_to_file(model.coefficients(),
		   	"../output/final_coeff");
	return EXIT_SUCCESS;
}

int correct_option(std::string message, int min_lim, int max_lim){
	int res;
	while (true){
		std::cout<<message;
		std::cin>>res;
		if ((res > min_lim - 1) && (res < max_lim  + 1)){
			return res;
		}
		std::cout<<"Incorrect option, please press enter to continue"<<
			std::endl<<std::endl;
	}
}

int main(int argc, char** argv){
	std::string prompt_to_run = "Select action\n" 
							"0 - Run case\n" 
							"1 - Perform test\n";
//	std::string prompt2 = "Select mesh type:\n" 
//							"0 - msh\n" 
//							"1 - text\n";
	std::string prompt_file = "Select Input File: ";
	std::string prompt_test = "Select Test Case\n"
							"1 - msh Importer\n"
							"2 - MATLAB mesh Importer\n"
							"3 - Pure Diffusion\n"
							"4 - Arbitrary Constant Respiration\n"
							"5 - Linear Model (Initial Condition)\n"
							"6 - Test Timing Linear Model\n"
							"7 - Sparse Non Linear Solver\n"
							"8 - Non Linear Solver\n"
							"9 - Stepped Method\n"
							"10 - Linear Respiration\n"
							"11 - Mesh Convergence to Analytical Solution\n"
							"12 - No O2 model\n"
							"13 - Jacobian Values\n"
							"14 - Jacobian Convergence\n";
	int cond1;
	int cond2;
	std::string input_path;

	cond1 = correct_option(prompt_to_run, 0, 1);
	if (cond1 == 0){
//		cond2 = correct_option(prompt2, input_1_2);
		std::cout<<prompt_file<<std::endl;
		std::cin>>input_path;
		program(input_path);
	}
	else{
		cond2 = correct_option(prompt_test, 1, 14);
		switch (cond2){
			case 1:
				test_importer();
				break;
			case 2:
				test_new_importer();
				break;
			case 3:
				test_pure_diffusion();
				break;
			case 4:
				test_constant_resp();
				break;
			case 5:
				test_initial_cond();
				break;
			case 6:
				test_timing_linear_model();
				break;
			case 7:
				test_concentration_model_sparse_nonlinear_solver();
				break;
			case 8:
				test_concentration_model_nonlinear();
				break;
			case 9:
				test_concentration_model_stepped_nonlinear();
				break;
			case 10:
				test_linear_resp();
				break;
			case 11:
				test_linear_solver_analytical();
				break;
			case 12:
				test_no_o2();
				break;
			case 13:
				test_jacobian();
				break;
			case 14:
				test_jacobian_convergence();
				break;
		}
	}
	return EXIT_SUCCESS;
}
