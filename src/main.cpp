//Tests
#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>


#include <chrono>
#include <thread>
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
int NUM_NODES_G = 529;
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

int test2(){

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
int test_quadrature(){
	std::vector<std::vector<double>> coords;
	std::vector<size_t> nodes;
	gsl_vector* coeff = gsl_vector_alloc(6);
	gsl_vector* result = gsl_vector_alloc(6);
	gsl_vector_set_all(coeff, 1.0);
	for (size_t i = 0; i < 3; i++){
		nodes.push_back(i);
		if (i == 0){
			coords.push_back(std::vector<double>({0.34, 1.2}));
		}
		else if (i == 1){
			coords.push_back(std::vector<double>({1.70, 0.10}));
		}
		else{
			coords.push_back(std::vector<double>({3.66, 2.50}));
		}
	}
	FEM_module::ElementTriangular<double>::NUM_NODES = 3;
	FEM_module::ElementTriangular<double> test_element(coords, nodes);
	test_element.r_u(coeff, 0.5, 0.1);
	test_element.r_v(coeff, 0.5, 0.1);
	test_element.integrate_non_linear_term(coeff, coords, result);

	gsl_vector_free(coeff);
	gsl_vector_free(result);
	return EXIT_SUCCESS;
};

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

int test_linear_model(){
	std::ifstream file_stream;
	std::stringstream str_to_num;
	std::string line;
	std::vector<std::string> line_data;
	gsl_vector* helper;
	std::vector<double> interior_point{0.01, 0.06};
	size_t count = 0;
	double help_num = 0;
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	helper = gsl_vector_alloc((size_t)(2*model.number_nodes()));
	
	file_stream.open("../output/f_calc");
    while (std::getline(file_stream, line)) {
		if (line == "4"){}
		else {
			line_data = FEM_module::split(line, ' ');
			for (auto val : line_data){
				std::cout<<val<<std::endl;
				str_to_num << val;
				str_to_num >> help_num;
				gsl_vector_set(helper, count, help_num);
				str_to_num.clear();
				count++;
			}
		}
    }
	count = 0;
	file_stream.close();
	//model.set_f_vector(helper);

	FEM_module::write_vector_to_file(helper, "f_cal_loaded");
	file_stream.open("../output/calculated_coeff");
	std::getline(file_stream, line);
    while (std::getline(file_stream, line)) {
		if (line == "4"){}
		else {
			line_data = FEM_module::split(line, ' ');
			for (auto val : line_data){
				std::cout<<val<<std::endl;
				str_to_num << val;
				str_to_num >> help_num;
				gsl_vector_set(helper, count, help_num);
				str_to_num.clear();
				count++;
			}
		}
    }
	
	model.set_f_vector(helper);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_linear_model_LU();
	FEM_module::write_vector_to_file(model.coefficients(), 
			"../output/initial_coeff");
	gsl_vector_sub(helper, model.coefficients());
	std::cout<<FEM_module::vector_to_string(helper);
	return EXIT_SUCCESS;
}

int test_concentration_model_1(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	model.solve_linear_model();
	return EXIT_SUCCESS;
}

int test_concentration_model_2(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_nonlinear_model();
	return EXIT_SUCCESS;
}

int test_concentration_model_3(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_stepped_nonlinear_model(14);
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
		model.solve_linear_model_LU();
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
			rel_err_u += pow(gsl_vector_get(model.coefficients(), idx) - val_u, 2);
			rel_err_v += pow(gsl_vector_get(model.coefficients(), idx + NUM_NODES_G) - val_v, 2);
			norm_u +=  pow(val_u, 2);
			norm_v +=  pow(val_v, 2);
		}
		myfile<<NUM_NODES_G<<" "<<rel_err_u/norm_u<<" "<<rel_err_v/norm_v<<std::endl;
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

int test_linear_solver(){
	NUM_NODES_G = 779;
	FEM_module::ElementTriangular<double>::NUM_NODES = NUM_NODES_G;
	FEM_module::ElementBoundary<double>::NUM_NODES = NUM_NODES_G;
	FEM_module::ElementTriangular<double>::SIGMA_UZ = 
		FEM_module::ElementTriangular<double>::SIGMA_UR;
	FEM_module::ElementTriangular<double>::SIGMA_VZ = 
		FEM_module::ElementTriangular<double>::SIGMA_VR;
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double> mesh_importer("../Input/circle.msh");
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.generate_f_vector();
	model.generate_stiffness_matrix();
	double radius = 0.06;
	double r;
	double z;
	double val_u = 0;
	double val_v = 0;
	double c0_u = 2.0;
	double c2_u = (c0_u - FEM_module::ElementTriangular<double>::C_U_AMB)/
		(-2*FEM_module::ElementTriangular<double>::SIGMA_UR*radius/
		 FEM_module::ElementBoundary<double>::RHO_U - pow(radius, 2));
	double c0_v = 25.0;
	double c2_v = (c0_v - FEM_module::ElementTriangular<double>::C_V_AMB)/
		(-2*FEM_module::ElementTriangular<double>::SIGMA_VR*radius/
		 FEM_module::ElementBoundary<double>::RHO_V - pow(radius, 2));
	gsl_vector* wanted_coeff;
	wanted_coeff = gsl_vector_alloc(2*model.number_nodes());
	gsl_vector* wanted_f;
	wanted_f = gsl_vector_calloc(2*model.number_nodes());
	for (size_t idx = 0; idx < model.number_nodes(); idx++){
		r = model.coords()[idx][0];
		z = model.coords()[idx][1];
		val_u = c0_u + c2_u*(pow(r, 2) + pow(z - 0.07, 2));
		val_v = c0_v + c2_v*(pow(r, 2) + pow(z - 0.07, 2));
		gsl_vector_set(wanted_coeff, idx, val_u);
		gsl_vector_set(wanted_coeff, idx + model.number_nodes(), val_v);
	}
	gsl_spblas_dgemv(CblasNoTrans, 1.0, model.stiffness_matrix(), 
			wanted_coeff, 0.0, wanted_f);
	gsl_vector_add(wanted_f, model.f_vector());
	gsl_vector_scale(wanted_f, -1.0);;
	FEM_module::write_vector_to_file(wanted_f, 
					"../output/f_vector");
	gsl_vector_add(wanted_f, model.f_vector());
	gsl_vector_scale(wanted_f, -1.0);
	model.set_f_vector(wanted_f);
	model.solve_linear_model();
	FEM_module::write_vector_to_file(model.coefficients(), 
					"../output/initial_coeff");
	std::cout<<"C0 u: "<<c0_u<<" C0 v: "<<c0_v<<std::endl;
	std::cout<<"C2 u: "<<c2_u<<" C2 v: "<<c2_v<<std::endl;
	std::cout<<"C_r u: "<<c0_u + c2_u*radius*radius<<std::endl;
	std::cout<<"C_r v: "<<c0_v + c2_v*radius*radius<<std::endl;
	std::cout<<"C_V_AMB: "<<FEM_module::ElementTriangular<double>::C_V_AMB<<std::endl;
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
	gsl_vector* f_val;
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
	f_val = gsl_vector_alloc(2*model.number_nodes());
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
	for (int iter = 0; iter < eps_vect.size(); iter++){
		fd_functor.change_epsilon(eps_vect[iter]);
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
		myfile<<eps_vect[iter]<<" "<<norm_delta/norm_jac<<std::endl;
	}
	gsl_vector_free(x);
	gsl_matrix_free(jac);
	gsl_matrix_free(jac_fd);
	myfile.close();
	return EXIT_SUCCESS;
}

int main(){
	test_jacobian_convergence();
	return EXIT_SUCCESS;
}
