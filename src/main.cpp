//Tests
#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#include <cmath>
#include <string>
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
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_UR = 2.8e-10;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_UZ = 1.10e-9;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_VR = 2.32e-9;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_VZ = 6.97e-9;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::RESP_Q = 0.97;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::MAX_FERM_CO2 = 1.61e-4*exp(
			(56700/8.314)*(1/293.15 - 1/TEMP));
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::K_MFU = 0.1149;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::K_MU = 0.4103;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::K_MV = 27.2438;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::V_MU = 2.39e-4*exp(
			(80200/8.314)*(1/293.15 - 1/TEMP));
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementTriangular<P, I>::NUM_ELM = 0;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementTriangular<P, I>::NUM_NODES = NUM_NODES_G;

template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::C_U_AMB = 101300*CON_O2/(8.314*TEMP);
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::C_V_AMB = 101300*CON_CO2/(8.314*TEMP);

// boundary values configuration
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::C_U_AMB = 101300*CON_O2/(8.314*TEMP);
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::C_V_AMB = 101300*CON_CO2/(8.314*TEMP);
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::RHO_U = 7e-7;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::RHO_V = 7.5e-7;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementBoundary<P, I>::NUM_NODES = NUM_NODES_G;

int test2(){
	FEM_module::ImporterMsh<double, long> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	const std::vector<std::vector<double>> coords = 
		mesh_importer.node_matrix();
	std::vector<long> nodes_from_element;
	std::cout<<"Node Num - x - y"<<std::endl;
	long count = 1;
    for (auto elem : coords) {
		std::cout<<"Node "<<count<<" - "<<elem[0]<<" - "<<elem[1]<<std::endl;
		count++;
    }
	std::cout<<std::endl;
	for (auto elem : mesh_importer.element_matrix()){
		FEM_module::ElementTriangular<double, long> temp_element(coords, elem);
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
	std::vector<int> nodes;
	gsl_vector* coeff = gsl_vector_alloc(6);
	gsl_vector* result = gsl_vector_alloc(6);
	gsl_vector_set_all(coeff, 1.0);
	for (int i = 0; i < 3; i++){
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
	FEM_module::ElementTriangular<double, int>::NUM_NODES = 3;
	FEM_module::ElementTriangular<double, int> test_element(coords, nodes);
	test_element.r_u(coeff, 0.5, 0.1);
	test_element.r_v(coeff, 0.5, 0.1);
	test_element.integrate_non_linear_term(coeff, coords, 15, result);

	gsl_vector_free(coeff);
	gsl_vector_free(result);
	return EXIT_SUCCESS;
};

int test_linear_model(){
	std::ifstream file_stream;
	std::stringstream str_to_num;
	std::string line;
	std::vector<std::string> line_data;
	gsl_vector* helper;
	std::vector<double> interior_point{0.01, 0.06};
	size_t count = 0;
	double help_num = 0;
	FEM_module::ImporterMsh<double, int> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	helper = gsl_vector_alloc((size_t)(2*model.number_nodes()));
	
	file_stream.open("../output/f_calc");
	std::getline(file_stream, line);
	line_data = FEM_module::split(line, ' ');
	for (auto val : line_data){
		std::cout<<val<<std::endl;
		str_to_num << val;
		str_to_num >> help_num;
		gsl_vector_set(helper, count, help_num);
		str_to_num.clear();
		count++;
	}
	count = 0;
	file_stream.close();
	model.set_f_vector(helper);
	FEM_module::write_vector_to_file(helper, "f_cal_loaded");
	file_stream.open("../output/calculated_coeff");
	std::getline(file_stream, line);
	line_data = FEM_module::split(line, ' ');
	for (auto val : line_data){
		str_to_num << val;
		str_to_num >> help_num;
		gsl_vector_set(helper, count, help_num);
		str_to_num.clear();
		count++;
	}
	file_stream.close();
	
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
	FEM_module::ImporterMsh<double, int> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	model.solve_linear_model();
	return EXIT_SUCCESS;
}

int test_concentration_model_2(){
	std::vector<double> interior_point{0.01, 0.06};
	FEM_module::ImporterMsh<double, int> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
			interior_point);
	model.write_elements_to_file("../output/elements");
	model.write_boundaries_to_file("../output/boundaries");
	model.write_coordinates_to_file("../output/coords");
	model.solve_nonlinear_model();
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
	FEM_module::ImporterMsh<double, int> mesh_importer(FILEPATH);
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
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
	FEM_module::NonLinearSystemFunctor<double, int> 
		nls_functor(model, 3);
	FEM_module::JacobianFunctor<double, int> 
		jac_functor(model);
	FEM_module::FiniteDifferenceFunctor<double, int> 
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

int main(){
	test_linear_model();
	return EXIT_SUCCESS;
}
