//Tests
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

// Element values configuation
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_UR = 2.8e-10;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementTriangular<P, I>::SIGMA_UZ = 1.1e-9;
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
	FEM_module::ElementTriangular<P, I>::MAX_FERM_CO2 = 1.61e-4;
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
	FEM_module::ElementTriangular<P, I>::V_MU = 2.39e-4;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementTriangular<P, I>::NUM_ELM = 0;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementTriangular<P, I>::NUM_NODES = 42;

// boundary values configuration
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::C_U_AMB = 101300*0.208/(8.32*293);
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::C_V_AMB = 0;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::RHO_U = 7e-7;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t 
	FEM_module::ElementBoundary<P, I>::RHO_V = 7.5e-7;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t 
	FEM_module::ElementBoundary<P, I>::NUM_NODES = 42;

int test2(){
	FEM_module::ImporterMsh<double, long> mesh_importer("../Input/pear.msh");
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
	gsl_vector_set_all(coeff, 1.1);
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
	test_element.integrate_non_linear_term(coeff, coords, 5, result);

	gsl_vector_free(coeff);
	gsl_vector_free(result);
	return EXIT_SUCCESS;
};

int test_concentration_model_1(){
	std::vector<double> interior_point{60, 5};
	FEM_module::ImporterMsh<double, int> mesh_importer("../Input/pear.msh");
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	model.generate_f_vector();
	model.solve_linear_model();
	return EXIT_SUCCESS;
}

int test_concentration_model_2(){
	std::vector<double> interior_point{60, 5};
	FEM_module::ImporterMsh<double, int> mesh_importer("../Input/pear.msh");
	mesh_importer.process_file();
	FEM_module::ConcentrationModel<double, int> model(mesh_importer, 
			interior_point);
	model.generate_stiffness_matrix();
	model.generate_f_vector();
	model.solve_nonlinear_model();
	std::cout<<model<<std::endl;
	return EXIT_SUCCESS;
}

int main(){
	//Create ConcentrationModel
	//Initialize Concentration model wjit the elements
	//ConcentrationModel.create_non_linear()
	//ConentrationModel.solve()
	test_concentration_model_2();
	return EXIT_SUCCESS;
}
