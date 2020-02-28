//Tests
#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#endif

#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include "element.hpp"
#endif
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t FEM_module::ElementTriangular<P, I>::SIGMA_u_r = 1;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t FEM_module::ElementTriangular<P, I>::SIGMA_u_z = 1;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t FEM_module::ElementTriangular<P, I>::SIGMA_v_r = 1;
template <typename P, typename I>
typename FEM_module::Element<P, I>::precision_t FEM_module::ElementTriangular<P, I>::SIGMA_v_z = 1;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t FEM_module::ElementTriangular<P, I>::NUM_ELM = 1;
template <typename P, typename I>
typename FEM_module::Element<P, I>::node_t FEM_module::ElementTriangular<P, I>::NUM_NODES = 1;
//FEM_module::ElementTriangular<double, long>::SIGMA_u_z = 1;
//FEM_module::ElementTriangular<double, long>::SIGMA_v_r = 1;
//FEM_module::ElementTriangular<double, long>::SIGMA_v_z = 1;
//FEM_module::ElementTriangular<double, long>::NUM_NODES = 5;
//FEM_module::ElementTriangular<double, long>::NUM_ELM = 0;
//FEM_module::ElementTriangular<double, long>::INITIALIZED = true;
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
int main(){
	//Create ConcentrationModel
	//Initialize Concentration model wjit the elements
	//ConcentrationModel.create_non_linear()
	//ConentrationModel.solve()
	test2();
	return EXIT_SUCCESS;
}
