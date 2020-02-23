#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>

#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#endif


namespace FEM_module{

template <typename P, typename I>
class ElementTriangular{
	public:
		typedef P precision_t;
		typedef I node_t;
	protected:
		std::vector<node_t> nodes_;
	public:
		precision_t static SIGMA_u_r = 0;
		precision_t static SIGMA_u_z = 0;
		precision_t static SIGMA_v_r = 0;
		precision_t static SIGMA_v_z = 0;
		node_t static NUM_ELM = 0;
		/*
		 *
		 *
		 *
		Ask if the derived constructor is going to update the static member 
		number_elements
		 *
		 *
		 *
		*/
		node_t static NUM_NODES = 0;

		ElementTriangular(const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& nodes)
		: nodes_(std::vector<node_t>())
		{
			bool current_rotation = get_rotation(coordinates, nodes);
			if (current_rotation){
				for (int i = 0; i < 3; i++){
					nodes_.push_back(nodes[i]);
				}
			}
			else{
				for (int i = 0; i < 3; i++){
					nodes_.push_back(nodes[2 - i]);
				}
			}
			NUM_ELM++;
		}

		int get_global_coords(const std::vector<precision_t>& coordinates, 
				std::vector<precision_t>& output_coords){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 2; j++){
					output_coords[i][j] = coordinates[i][j];
				}
			}
			return EXIT_SUCCESS;
		}
		
		// Access
		const std::vector<node_t>& nodes(){
			return nodes_;
		}

	private:
		bool get_rotation(const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& nodes){
		/* Function use to get the orientation of the nodes, in order to order
		them always counter-clockwise.
		Output:
			True: Counter-Clockwise order.
			False: Clockwise order.
		*/
		precision_t slope_1 = 
			(coordinates[nodes[1] - 1][2] - coordinates[nodes[0] - 1][2])/
			(coordinates[nodes[1] - 1][1] - coordinates[nodes[0] - 1][1]);
		precision_t slope_2 = 
			(coordinates[nodes[2] - 1][2] - coordinates[nodes[1] - 1][2])/
			(coordinates[nodes[2] - 1][1] - coordinates[nodes[1] - 1][1]);
		return (slope_1 < slope_2);
		}

	public:
		precision_t jacobian(
				const std::vector<std::vector<precision_t>>& coordinates){
			precision_t r1 = coordinates[nodes_[1]][0];
			precision_t r2 = coordinates[nodes_[2]][0];
			precision_t r3 = coordinates[nodes_[3]][0];
			precision_t z1 = coordinates[nodes_[1]][1];
			precision_t z2 = coordinates[nodes_[2]][1];
			precision_t z3 = coordinates[nodes_[3]][1];
			return (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1);
		}

		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<std::vector<precision_t>>& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			precision_t r1 = coordinates[nodes_[1]][0];
			precision_t r2 = coordinates[nodes_[2]][0];
			precision_t r3 = coordinates[nodes_[3]][0];
			precision_t z1 = coordinates[nodes_[1]][1];
			precision_t z2 = coordinates[nodes_[2]][1];
			precision_t z3 = coordinates[nodes_[3]][1];
			/*
			 *
			 *
			 Ask about casting
			 *
			 *
			 */
			precision_t k = (r1 + (r2 - r1 + r3 - r1)/
					((precision_t) 3))/((precision_t) 2);
			k = k*jacobian();
			
			/*
			Here we update the matrix. We still havent choose a data structure 
			for the stiffness matrix, but it will probably accept the matrix
			interface
			*/

			global_stiffness[nodes_[1]-1][nodes_[1]-1] += (
					SIGMA_u_r + SIGMA_u_z)*k;
			global_stiffness[nodes_[1]-1][nodes_[2]-1] -= SIGMA_u_r*k;
			global_stiffness[nodes_[1]-1][nodes_[3]-1] -= SIGMA_u_z*k;
			global_stiffness[nodes_[2]-1][nodes_[1]-1] -= SIGMA_u_r*k;
			global_stiffness[nodes_[2]-1][nodes_[2]-1] += SIGMA_u_r*k;
			//global_stiffness[nodes_[2]-1][nodes_[3]-1] += 0;
			global_stiffness[nodes_[3]-1][nodes_[1]-1] -= SIGMA_u_z*k;
			//global_stiffness[nodes_[3]-1][nodes_[2]-1] += 0;
			global_stiffness[nodes_[3]-1][nodes_[3]-1] += SIGMA_u_z*k;

			global_stiffness[nodes_[1]-1+NUM_NODES][nodes_[1]-1+NUM_NODES] += 
				(SIGMA_v_r + SIGMA_v_z)*k;
			global_stiffness[nodes_[1]-1+NUM_NODES][nodes_[2]-1+NUM_NODES] -= 
				SIGMA_v_r*k;
			global_stiffness[nodes_[1]-1+NUM_NODES][nodes_[3]-1+NUM_NODES] -=
			   	SIGMA_v_z*k;
			global_stiffness[nodes_[2]-1+NUM_NODES][nodes_[1]-1+NUM_NODES] -= 
				SIGMA_v_r*k;
			global_stiffness[nodes_[2]-1+NUM_NODES][nodes_[2]-1+NUM_NODES] +=
			   	SIGMA_v_r*k;
			//global_stiffness[nodes_[2]-1+NUM_NODES][nodes_[3]-1+NUM_NODES] += 
			//	0;
			global_stiffness[nodes_[3]-1+NUM_NODES][nodes_[1]-1+NUM_NODES] -= 
				SIGMA_v_z*k;
			//global_stiffness[nodes_[3]-1+NUM_NODES][nodes_[2]-1+NUM_NODES] += 
			//	0;
			global_stiffness[nodes_[3]-1+NUM_NODES][nodes_[3]-1+NUM_NODES] += 
				SIGMA_v_z*k;
			return EXIT_SUCCESS;
		}

};

template <typename P, typename I>
class ElementBoundary : public ElementTriangular<P, I>{
	public:
		using typename ElementTriangular<P, I>::precision_t;
		using typename ElementTriangular<P, I>::node_t;
	private:
		std::vector<node_t> boundary_nodes;
	public:
		ElementBoundary(
				const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& nodes)
		: ElementTriangular<precision_t, node_t>(coordinates, nodes)
		{}

		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<std::vector<precision_t>>& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			return EXIT_SUCCESS;
		}

		int update_vector_f(
				const std::vector<std::vector<precision_t>>& coordinates, 
				std::vector<precision_t>& vector_f){
				
			return EXIT_SUCCESS;
		}
};
}

int test2(){
	FEM_module::ImporterMsh<double, long> mesh_importer("../Input/pear.msh");
	mesh_importer.process_file();
	const std::vector<std::vector<double>> coords = 
		mesh_importer.node_matrix();
	std::vector<long> nodes_from_element;
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
