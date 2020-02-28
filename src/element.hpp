#include <cmath>
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
class Element{
	public:
		typedef P precision_t;
		typedef I node_t;
	protected:
		std::vector<node_t> nodes_;
	public:
		Element()
		: nodes_(std::vector<node_t>())
		{}
		virtual int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<std::vector<precision_t>>& global_stiffness) = 0;
};

template <typename P, typename I>
class ElementTriangular : public Element<P, I>{
	public:
		using typename Element<P, I>::precision_t;
		using typename Element<P, I>::node_t;
	public:
		// static class members
		static precision_t SIGMA_u_r;
		static precision_t SIGMA_u_z;
		static precision_t SIGMA_v_r;
		static precision_t SIGMA_v_z;
		static node_t NUM_ELM;
		static node_t NUM_NODES;
		//static bool INITIALIZED;
		
		// Constructord
		ElementTriangular(const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& p_nodes)
		{
			bool current_rotation = get_orientation(coordinates, p_nodes);
			if (current_rotation){
				for (int i = 0; i < 3; i++){
					this->nodes_.push_back(p_nodes[i]);
				}
			}
			else{
				for (int i = 0; i < 3; i++){
					this->nodes_.push_back(p_nodes[2 - i]);
				}
			}
			ElementTriangular<precision_t, node_t>::NUM_ELM++;
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
			return this->nodes_;
		}

	private:
		bool get_orientation(
				const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& p_nodes){
		/* 
		 * Function use to get the orientation of the nodes, in order to order 
		 * them always counter-clockwise.
		 * Output:
		 *	True: Counter-Clockwise order.
		 *	False: Clockwise order.
		*/
		precision_t slope_1 = 
			(coordinates[p_nodes[1]][2] - coordinates[p_nodes[0]][2])/
			(coordinates[p_nodes[1]][1] - coordinates[p_nodes[0]][1]);
		precision_t slope_2 = 
			(coordinates[p_nodes[2]][2] - coordinates[p_nodes[1]][2])/
			(coordinates[p_nodes[2]][1] - coordinates[p_nodes[1]][1]);
		return (slope_1 < slope_2);
		}
		
		precision_t phi_1(precision_t epsilon, precision_t eta){
			return 1 - epsilon - eta;
		}
		
		precision_t phi_2(precision_t epsilon, precision_t eta){
			return epsilon;
		}
		
		precision_t phi_3(precision_t epsilon, precision_t eta){
			return eta;
		}

	public:
		inline precision_t jacobian(
				const std::vector<std::vector<precision_t>>& coordinates){
			/*
			 * returns the jacobian of the element, when transformed to a
			 * rectangle with minnor sides euqal to 1. The jacobian is related
			 * to the area of the element
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			precision_t z1 = coordinates[this->nodes_[0]][1];
			precision_t z2 = coordinates[this->nodes_[1]][1];
			precision_t z3 = coordinates[this->nodes_[2]][1];
			return (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1);
		}

		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<std::vector<precision_t>>& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			precision_t z1 = coordinates[this->nodes_[0]][1];
			precision_t z2 = coordinates[this->nodes_[1]][1];
			precision_t z3 = coordinates[this->nodes_[2]][1];
			/*
			 *
			 *
			 Ask about casting
			 *
			 *
			 */
			precision_t k = (r1 + (r2 - r1 + r3 - r1)/
					((precision_t) 3))/((precision_t) 2);
			k = k*jacobian(coordinates);
			
			/*
			 * Here we update the matrix. We still havent choose a data 
			 * structure for the stiffness matrix, but it will probably accept 
			 * the matrix interface
			*/

			global_stiffness[this->nodes_[0]][this->nodes_[0]] += (
					SIGMA_u_r + SIGMA_u_z)*k;
			global_stiffness[this->nodes_[0]][this->nodes_[1]] -= SIGMA_u_r*k;
			global_stiffness[this->nodes_[0]][this->nodes_[2]] -= SIGMA_u_z*k;
			global_stiffness[this->nodes_[1]][this->nodes_[0]] -= SIGMA_u_r*k;
			global_stiffness[this->nodes_[1]][this->nodes_[1]] += SIGMA_u_r*k;
			//global_stiffness[nodes_[1]][nodes_[2]] += 0;
			global_stiffness[this->nodes_[2]][this->nodes_[0]-1] -= SIGMA_u_z*k;
			//global_stiffness[nodes_[2]][nodes_[1]] += 0;
			global_stiffness[this->nodes_[2]][this->nodes_[2]] += SIGMA_u_z*k;

			global_stiffness[this->nodes_[0]+NUM_NODES][this->nodes_[0]+NUM_NODES] += 
				(SIGMA_v_r + SIGMA_v_z)*k;
			global_stiffness[this->nodes_[0]+NUM_NODES][this->nodes_[1]+NUM_NODES] -= 
				SIGMA_v_r*k;
			global_stiffness[this->nodes_[0]+NUM_NODES][this->nodes_[2]+NUM_NODES] -=
			   	SIGMA_v_z*k;
			global_stiffness[this->nodes_[1]+NUM_NODES][this->nodes_[0]+NUM_NODES] -= 
				SIGMA_v_r*k;
			global_stiffness[this->nodes_[1]+NUM_NODES][this->nodes_[1]+NUM_NODES] +=
			   	SIGMA_v_r*k;
			//global_stiffness[nodes_[1]+NUM_NODES][nodes_[2]+NUM_NODES] += 
			//	0;
			global_stiffness[this->nodes_[2]+NUM_NODES][this->nodes_[0]+NUM_NODES] -= 
				SIGMA_v_z*k;
			//global_stiffness[nodes_[2]+NUM_NODES][nodes_[1]+NUM_NODES] += 
			//	0;
			global_stiffness[this->nodes_[2]+NUM_NODES][this->nodes_[2]+NUM_NODES] += 
				SIGMA_v_z*k;
			return EXIT_SUCCESS;
		}
};

template <typename P, typename I>
class ElementBoundary : public Element<P, I>{
	public:
		using typename Element<P, I>::precision_t;
		using typename Element<P, I>::node_t;
	public:
		static precision_t C_AMB;
		ElementBoundary(
				const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& p_nodes, 
				const std::vector<precision_t> interior_point)
		{
			if (get_orientation(coordinates, p_nodes, interior_point)){
				for (int i = 0; i < 2; i++){
					this->nodes_.push_back(p_nodes[i]);
				}
			}
			else{
				for (int i = 0; i < 2; i++){
					this->nodes_.push_back(p_nodes[1 - i]);
				}
			}
		
		}

	private:
		inline precision_t length(
				std::vector<std::vector<precision_t>>& coordinates){
			return sqrt(pow((coordinates[this->nodes_[0]][0] - 
						coordinates[this->nodes_[1]][0]), 2) +
				pow((coordinates[this->nodes_[0]][1] - 
							coordinates[this->nodes_[1]][1]), 2));
		}
		
		bool get_orientation(
				std::vector<std::vector<precision_t>>& coordinates,
				std::vector<precision_t>& p_nodes,
				std::vector<precision_t>& interior_point){
			/*
			 * Computes the orientation of the segment that defines the 
			 * boundary element, when given an interior point of the closed
			 * surface intetior_point, returns true if the element is
			 * counterclockwise oriented, and false if its clockwise
			*/
			precision_t ang1 = atan2(
					coordinates[p_nodes[0]][1] - interior_point[1], 
					coordinates[p_nodes[0]][0] - interior_point[0]);
			precision_t ang2 = atan2(
					coordinates[p_nodes[1]][1] - interior_point[1], 
					coordinates[p_nodes[1]][0] - interior_point[0]);
			return ang2 > ang1;
		}
		
	public:
		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<std::vector<precision_t>>& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t z1 = coordinates[this->nodes_[0]][1];
			precision_t z2 = coordinates[this->nodes_[1]][1];
			
			global_stiffness[this->nodes_[0]][this->nodes_[0]] += 
				length(coordinates)*(r2/4  + r1/12);
			global_stiffness[this->nodes_[0]][this->nodes_[1]] += 
				length(coordinates)*(r1+r2)/12;
			global_stiffness[this->nodes_[1]][this->nodes_[0]] += 
				length(coordinates)*(r1+r2)/12;
			global_stiffness[this->nodes_[1]][this->nodes_[1]] +=
				length(coordinates)*(r1/4  + r2/12);
			return EXIT_SUCCESS;
		}

		int update_vector_f(
				const std::vector<std::vector<precision_t>>& coordinates, 
				std::vector<precision_t>& vector_f){
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			vector_f[this->nodes_[0]] += length(coordinates)*
				C_AMB*(r1/6 + r2/3);
			vector_f[this->nodes_[1]] += length(coordinates)*
				C_AMB*(r2/6 + r1/3);
			return EXIT_SUCCESS;
		}
};
}

