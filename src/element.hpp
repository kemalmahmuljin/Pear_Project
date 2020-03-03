#include <cmath>
#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_integration.h>

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
		typedef gsl_spmatrix stiff_mat_t;
		typedef gsl_vector global_vect_t;
	protected:
		std::vector<node_t> nodes_;
	public:
		Element()
		: nodes_(std::vector<node_t>())
		{}
		virtual int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				stiff_mat_t& global_stiffness) = 0;
};

template <typename P, typename I>
class ElementTriangular : public Element<P, I>{
	public:
		using typename Element<P, I>::precision_t;
		using typename Element<P, I>::node_t;
		using typename Element<P, I>::stiff_mat_t;
		using typename Element<P, I>::global_vect_t;
	public:
		// static class members
		static precision_t SIGMA_UR;
		static precision_t SIGMA_UZ;
		static precision_t SIGMA_VR;
		static precision_t SIGMA_VZ;
		static precision_t RESP_Q;
		static precision_t MAX_FERM_CO2;
		static precision_t K_MFU;
		static precision_t K_MU;
		static precision_t K_MV;
		static precision_t V_MU;
		static node_t NUM_ELM;
		static node_t NUM_NODES;
		//static bool INITIALIZED;
		
		// Constructor
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
			(coordinates[p_nodes[1]][1] - coordinates[p_nodes[0]][1])/
			(coordinates[p_nodes[1]][0] - coordinates[p_nodes[0]][0]);
		precision_t slope_2 = 
			(coordinates[p_nodes[2]][1] - coordinates[p_nodes[1]][1])/
			(coordinates[p_nodes[2]][0] - coordinates[p_nodes[1]][0]);
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

		precision_t r_u_num(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			precision_t p_1 = phi_1(epsilon, eta);
			precision_t p_2 = phi_2(epsilon, eta);
			precision_t p_3 = phi_3(epsilon, eta);

			precision_t res = 
				p_1*gsl_vector_get(coefficients, (size_t)this->nodes_[0]) +
				p_2*gsl_vector_get(coefficients, (size_t)this->nodes_[1]) +
				p_3*gsl_vector_get(coefficients, (size_t)this->nodes_[2]);
			res *= V_MU;
			return res;
		}

		precision_t r_u_den(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			precision_t p_1 = phi_1(epsilon, eta);
			precision_t p_2 = phi_2(epsilon, eta);
			precision_t p_3 = phi_3(epsilon, eta);
			precision_t c_u_1 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[0]);
			precision_t c_u_2 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[1]);
			precision_t c_u_3 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[2]);
			precision_t c_v_1 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[0] + NUM_NODES);
			precision_t c_v_2 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[1] + NUM_NODES);
			precision_t c_v_3 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[2] + NUM_NODES);
			precision_t res = 0;
			res += K_MU*(1 + (c_v_1/K_MV + c_u_1)*p_1 +
				(c_v_2/K_MV + c_u_2)*p_2 +
				(c_v_3/K_MV + c_u_3)*p_3);
			res += (c_u_1*c_v_1*pow(p_1, 2) + 
				c_u_2*c_v_2*pow(p_2, 2) +
				c_u_3*c_v_3*pow(p_3, 2))/K_MV;
			res += (c_u_1*p_1*(c_v_2*p_2 + c_v_3*p_3) +
				c_u_2*p_2*(c_v_1*p_1 + c_v_3*p_3) +
				c_v_3*p_3*(c_v_1*p_1 + c_v_2*p_2))/K_MV;
			return res;
		}

	public:
		precision_t r_u(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			return r_u_num(coefficients, epsilon, eta)/
				r_u_den(coefficients, epsilon, eta);
		}

		precision_t r_v(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			precision_t p_1 = phi_1(epsilon, eta);
			precision_t p_2 = phi_2(epsilon, eta);
			precision_t p_3 = phi_3(epsilon, eta);
			precision_t c_u_1 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[0]);
			precision_t c_u_2 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[1]);
			precision_t c_u_3 = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[2]);
			precision_t c_u = p_1*c_u_1 + p_2*c_u_2 + p_3*c_u_3; 
			return RESP_Q*r_u(coefficients, epsilon, eta) + 
				MAX_FERM_CO2/(1 + c_u/K_MFU);
		}

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

		precision_t integrand_u(precision_t u, precision_t v, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx){
			/*
			 * returns the integrand r*R_u(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a rectangular domain with boundaries [0,1]
			 * and [0,1], it saves the results for phi_j (j=1:3) inside 
			 * integrand_res.
			 * Note that integrand_v hast 3 components, it's not a gloval ector
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			if (node_idx == 1){
				return (r1 + (r2 - r1))*u*r_u(coefficients, u, (1-u)*v)*
					phi_1(u, (1-u)*v)*jacobian(coordinates);
			}else if (node_idx == 2){
				return (r1 + (r2 - r1))*u*r_u(coefficients, u, (1-u)*v)*
					phi_2(u, (1-u)*v)*jacobian(coordinates);
			}else if (node_idx == 3){
				return (r1 + (r2 - r1))*u*r_u(coefficients, u, (1-u)*v)*
					phi_2(u, (1-u)*v)*jacobian(coordinates);
			}
		}
		
		precision_t integrand_v(precision_t u, precision_t v, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx){
			/*
			 * returns the integrand r*R_u(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a rectangular domain with boundaries [0,1]
			 * and [0,1], it saves the results for phi_j (j=1:3) inside 
			 * integrand_res.
			 * Note that integrand_v hast 3 components, it's not a gloval ector
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			if (node_idx == 1){
				return (r1 + (r2 - r1))*u*r_v(coefficients, u, (1-u)*v)*
					phi_1(u, (1-u)*v)*jacobian(coordinates);
			}else if (node_idx == 2){
				return (r1 + (r2 - r1))*u*r_v(coefficients, u, (1-u)*v)*
					phi_2(u, (1-u)*v)*jacobian(coordinates);
			}else if (node_idx == 3){
				return (r1 + (r2 - r1))*u*r_v(coefficients, u, (1-u)*v)*
					phi_2(u, (1-u)*v)*jacobian(coordinates);
			}
		}

		int integrate_non_linear_term(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				size_t points, gsl_vector* result_vector){
			gsl_integration_glfixed_table* table = 
				gsl_integration_glfixed_table_alloc(points);
			precision_t u;
			precision_t v;
			precision_t w_i;
			precision_t w_j;
			precision_t result_u;
			precision_t result_v;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				result_u = 0;
				result_v = 0;
				for (size_t i = 0; i < points; i++){
					gsl_integration_glfixed_point(0, 1, i, &u, &w_i, table);
					for (size_t j = 0; j < points; j++){
						gsl_integration_glfixed_point(0, 1, j, &v, &w_j, table);
						result_u += w_j*w_i*integrand_u(u, v, coefficients, 
								coordinates, node_idx);
						result_v += w_j*w_i*integrand_v(u, v, coefficients, 
								coordinates, node_idx);
					}
				}
				gsl_vector_set(result_vector, 
						(size_t)(this->nodes_[node_idx]), 
						gsl_vector_get(result_vector, 
							(size_t)(this->nodes_[node_idx])) + result_u);
				gsl_vector_set(result_vector, 
						(size_t)(this->nodes_[node_idx] + NUM_NODES), 
						gsl_vector_get(result_vector, 
							(size_t)(this->nodes_[node_idx] + NUM_NODES)) - 
						result_v);
			}
			return EXIT_SUCCESS;
		}	

		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				stiff_mat_t& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];
			precision_t z1 = coordinates[n_1][1];
			precision_t z2 = coordinates[n_2][1];
			precision_t z3 = coordinates[n_3][1];
			precision_t k = (r1 + (r2 - r1 + r3 - r1)/3)/2;
			k = k*jacobian(coordinates);
			precision_t s11 = (SIGMA_UR + SIGMA_UZ)*k;
			precision_t s12 = -SIGMA_UR*k;
			precision_t s13 = -SIGMA_UZ*k;
			// s21 = s12;
			precision_t s22 = SIGMA_UR*k;
			// s23 = 0;
			// s31 = s13;
			// s32 = 0;
			precision_t s33 = SIGMA_UZ*k;
			
			gsl_spmatrix_set(&global_stiffness, n_1, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(&global_stiffness, n_1, n_2, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(&global_stiffness, n_1, n_3, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_3) + s13);
			gsl_spmatrix_set(&global_stiffness, n_2, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_2, 
					gsl_spmatrix_get(&global_stiffness, n_2, n_2) + s22);
			gsl_spmatrix_set(&global_stiffness, n_3, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_3, n_1) + s13);
			gsl_spmatrix_set(&global_stiffness, n_3, n_3, 
					gsl_spmatrix_get(&global_stiffness, n_3, n_3) + s33);
			
			n_1 += NUM_NODES;
			n_2 += NUM_NODES;			
			n_3 += NUM_NODES;
			s11 = (SIGMA_VR + SIGMA_VZ)*k;
			s12 = -SIGMA_VR*k;
			s13 = -SIGMA_VZ*k;
			// s21 = s12;
			s22 = SIGMA_VR*k;
			// s23 = 0;
			// s31 = s13;
			// s32 = 0;
			s33 = SIGMA_VZ*k;

			gsl_spmatrix_set(&global_stiffness, n_1, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(&global_stiffness, n_1, n_2, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(&global_stiffness, n_1, n_3, 
					gsl_spmatrix_get(&global_stiffness, n_1, n_3) + s13);
			gsl_spmatrix_set(&global_stiffness, n_2, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_2, 
					gsl_spmatrix_get(&global_stiffness, n_2, n_2) + s22);
			gsl_spmatrix_set(&global_stiffness, n_3, n_1, 
					gsl_spmatrix_get(&global_stiffness, n_3, n_1) + s13);
			gsl_spmatrix_set(&global_stiffness, n_3, n_3, 
					gsl_spmatrix_get(&global_stiffness, n_3, n_3) + s33);

			return EXIT_SUCCESS;
		}
};

template <typename P, typename I>
class ElementBoundary : public Element<P, I>{
	public:
		using typename Element<P, I>::precision_t;
		using typename Element<P, I>::node_t;
		using typename Element<P, I>::stiff_mat_t;
		using typename Element<P, I>::global_vect_t;
	public:
		static node_t NUM_NODES;
		static precision_t RHO_U;
		static precision_t RHO_V;
		static precision_t C_U_AMB;
		static precision_t C_V_AMB;
		bool axis_flag = false;
		ElementBoundary(
				const std::vector<std::vector<precision_t>>& coordinates, 
				const std::vector<node_t>& p_nodes, 
				const std::vector<precision_t>& interior_point)
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
			axis_flag = (coordinates[this->nodes_[0]][0] -
				coordinates[this->nodes_[1]][0]) < 1e-2;
		
		}

	private:
		inline precision_t length(
				const std::vector<std::vector<precision_t>>& coordinates){
			return sqrt(pow((coordinates[this->nodes_[0]][0] - 
						coordinates[this->nodes_[1]][0]), 2) +
				pow((coordinates[this->nodes_[0]][1] - 
							coordinates[this->nodes_[1]][1]), 2));
		}
		
		bool get_orientation(
				const std::vector<std::vector<precision_t>>& coordinates,
				const std::vector<node_t>& p_nodes,
				const std::vector<precision_t>& interior_point){
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
				stiff_mat_t& global_stiffness){
			/*
			Updates the global stiffness matrix
			*/
			if (axis_flag){
				return EXIT_SUCCESS;
			}
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t z1 = coordinates[n_1][1];
			precision_t z2 = coordinates[n_2][1];
			precision_t s11 = RHO_U*length(coordinates)*(r2/4  + r1/12);
			precision_t s12 = RHO_U*length(coordinates)*(r1+r2)/12;
			precision_t s22 = RHO_U*length(coordinates)*(r1/4  + r2/12);

			gsl_spmatrix_set(&global_stiffness, n_1, n_1, gsl_spmatrix_get(
						&global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(&global_stiffness, n_1, n_2, gsl_spmatrix_get(
						&global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_1, gsl_spmatrix_get(
						&global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_2, gsl_spmatrix_get(
						&global_stiffness, n_2, n_2) + s22);

			n_1 += NUM_NODES;
			n_2 += NUM_NODES;
			s11 = RHO_V*length(coordinates)*(r2/4  + r1/12);
			s12 = RHO_V*length(coordinates)*(r1+r2)/12;
			s22 = RHO_V*length(coordinates)*(r1/4  + r2/12);

			gsl_spmatrix_set(&global_stiffness, n_1, n_1, gsl_spmatrix_get(
						&global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(&global_stiffness, n_1, n_2, gsl_spmatrix_get(
						&global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_1, gsl_spmatrix_get(
						&global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(&global_stiffness, n_2, n_2, gsl_spmatrix_get(
						&global_stiffness, n_2, n_2) + s22);
			return EXIT_SUCCESS;
		}

		int update_vector_f(
				const std::vector<std::vector<precision_t>>& coordinates, 
				global_vect_t& vector_f){
			if (axis_flag){
				return EXIT_SUCCESS;
			}
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			gsl_vector_set(&vector_f, n_1, gsl_vector_get(&vector_f, n_1) - 
					RHO_U*length(coordinates)*C_U_AMB*(r1/6 + r2/3));
			gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) - 
					RHO_U*length(coordinates)*C_U_AMB*(r2/6 + r1/3));

			n_1 += NUM_NODES;
			n_2 += NUM_NODES;
			gsl_vector_set(&vector_f, n_1, gsl_vector_get(&vector_f, n_1) - 
					RHO_V*length(coordinates)*C_V_AMB*(r1/6 + r2/3));
			gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) - 
					RHO_V*length(coordinates)*C_V_AMB*(r2/6 + r1/3));
			return EXIT_SUCCESS;
		}
};
}

