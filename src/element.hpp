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
		static precision_t C_U_AMB;
		static precision_t C_V_AMB;
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
		//precision_t slope_1 = 
		//	(coordinates[p_nodes[1]][1] - coordinates[p_nodes[0]][1])/
		//	(coordinates[p_nodes[1]][0] - coordinates[p_nodes[0]][0]);
		//precision_t slope_2 = 
		//	(coordinates[p_nodes[2]][1] - coordinates[p_nodes[1]][1])/
		//	(coordinates[p_nodes[2]][0] - coordinates[p_nodes[1]][0]);
		//return (slope_1 < slope_2);
		return true;
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

		precision_t phi(precision_t epsilon, precision_t eta, int idx){
			assert((idx > 0) && (idx < 4));
			if (idx == 1){
				return phi_1(epsilon, eta);
			}
			else if (idx == 2){
				return phi_2(epsilon, eta);
			}
			else if (idx == 3){
				return phi_3(epsilon, eta);
			}
		}
		
	public:
		precision_t r_u(const gsl_vector* coefficients, 
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
			precision_t c_u = c_u_1*p_1 + c_u_2*p_2+c_u_3*p_3;
			precision_t c_v = c_v_1*p_1 + c_v_2*p_2+c_v_3*p_3;
			return V_MU*c_u/((K_MU + c_u)*(1 + c_v/K_MV));
		}

		precision_t diff_r_u(const gsl_vector* coefficients, int coeff_idx){
			assert((coeff_idx > 0) && (coeff_idx < 7));
			node_t node_num = (coeff_idx - 1)%3;
			precision_t c_u = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[node_num]);
			precision_t c_v = gsl_vector_get(coefficients, 
					(size_t)(this->nodes_[node_num] + NUM_NODES));
			precision_t k = 0;
			if (coeff_idx < 4){
				return ((V_MU*(K_MU + c_u)*(1 + c_v/K_MV) - 
						V_MU*c_u*(1 + c_v/K_MV))/
					pow((K_MU + c_u)*(1 + c_v/K_MV), 2));
			}
			else{
				return (-(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1 + c_v/K_MV), 2)));
			}
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

		precision_t diff_r_v(const gsl_vector* coefficients, int coeff_idx){
			assert((coeff_idx > 0) && (coeff_idx < 7));
			node_t node_num = (coeff_idx - 1)%3;
			precision_t c_u = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[node_num]);
			precision_t c_v = gsl_vector_get(coefficients, 
					(size_t)(this->nodes_[node_num] + NUM_NODES));
			if (coeff_idx < 4){
				return (RESP_Q*(V_MU*(K_MU + c_u)*(1 + c_v/K_MV) - 
						V_MU*c_u*(1 + c_v/K_MV))/
					pow((K_MU + c_u)*(1 + c_v/K_MV), 2) - 
					MAX_FERM_CO2/(K_MFU*pow((1 + c_u/K_MFU), 2)));
			}
			else{
				return (-RESP_Q*(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1 + c_v/K_MV), 2)));
				
			}
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
			assert((node_idx > 0) && (node_idx < 4));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			if (node_idx == 1){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_u(coefficients, 
						u, v)*phi_1(u, v)*jacobian(coordinates);
			}else if (node_idx == 2){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_u(coefficients, 
						u, v)*phi_2(u, v)*jacobian(coordinates);
			}else if (node_idx == 3){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_u(coefficients, 
						u, v)*phi_3(u, v)*jacobian(coordinates);
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
			assert((node_idx > 0) && (node_idx < 4));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			if (node_idx == 1){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_v(coefficients, 
						u, v)*phi_1(u, v)*jacobian(coordinates);
			}else if (node_idx == 2){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_v(coefficients, 
						u, v)*phi_2(u, v)*jacobian(coordinates);
			}else if (node_idx == 3){
				return (r1 + (r2 - r1)*u + (r3 - r1)*v)*r_v(coefficients, 
						u, v)*phi_3(u, v)*jacobian(coordinates);
			}
		}

		precision_t diff_integrand_u(precision_t epsilon, precision_t eta, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx, int coef_idx){
			assert((coef_idx > 0) && (coef_idx < 7));
			assert((node_idx >= 0) && (node_idx < 3));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
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
			precision_t c_u = c_u_1*p_1 + c_u_2*p_2 + c_u_3*p_3;
			precision_t c_v = c_v_1*p_1 + c_v_2*p_2 + c_v_3*p_3;
			precision_t k = 0;
			if (coef_idx < 4){
				k = (V_MU*(K_MU + c_u)*(1 + c_v/K_MV) - 
						V_MU*c_u*(1 + c_v/K_MV))/
					pow((K_MU + c_u)*(1 + c_v/K_MV), 2);
			}
			else{
				k = -(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1 + c_v/K_MV), 2));
				
			}
			return k*phi(epsilon, eta, node_idx+1)*phi(epsilon, eta, 
					(coef_idx - 1)%3 + 1)*jacobian(coordinates)*(
						r1 + (r2 - r1)*epsilon + (r3 - r1)*eta);
		}

		precision_t diff_integrand_v(precision_t epsilon, precision_t eta, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx, int coef_idx){
			assert((coef_idx > 0) && (coef_idx < 7));
			assert((node_idx >= 0) && (node_idx < 3));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
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
			precision_t c_u = c_u_1*p_1 + c_u_2*p_2 + c_u_3*p_3;
			precision_t c_v = c_v_1*p_1 + c_v_2*p_2 + c_v_3*p_3;
			precision_t k = 0;
			if (coef_idx < 4){
				k = RESP_Q*(V_MU*(K_MU + c_u)*(1 + c_v/K_MV) - 
						V_MU*c_u*(1 + c_v/K_MV))/
					pow((K_MU + c_u)*(1 + c_v/K_MV), 2) - 
					MAX_FERM_CO2/(K_MFU*pow((1 + c_u/K_MFU), 2));
			}
			else{
				k = -RESP_Q*(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1 + c_v/K_MV), 2));
				
			}
			return k*phi(epsilon, eta, node_idx+1)*phi(epsilon, eta, 
					(coef_idx - 1)%3 + 1)*jacobian(coordinates)*(
						r1 + (r2 - r1)*epsilon + (r3 - r1)*eta);
		}

		int integrate_non_linear_term_3(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				size_t points, gsl_vector* result_vector){
			precision_t result_u;
			precision_t result_v;
			int phi_idx = 1;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				result_u = (
					integrand_u(0.5, 0, coefficients, coordinates, phi_idx) +
					integrand_u(0, 0.5, coefficients, coordinates, phi_idx) +
					integrand_u(0.5, 0.5, coefficients, coordinates, phi_idx)
					)/6;

				result_v = (
					integrand_v(0.5, 0, coefficients, coordinates, phi_idx) +
					integrand_v(0, 0.5, coefficients, coordinates, phi_idx) +
					integrand_v(0.5, 0.5, coefficients, coordinates, phi_idx)
					)/6;

				gsl_vector_set(result_vector, 
						(size_t)(this->nodes_[node_idx]), 
						gsl_vector_get(result_vector, 
							(size_t)(this->nodes_[node_idx])) + result_u);
				gsl_vector_set(result_vector, 
						(size_t)(this->nodes_[node_idx] + NUM_NODES), 
						gsl_vector_get(result_vector, 
							(size_t)(this->nodes_[node_idx] + NUM_NODES)) - 
						result_v);
				phi_idx++;
			}
			return EXIT_SUCCESS;
		}
		
		int update_with_jacobian(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_matrix* result_matrix){
			precision_t result_uu;
			precision_t result_uv;
			precision_t result_vu;
			precision_t result_vv;
			size_t node_1;
			size_t node_2;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				node_1 = this->nodes_[node_idx];
				for (int coeff_idx = 1; coeff_idx < 4; coeff_idx++){
					node_2 = this->nodes_[coeff_idx - 1];	
					result_uu = (
						diff_integrand_u(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx) + 
						diff_integrand_u(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx) +
						diff_integrand_u(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx)
						)/6;

					result_uv = (
						diff_integrand_u(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx + 3) + 
						diff_integrand_u(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3) +
						diff_integrand_u(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3)
						)/6;

					result_vu = (
						diff_integrand_v(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx) + 
						diff_integrand_v(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx) +
						diff_integrand_v(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx)
						)/6;

					result_vv = (
						diff_integrand_v(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx + 3) + 
						diff_integrand_v(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3) +
						diff_integrand_v(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3)
						)/6;

					gsl_matrix_set(result_matrix, 
							node_1, node_2, gsl_matrix_get(result_matrix, 
								node_1, node_2) + result_uu);
					gsl_matrix_set(result_matrix, 
							node_1, node_2 + NUM_NODES, 
							gsl_matrix_get(result_matrix, node_1, 
								node_2 + NUM_NODES) + result_uv);
					gsl_matrix_set(result_matrix, 
							node_1 + NUM_NODES, node_2, gsl_matrix_get(
								result_matrix, node_1 + NUM_NODES, node_2) - 
							result_vu);
					gsl_matrix_set(result_matrix, 
							node_1 + NUM_NODES, node_2 + NUM_NODES, 
							gsl_matrix_get(result_matrix, node_1 + NUM_NODES, 
								node_2 + NUM_NODES) - result_vv);
				} 
			}
			return EXIT_SUCCESS;
		}
	

/*
 *
 *BEGINING OF NEW APPROACH
 *
 *
 * */	
		
		int update_sp_with_linearized_int_2(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_spmatrix* global_stiffness){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];

			precision_t r_u1 = diff_r_u(coefficients, 1);
			precision_t r_u2 = diff_r_u(coefficients, 2);
			precision_t r_u3 = diff_r_u(coefficients, 3);
			precision_t r_u4 = diff_r_u(coefficients, 4);
			precision_t r_u5 = diff_r_u(coefficients, 5);
			precision_t r_u6 = diff_r_u(coefficients, 6);

			precision_t r_v1 = diff_r_v(coefficients, 1);
			precision_t r_v2 = diff_r_v(coefficients, 2);
			precision_t r_v3 = diff_r_v(coefficients, 3);
			precision_t r_v4 = diff_r_v(coefficients, 4);
			precision_t r_v5 = diff_r_v(coefficients, 5);
			precision_t r_v6 = diff_r_v(coefficients, 6);

			precision_t	s11 = jacobian(coordinates)*(6*r1 + 2*r2 + 2*r3)/120;
			precision_t	s12 = jacobian(coordinates)*(2*r1 + 2*r2 + r3)/120;
			precision_t	s13 = jacobian(coordinates)*(2*r1 + r2 + 2*r3)/120;
			precision_t	s22 = jacobian(coordinates)*(2*r1 + 6*r2 + 2*r3)/120;
			precision_t	s23 = jacobian(coordinates)*(r1 + 2*r2 + 2*r3)/120;
			precision_t	s33 = jacobian(coordinates)*(2*r1 + 2*r2 + 6*r3)/120;
			
			gsl_spmatrix_set(global_stiffness, n_1, n_1, 
					gsl_spmatrix_get(global_stiffness, n_1, n_1) + 
					r_u1*s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2, 
					gsl_spmatrix_get(global_stiffness, n_1, n_2) + 
					r_u2*s12);
			gsl_spmatrix_set(global_stiffness, n_1, n_3, 
					gsl_spmatrix_get(global_stiffness, n_1, n_3) + 
					r_u3*s13);
			gsl_spmatrix_set(global_stiffness, n_2, n_1, 
					gsl_spmatrix_get(global_stiffness, n_2, n_1) + 
					r_u1*s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2, 
					gsl_spmatrix_get(global_stiffness, n_2, n_2) + 
					r_u2*s22);
			gsl_spmatrix_set(global_stiffness, n_2, n_3, 
					gsl_spmatrix_get(global_stiffness, n_2, n_3) + 
					r_u3*s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_1, 
					gsl_spmatrix_get(global_stiffness, n_3, n_1) + 
					r_u1*s13);
			gsl_spmatrix_set(global_stiffness, n_3, n_2, 
					gsl_spmatrix_get(global_stiffness, n_3, n_2) + 
					r_u2*s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_3, 
					gsl_spmatrix_get(global_stiffness, n_3, n_3) + 
					r_u3*s33);
			
			gsl_spmatrix_set(global_stiffness, n_1, n_1 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_1, n_1 + NUM_NODES)
					+ r_u4*s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_1, n_2 + NUM_NODES)
					+ r_u5*s12);
			gsl_spmatrix_set(global_stiffness, n_1, n_3 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_1, n_3 + NUM_NODES)
					+ r_u6*s13);
			gsl_spmatrix_set(global_stiffness, n_2, n_1 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_2, n_1 + NUM_NODES)
					+ r_u4*s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_2, n_2 + NUM_NODES)
					+ r_u5*s22);
			gsl_spmatrix_set(global_stiffness, n_2, n_3 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_2, n_3 + NUM_NODES)
					+ r_u6*s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_1 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_3, n_1 + NUM_NODES)
					+ r_u4*s13);
			gsl_spmatrix_set(global_stiffness, n_3, n_2 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_3, n_2 + NUM_NODES)
					+ r_u5*s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_3 + NUM_NODES, 
					gsl_spmatrix_get(global_stiffness, n_3, n_3 + NUM_NODES)
					+ r_u6*s33);
			
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_1, 
					gsl_spmatrix_get(global_stiffness, n_1 + NUM_NODES, n_1)
					- r_v1*s11);
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_2, 
					gsl_spmatrix_get(global_stiffness, n_1 + NUM_NODES, n_2)
					- r_v2*s12);
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_3, 
					gsl_spmatrix_get(global_stiffness, n_1 + NUM_NODES, n_3)
					- r_v3*s13);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_1, 
					gsl_spmatrix_get(global_stiffness, n_2 + NUM_NODES, n_1)
					- r_v1*s12);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_2, 
					gsl_spmatrix_get(global_stiffness, n_2 + NUM_NODES, n_2)
					- r_v2*s22);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_3, 
					gsl_spmatrix_get(global_stiffness, n_2 + NUM_NODES, n_3)
					- r_v3*s23);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_1, 
					gsl_spmatrix_get(global_stiffness, n_3 + NUM_NODES, n_1)
					- r_v1*s13);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_2, 
					gsl_spmatrix_get(global_stiffness, n_3 + NUM_NODES, n_2)
					- r_v2*s23);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_3, 
					gsl_spmatrix_get(global_stiffness, n_3 + NUM_NODES, n_3)
					- r_v3*s33);
			
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_1 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_1 + 
						NUM_NODES, n_1 + NUM_NODES) - r_v4*s11);
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_2 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_1 + 
						NUM_NODES, n_2 + NUM_NODES)	- r_v5*s12);
			gsl_spmatrix_set(global_stiffness, n_1 + NUM_NODES, n_3 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_1 + 
						NUM_NODES, n_3 + NUM_NODES)	- r_v6*s13);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_1 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_2 + 
						NUM_NODES, n_1 + NUM_NODES) - r_v4*s12);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_2 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_2 + 
						NUM_NODES, n_2 + NUM_NODES)	- r_v5*s22);
			gsl_spmatrix_set(global_stiffness, n_2 + NUM_NODES, n_3 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_2 + 
						NUM_NODES, n_3 + NUM_NODES)	- r_v6*s23);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_1 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_3 + 
						NUM_NODES, n_1 + NUM_NODES) - r_v4*s13);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_2 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_3 + 
						NUM_NODES, n_2 + NUM_NODES)	- r_v5*s23);
			gsl_spmatrix_set(global_stiffness, n_3 + NUM_NODES, n_3 + 
					NUM_NODES, gsl_spmatrix_get(global_stiffness, n_3 + 
						NUM_NODES, n_3 + NUM_NODES)	- r_v6*s33);
			return EXIT_SUCCESS;
		}


		
		int update_f_vector_with_linearized_integral_2(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates, 
				global_vect_t& vector_f){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];

			precision_t d_r_u1 = diff_r_u(coefficients, 1);
			precision_t d_r_u2 = diff_r_u(coefficients, 2);
			precision_t d_r_u3 = diff_r_u(coefficients, 3);
			precision_t d_r_u4 = diff_r_u(coefficients, 4);
			precision_t d_r_u5 = diff_r_u(coefficients, 5);
			precision_t d_r_u6 = diff_r_u(coefficients, 6);

			precision_t d_r_v1 = diff_r_v(coefficients, 1);
			precision_t d_r_v2 = diff_r_v(coefficients, 2);
			precision_t d_r_v3 = diff_r_v(coefficients, 3);
			precision_t d_r_v4 = diff_r_v(coefficients, 4);
			precision_t d_r_v5 = diff_r_v(coefficients, 5);
			precision_t d_r_v6 = diff_r_v(coefficients, 6);

			precision_t r_u1 = r_u(coefficients, 0, 0);
			precision_t r_u2 = r_u(coefficients, 1, 0);
			precision_t r_u3 = r_u(coefficients, 0, 1);
			precision_t r_v1 = r_v(coefficients, 0, 0);
			precision_t r_v2 = r_v(coefficients, 1, 0);
			precision_t r_v3 = r_v(coefficients, 0, 1);

			precision_t	s11 = jacobian(coordinates)*(6*r1 + 2*r2 + 2*r3)/120;
			precision_t	s12 = jacobian(coordinates)*(2*r1 + 2*r2 + r3)/120;
			precision_t	s13 = jacobian(coordinates)*(2*r1 + r2 + 2*r3)/120;
			precision_t	s22 = jacobian(coordinates)*(2*r1 + 6*r2 + 2*r3)/120;
			precision_t	s23 = jacobian(coordinates)*(r1 + 2*r2 + 2*r3)/120;
			precision_t	s33 = jacobian(coordinates)*(2*r1 + 2*r2 + 6*r3)/120;
			
			precision_t f_u_1 = 
				s11*(r_u1 - d_r_u1*C_U_AMB - d_r_u4*C_V_AMB) + 
				s12*(r_u2 - d_r_u2*C_U_AMB - d_r_u5*C_V_AMB) + 
				s13*(r_u3 - d_r_u3*C_U_AMB - d_r_u6*C_V_AMB);
			
			precision_t f_u_2 = 
				s12*(r_u1 - d_r_u1*C_U_AMB - d_r_u4*C_V_AMB) + 
				s22*(r_u2 - d_r_u2*C_U_AMB - d_r_u5*C_V_AMB) + 
				s23*(r_u3 - d_r_u3*C_U_AMB - d_r_u6*C_V_AMB);
			
			precision_t f_u_3 = 
				s13*(r_u1 - d_r_u1*C_U_AMB - d_r_u4*C_V_AMB) + 
				s23*(r_u2 - d_r_u2*C_U_AMB - d_r_u5*C_V_AMB) + 
				s33*(r_u3 - d_r_u3*C_U_AMB - d_r_u6*C_V_AMB);
			
			precision_t f_v_1 = 
				s11*(r_v1 - d_r_v1*C_U_AMB - d_r_v4*C_V_AMB) + 
				s12*(r_v2 - d_r_v2*C_U_AMB - d_r_v5*C_V_AMB) + 
				s13*(r_v3 - d_r_v3*C_U_AMB - d_r_v6*C_V_AMB);
			
			precision_t f_v_2 = 
				s12*(r_v1 - d_r_v1*C_U_AMB - d_r_v4*C_V_AMB) + 
				s22*(r_v2 - d_r_v2*C_U_AMB - d_r_v5*C_V_AMB) + 
				s23*(r_v3 - d_r_v3*C_U_AMB - d_r_v6*C_V_AMB);
			
			precision_t f_v_3 = 
				s13*(r_v1 - d_r_v1*C_U_AMB - d_r_v4*C_V_AMB) + 
				s23*(r_v2 - d_r_v2*C_U_AMB - d_r_v5*C_V_AMB) + 
				s33*(r_v3 - d_r_v3*C_U_AMB - d_r_v6*C_V_AMB);
			
			gsl_vector_set(&vector_f, n_1,
					gsl_vector_get(&vector_f, n_1) 
					+ f_u_1);
			gsl_vector_set(&vector_f, n_2,
					gsl_vector_get(&vector_f, n_2) 
					+ f_u_2);
			gsl_vector_set(&vector_f, n_3,
					gsl_vector_get(&vector_f, n_3) 
					+ f_u_3);
			gsl_vector_set(&vector_f, n_1 + NUM_NODES,
					gsl_vector_get(&vector_f, n_1 + NUM_NODES) 
					- f_v_1);
			gsl_vector_set(&vector_f, n_2 + NUM_NODES,
					gsl_vector_get(&vector_f, n_2 + NUM_NODES) 
					- f_v_2);
			gsl_vector_set(&vector_f, n_3 + NUM_NODES,
					gsl_vector_get(&vector_f, n_3 + NUM_NODES) 
					- f_v_3);
			
			return EXIT_SUCCESS;	
		}

/*
 *
 *END OF NEW APPROACH
 *
 *
 * */





		int update_sp_with_jacobian(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_spmatrix* result_matrix){
			precision_t result_uu;
			precision_t result_uv;
			precision_t result_vu;
			precision_t result_vv;
			size_t node_1;
			size_t node_2;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				for (int coeff_idx = 1; coeff_idx < 4; coeff_idx++){
				node_1 = this->nodes_[node_idx];	
				node_2 = this->nodes_[coeff_idx - 1];	
				result_uu = (
					diff_integrand_u(0.5, 0, coefficients, coordinates,
						node_idx, coeff_idx) + 
					diff_integrand_u(0, 0.5, coefficients, coordinates,
						node_idx, coeff_idx) +
					diff_integrand_u(0.5, 0.5, coefficients, coordinates,
						node_idx, coeff_idx)
					)/6;

				result_uv = (
					diff_integrand_u(0.5, 0, coefficients, coordinates,
						node_idx, coeff_idx + 3) + 
					diff_integrand_u(0, 0.5, coefficients, coordinates,
						node_idx, coeff_idx + 3) +
					diff_integrand_u(0.5, 0.5, coefficients, coordinates,
						node_idx, coeff_idx + 3)
					)/6;

				result_vu = (
					diff_integrand_v(0.5, 0, coefficients, coordinates,
						node_idx, coeff_idx) + 
					diff_integrand_v(0, 0.5, coefficients, coordinates,
						node_idx, coeff_idx) +
					diff_integrand_v(0.5, 0.5, coefficients, coordinates,
						node_idx, coeff_idx)
					)/6;

				result_vv = (
					diff_integrand_v(0.5, 0, coefficients, coordinates,
						node_idx, coeff_idx + 3) + 
					diff_integrand_v(0, 0.5, coefficients, coordinates,
						node_idx, coeff_idx + 3) +
					diff_integrand_v(0.5, 0.5, coefficients, coordinates,
						node_idx, coeff_idx + 3)
					)/6;

				gsl_spmatrix_set(result_matrix, 
						node_1, node_2, gsl_spmatrix_get(result_matrix, 
							node_1, node_2) + result_uu);
				gsl_spmatrix_set(result_matrix, 
						node_1, node_2 + NUM_NODES, 
						gsl_spmatrix_get(result_matrix, node_1, 
							node_2 + NUM_NODES) + result_uv);
				gsl_spmatrix_set(result_matrix, 
						node_1 + NUM_NODES, node_2, gsl_spmatrix_get(
							result_matrix, node_1 + NUM_NODES, node_2) - 
						result_vu);
				gsl_spmatrix_set(result_matrix, 
						node_1 + NUM_NODES, node_2 + NUM_NODES, 
						gsl_spmatrix_get(result_matrix, node_1 + NUM_NODES, 
							node_2 + NUM_NODES) - result_vv);
				} 
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
			precision_t k = (r1 + r2 + r3)/6;
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
		int update_stiffness_with_linearized_integral(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				stiff_mat_t& global_stiffness){
			update_sp_with_jacobian(coefficients, 
					coordinates, &global_stiffness);
		
		}
		//int update_stiffness_with_linearized_integral(
		//		const std::vector<std::vector<precision_t>>& coordinates,
		//		stiff_mat_t& global_stiffness){

		//	node_t n_1 = this->nodes_[0];
		//	node_t n_2 = this->nodes_[1];
		//	node_t n_3 = this->nodes_[2];
		//	precision_t r1 = coordinates[n_1][0];
		//	precision_t r2 = coordinates[n_2][0];
		//	precision_t r3 = coordinates[n_3][0];
		//	


// THIS //IS FOR Cu's

		//	precision_t k = K_MU*K_MV*V_MU/( std::pow(C_U_AMB + K_MU,2)*(C_V_AMB + K_MV) );
		//	k = k*jacobian(coordinates);
		//	precision_t s11 = (3*r1 + r2 + r3)*k/60;
		//	precision_t s12 = (2*r1 + 2*r2 + r3)*k/120;
		//	precision_t s13 = (2*r1 + r2 + 2*r3)*k/120;

		//	precision_t s21 = (2*r1 + 2*r2 + r3)*k/120;
		//	precision_t s22 = (r1 + 3*r2 + r3)*k/60;
		//	precision_t s23 = (r1 + 2*r2 + 2*r3)*k/120;
		//	
		//	precision_t s31 = (2*r1 + r2 + 2*r3)*k/120;
		//	precision_t s32 = (r1 + 2*r2 + 2*r3)*k/120;
		//	precision_t s33 = (r1 + r2 + 3*r3)*k/60;

		//	gsl_spmatrix_set(&global_stiffness, n_1, n_1, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_1)
		//			+ s11);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_2, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_2)
		//			+ s12);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_3, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_3)
		//			+ s13);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_1, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_1)
		//			+ s21);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_2, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_2)
		//			+ s22);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_3, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_3)
		//			+ s23);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_1, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_1)
		//			+ s31);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_2, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_2)
		//			+ s32);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_3, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_3)
		//			+ s33);


// THIS //IS FOR COUPLING BETWEEN  Cv's against Cu's
		//	k = -C_U_AMB*V_MU*K_MV/( (C_U_AMB + K_MU)*std::pow(C_V_AMB + K_MV,2) );
		//	k = k*jacobian(coordinates);

		//	s11 = (3*r1 + r2 + r3)*k/60;
		//	s12 = (2*r1 + 2*r2 + r3)*k/120;
		//	s13 = (2*r1 + r2 + 2*r3)*k/120;

		//	s21 = (2*r1 + 2*r2 + r3)*k/120;
		//	s22 = (r1 + 3*r2 + r3)*k/60;
		//	s23 = (r1 + 2*r2 + 2*r3)*k/120;
		//	
		//	s31 = (2*r1 + r2 + 2*r3)*k/120;
		//	s32 = (r1 + 2*r2 + 2*r3)*k/120;
		//	s33 = (r1 + r2 + 3*r3)*k/60;

		//	gsl_spmatrix_set(&global_stiffness, n_1, n_1 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_1 + NUM_NODES)
		//			+ s11);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_2 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_2 + NUM_NODES)
		//			+ s12);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_3 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_3 + NUM_NODES)
		//			+ s13);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_1 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_1 + NUM_NODES)
		//			+ s21);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_2 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_2 + NUM_NODES)
		//			+ s22);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_3 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_3 + NUM_NODES)
		//			+ s23);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_1 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_1 + NUM_NODES)
		//			+ s31);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_2 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_2 + NUM_NODES)
		//			+ s32);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_3 + NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_3 + NUM_NODES)
		//			+ s33);

// THIS //IS FOR Cv's 
		//	n_1 += NUM_NODES;
		//	n_2 += NUM_NODES;
		//	n_3 += NUM_NODES;

		//	k = -C_U_AMB*K_MV*V_MU*RESP_Q/( (C_U_AMB + K_MU)*std::pow(C_V_AMB + K_MV,2) );
		//	k = k*jacobian(coordinates);
		//	s11 = (3*r1 + r2 + r3)*k/60;
		//	s12 = (2*r1 + 2*r2 + r3)*k/120;
		//	s13 = (2*r1 + r2 + 2*r3)*k/120;

		//	s21 = (2*r1 + 2*r2 + r3)*k/120;
		//	s22 = (r1 + 3*r2 + r3)*k/60;
		//	s23 = (r1 + 2*r2 + 2*r3)*k/120;
		//	
		//	s31 = (2*r1 + r2 + 2*r3)*k/120;
		//	s32 = (r1 + 2*r2 + 2*r3)*k/120;
		//	s33 = (r1 + r2 + 3*r3)*k/60;

		//	gsl_spmatrix_set(&global_stiffness, n_1, n_1, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_1)
		//			- s11);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_2,
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_2)
		//			- s12);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_3,
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_3)
		//			- s13);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_1,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_1)
		//			- s21);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_2,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_2)
		//			- s22);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_3,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_3)
		//			- s23);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_1,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_1)
		//			- s31);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_2,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_2)
		//			- s32);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_3,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_3)
		//			- s33);

// THIS //IS FOR COUPLING BETWEEN  Cu's against Cv's

		//	// Same as first part
		//	
		//	precision_t numerator = K_MV*(-std::pow(C_U_AMB,2)*K_MFU*MAX_FERM_CO2 + std::pow(C_U_AMB,2)*K_MU*V_MU*RESP_Q - 2*C_U_AMB*K_MFU*K_MU*MAX_FERM_CO2 + 2*C_U_AMB*K_MFU*K_MU*V_MU*RESP_Q + std::pow(K_MFU,2)*K_MU*V_MU*RESP_Q - K_MFU*std::pow(K_MU,2)*MAX_FERM_CO2) ;
		//	//precision_t numerator = -C_U_AMB*K_MFU*MAX_FERM_CO2 + C_U_AMB*K_MU*V_MU*RESP_Q \
		//						 - 2*K_MFU*K_MU*MAX_FERM_CO2 + 2*K_MFU*K_MU*V_MU*RESP_Q;
		//	
		//	precision_t denominator = C_U_AMB*(std::pow(C_U_AMB,3)*C_V_AMB + std::pow(C_U_AMB,3)*K_MV + 2*std::pow(C_U_AMB,2)*K_MFU*K_MV + 2*std::pow(C_U_AMB,2)*K_MU*K_MV + C_U_AMB*std::pow(K_MFU,2)*K_MV + 4*C_U_AMB*K_MFU*K_MU*K_MV + C_U_AMB*std::pow(K_MU,2)*K_MV + 2*std::pow(K_MFU,2)*K_MU*K_MV + 2*K_MFU*std::pow(K_MU,2)*K_MV);
		//	//precision_t denominator = C_U_AMB*( std::pow(C_U_AMB,2) + 2*C_U_AMB*K_MFU + 2*C_U_AMB*K_MU + 4*K_MFU*K_MU + std::pow(K_MU,2) );
		//	k = numerator/denominator;
		//	k = k*jacobian(coordinates);
		//	s11 = (3*r1 + r2 + r3)*k/60;
		//	s12 = (2*r1 + 2*r2 + r3)*k/120;
		//	s13 = (2*r1 + r2 + 2*r3)*k/120;

		//	s21 = (2*r1 + 2*r2 + r3)*k/120;
		//	s22 = (r1 + 3*r2 + r3)*k/60;
		//	s23 = (r1 + 2*r2 + 2*r3)*k/120;
		//	
		//	s31 = (2*r1 + r2 + 2*r3)*k/120;
		//	s32 = (r1 + 2*r2 + 2*r3)*k/120;
		//	s33 = (r1 + r2 + 3*r3)*k/60;

		//	gsl_spmatrix_set(&global_stiffness, n_1, n_1 - NUM_NODES, 
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_1 - NUM_NODES)
		//			- s11);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_2 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_2 - NUM_NODES)
		//			- s12);
		//	gsl_spmatrix_set(&global_stiffness, n_1, n_3 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_1, n_3 - NUM_NODES)
		//			- s13);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_1 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_1 - NUM_NODES)
		//			- s21);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_2 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_2 - NUM_NODES)
		//			- s22);
		//	gsl_spmatrix_set(&global_stiffness, n_2, n_3 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_2, n_3 - NUM_NODES)
		//			- s23);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_1 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_1 - NUM_NODES)
		//			- s31);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_2 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_2 - NUM_NODES)
		//			- s32);
		//	gsl_spmatrix_set(&global_stiffness, n_3, n_3 - NUM_NODES,
		//			gsl_spmatrix_get(&global_stiffness, n_3, n_3 - NUM_NODES)
		//			- s33);



		//	return EXIT_SUCCESS;	
		//}
		
		//int update_f_vector_with_linearized_integral(
		//		const std::vector<std::vector<precision_t>>& coordinates, 
		//		global_vect_t& vector_f){
		//	node_t n_1 = this->nodes_[0];
		//	node_t n_2 = this->nodes_[1];
		//	node_t n_3 = this->nodes_[2];
		//	precision_t r1 = coordinates[n_1][0];
		//	precision_t r2 = coordinates[n_2][0];
		//	precision_t r3 = coordinates[n_3][0];



		//	precision_t numerator = C_U_AMB*K_MV*V_MU*(2*C_U_AMB*C_V_AMB + C_U_AMB*K_MV + C_V_AMB*K_MU);
		//	precision_t denominator = 24*std::pow(C_U_AMB*C_V_AMB + C_U_AMB*K_MV + C_V_AMB*K_MU + K_MU*K_MV,2);

		//	precision_t k = numerator/denominator;
		//	k = k*jacobian(coordinates);
		//	
		//	precision_t val_1 = k*(2*r1 + r2 + r3);
		//	precision_t val_2 = k*(r1 + 2*r2 + r3);
		//	precision_t val_3 = k*(r1 + r2 + 2*r3);

		//	gsl_vector_set(&vector_f, n_1, gsl_vector_get(&vector_f, n_1) + 
		//			val_1);
		//	gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) + 
		//			val_2);
		//	gsl_vector_set(&vector_f, n_3, gsl_vector_get(&vector_f, n_3) + 
		//			val_3);

		//	n_1 += NUM_NODES;
		//	n_2 += NUM_NODES;
		//	n_3 += NUM_NODES;


		//	numerator = 2*std::pow(C_U_AMB,2)*C_V_AMB*V_MU*RESP_Q + std::pow(C_U_AMB,2)*K_MV*V_MU*RESP_Q + 2*C_U_AMB*K_MFU*K_MV*MAX_FERM_CO2 \
		//				+ 2*C_U_AMB*K_MFU*K_MV*V_MU*RESP_Q + 4*K_MFU*K_MU*K_MV*MAX_FERM_CO2;

		//	denominator = 24*K_MV*(std::pow(C_U_AMB,2) + 2*C_U_AMB*K_MFU + 2*C_U_AMB*K_MU + 4*K_MFU*K_MU + std::pow(K_MU,2) );

		//	k = numerator/denominator;
		//	k= k*jacobian(coordinates);
		//
		//	val_1 = k*(2*r1 + r2 + r3);
		//	val_2 = k*(r1 + 2*r2 + r3);
		//	val_3 = k*(r1 + r2 + 2*r3);


		//	gsl_vector_set(&vector_f, n_1, gsl_vector_get(&vector_f, n_1) - 
		//			val_1);
		//	gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) - 
		//			val_2);
		//	gsl_vector_set(&vector_f, n_3, gsl_vector_get(&vector_f, n_3) - 
		//			val_3);
		//	return EXIT_SUCCESS;	
		//}
		
		int update_f_vector_with_linearized_integral(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates, 
				global_vect_t& vector_f){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];

			integrate_non_linear_term_3(coefficients, coordinates, 1, 
					&vector_f);
			
			precision_t result_uu;
			precision_t result_uv;
			precision_t result_vu;
			precision_t result_vv;
			size_t node_1;
			size_t node_2;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				for (int coeff_idx = 1; coeff_idx < 4; coeff_idx++){
					node_1 = this->nodes_[node_idx];	
					node_2 = this->nodes_[coeff_idx - 1];	
					result_uu = (
						diff_integrand_u(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx) + 
						diff_integrand_u(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx) +
						diff_integrand_u(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx)
						)/6;

					result_uv = (
						diff_integrand_u(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx + 3) + 
						diff_integrand_u(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3) +
						diff_integrand_u(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3)
						)/6;

					result_vu = (
						diff_integrand_v(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx) + 
						diff_integrand_v(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx) +
						diff_integrand_v(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx)
						)/6;

					result_vv = (
						diff_integrand_v(0.5, 0, coefficients, coordinates,
							node_idx, coeff_idx + 3) + 
						diff_integrand_v(0, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3) +
						diff_integrand_v(0.5, 0.5, coefficients, coordinates,
							node_idx, coeff_idx + 3)
						)/6;

					gsl_vector_set(&vector_f, node_1,
							gsl_vector_get(&vector_f, node_1) 
							- (result_uu*C_U_AMB + result_uv*C_V_AMB));
					gsl_vector_set(&vector_f, node_1 + NUM_NODES,
							gsl_vector_get(&vector_f, node_1 + NUM_NODES) 
							+ (result_vu*C_U_AMB + result_vv*C_V_AMB));
					} 
			}
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
			axis_flag = (abs((coordinates[this->nodes_[0]][0] -
				coordinates[this->nodes_[1]][0])) < 1e-6) && 
				(coordinates[this->nodes_[0]][0] < 1e-6);
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
			if (ang2 < 0 && ang1 > 0){
				return true;
			}
			return ang2 > ang1;
		}
		
	public:

		// Access
		const std::vector<node_t>& nodes(){
			return this->nodes_;
		}
		
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
			precision_t s11 = RHO_U*length(coordinates)*(r1/4  + r2/12);
			precision_t s12 = RHO_U*length(coordinates)*(r1 + r2)/12;
			precision_t s22 = RHO_U*length(coordinates)*(r2/4  + r1/12);

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
			s11 = RHO_V*length(coordinates)*(r1/4  + r2/12);
			s12 = RHO_V*length(coordinates)*(r1 + r2)/12;
			s22 = RHO_V*length(coordinates)*(r2/4  + r1/12);

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
					RHO_U*length(coordinates)*C_U_AMB*(r1/3 + r2/6));
			gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) - 
					RHO_U*length(coordinates)*C_U_AMB*(r2/3 + r1/6));

			n_1 += NUM_NODES;
			n_2 += NUM_NODES;
			gsl_vector_set(&vector_f, n_1, gsl_vector_get(&vector_f, n_1) - 
					RHO_V*length(coordinates)*C_V_AMB*(r1/3 + r2/6));
			gsl_vector_set(&vector_f, n_2, gsl_vector_get(&vector_f, n_2) - 
					RHO_V*length(coordinates)*C_V_AMB*(r2/3 + r1/6));
			return EXIT_SUCCESS;
		}
};
}
