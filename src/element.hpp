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

template <typename P>
class Element{
	public:
		typedef P precision_t;
		typedef size_t node_t;
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
				stiff_mat_t* global_stiffness) = 0;
};

template <typename P>
class ElementTriangular : public Element<P>{
	public:
		using typename Element<P>::precision_t;
		using typename Element<P>::node_t;
		using typename Element<P>::stiff_mat_t;
		using typename Element<P>::global_vect_t;
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
			ElementTriangular<precision_t>::NUM_ELM++;
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

		bool is_point_interior(precision_t r, precision_t z, 
				const std::vector<std::vector<precision_t>>& coordinates){
			/*
			 * Checks if a point r, z is interior to the element (function used
			 * for debugging purposes)
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			precision_t z1 = coordinates[this->nodes_[0]][1];
			precision_t z2 = coordinates[this->nodes_[1]][1];
			precision_t z3 = coordinates[this->nodes_[2]][1];
			precision_t area_1 = abs((r1 - r)*(z2 - z)-(r2 - r)*(z1 - z))/2;
			precision_t area_2 = abs((r1 - r)*(z3 - z)-(r3 - r)*(z1 - z))/2;
			precision_t area_3 = abs((r2 - r)*(z3 - z)-(r3 - r)*(z2 - z))/2;
			precision_t tot_a = abs(abs((r2 - r1)*(z3 - z1) - 
						(r3 - r1)*(z2 - z1)))/2;
			return abs(tot_a - (area_1 + area_2 + area_3))/tot_a < 1e-10;
		}

		int get_cons_at(precision_t r, precision_t z,
			const std::vector<std::vector<precision_t>>& coordinates,
			gsl_vector* coefficients, std::vector<precision_t>& res){
			/*
			 * gets the concentration at value r, z (function used
			 * for debugging purposes)
			*/
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t c_u_1 = gsl_vector_get(coefficients, n_1);
			precision_t c_u_2 = gsl_vector_get(coefficients, n_2);
			precision_t c_u_3 = gsl_vector_get(coefficients, n_3);
			precision_t c_v_1 = gsl_vector_get(coefficients, n_1);
			precision_t c_v_2 = gsl_vector_get(coefficients, n_2);
			precision_t c_v_3 = gsl_vector_get(coefficients, n_3);
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];
			precision_t z1 = coordinates[n_1][1];
			precision_t z2 = coordinates[n_2][1];
			precision_t z3 = coordinates[n_3][1];
			precision_t epsilon = (-r*(z1-z3) + z*(r1-r3))/
				((r1-r2)*(z1-z3) - (r1-r3)*(z1-z2));
			precision_t eta = (r*(z1-z2) - z*(r1-r2))/
				((r1-r2)*(z1-z3) - (r1-r3)*(z1-z2));
			res[0] = phi_1(epsilon, eta)*c_u_1 + phi_2(epsilon, eta)*c_u_2 +
				phi_3(epsilon, eta)*c_u_3;
			res[1] = phi_1(epsilon, eta)*c_v_1 + phi_2(epsilon, eta)*c_v_2 +
				phi_3(epsilon, eta)*c_v_3;
			return EXIT_SUCCESS;
		}
		
		bool has_nodes(std::vector<node_t>& nodes){
			bool has_flag = true;
			for (auto nod : nodes){
				if (std::find(this->nodes_.begin(), this->nodes_.end(), 
							nod) == this->nodes_.end()){
					has_flag = false;
					break;
				}
			}
			return has_flag;
		}

		precision_t phi_1(precision_t epsilon, precision_t eta){
			/*
			 * returns the value of phi_1 evaluated at 
			 * possition epsilon, eta of the local coordinates 
			*/
			return 1 - epsilon - eta;
		}
		
		precision_t phi_2(precision_t epsilon, precision_t eta){
			/*
			 * returns the value of phi_2 evaluated at 
			 * possition epsilon, eta of the local coordinates 
			*/
			return epsilon;
		}
		
		precision_t phi_3(precision_t epsilon, precision_t eta){
			/*
			 * returns the value of phi_3 evaluated at 
			 * possition epsilon, eta of the local coordinates 
			*/
			return eta;
		}

		precision_t phi(precision_t epsilon, precision_t eta, int idx){
			/*
			 * returns the value of phi_idx evaluated at 
			 * possition epsilon, eta of the local coordinates 
			*/
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
			else{
				assert(1 == 0);
				return EXIT_FAILURE;
			}
		}
		
	public:
		precision_t r_u(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			/*
			 * returns the value of Ru(Cu,Cv) evaluated at coefficients, and in
			 * possition epsilon, eta of the local coordinates 
			*/
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
			return V_MU*c_u/((K_MU + c_u)*(1.0 + c_v/K_MV));
		}

		precision_t diff_r_u(const gsl_vector* coefficients, int coeff_idx){
			/*
			 * returns the derivative of Ru(Cu,Cv) whit respect to 
			 * C_{coeff_idx}. Where coeff_idx maps is {1,2,3} corresponding to
			 * the coefficient C_u in node_1, node_2, node_3, and {4,5,6} for
			 * the coefficients C_v in node_1, node_2, node_3.
			*/
			assert((coeff_idx > 0) && (coeff_idx < 7));
			node_t node_num = (coeff_idx - 1)%3;
			precision_t c_u = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[node_num]);
			precision_t c_v = gsl_vector_get(coefficients, 
					(size_t)(this->nodes_[node_num] + NUM_NODES));
			precision_t k = 0;
			if (coeff_idx < 4){
				return ((V_MU*(K_MU + c_u)*(1.0 + c_v/K_MV) - 
						V_MU*c_u*(1.0 + c_v/K_MV))/
					pow((K_MU + c_u)*(1.0 + c_v/K_MV), 2));
			}
			else{
				return (-(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1.0 + c_v/K_MV), 2)));
			}
		}

		precision_t r_v(const gsl_vector* coefficients, 
				precision_t epsilon, precision_t eta){
			/*
			 * returns the valuue of Rv(Cu,Cv) evaluated at coefficients, and in
			 * possition epsilon, eta of the local coordinates 
			*/
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
				MAX_FERM_CO2/(1.0 + c_u/K_MFU);
		}

		precision_t diff_r_v(const gsl_vector* coefficients, int coeff_idx){
			/*
			 * returns the derivative of Rv(Cu,Cv) whit respect to 
			 * C_{coeff_idx}. Where coeff_idx maps is {1,2,3} corresponding to
			 * the coefficient C_u in node_1, node_2, node_3, and {4,5,6} for
			 * the coefficients C_v in node_1, node_2, node_3.
			*/
			assert((coeff_idx > 0) && (coeff_idx < 7));
			node_t node_num = (coeff_idx - 1)%3;
			precision_t c_u = gsl_vector_get(coefficients, 
					(size_t)this->nodes_[node_num]);
			precision_t c_v = gsl_vector_get(coefficients, 
					(size_t)(this->nodes_[node_num] + NUM_NODES));
			if (coeff_idx < 4){
				return (RESP_Q*(V_MU*(K_MU + c_u)*(1.0 + c_v/K_MV) - 
						V_MU*c_u*(1.0 + c_v/K_MV))/
					pow((K_MU + c_u)*(1.0 + c_v/K_MV), 2) - 
					MAX_FERM_CO2/(K_MFU*pow((1.0 + c_u/K_MFU), 2)));
			}
			else{
				return (-RESP_Q*(V_MU*c_u)*(K_MU + c_u)/
					(K_MV*pow((K_MU + c_u)*(1.0 + c_v/K_MV), 2)));
				
			}
		}

		inline precision_t jacobian(
				const std::vector<std::vector<precision_t>>& coordinates){
			/*
			 * returns the jacobian of the element, when transformed to a
			 * rectangle triangle with minnor sides euqal to 1. The jacobian is
			 * related to the area of the element as J/2 = area
			*/
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			precision_t z1 = coordinates[this->nodes_[0]][1];
			precision_t z2 = coordinates[this->nodes_[1]][1];
			precision_t z3 = coordinates[this->nodes_[2]][1];
			return abs((r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1));
		}	

		precision_t integrand_u(precision_t u, precision_t v, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx){
			/*
			 * returns the integrand r*R_u(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a triangular domain with boundaries [0,1],
			 * [0,1] and [0,0],
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
			else{
				assert(1 == 0);
				return EXIT_FAILURE;
			}
		}
		
		precision_t integrand_v(precision_t u, precision_t v, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx){
			/*
			 * returns the integrand r*R_v(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a triangular domain with boundaries [0,1],
			 * [0,1] and [0,0],
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
			else{
				assert(1 == 0);
				return EXIT_FAILURE;
			}
		}
		
		precision_t integrand_analytical_u(precision_t u, precision_t v, 
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx, precision_t cu_2){
			/*
			 * returns the integrand r*R_u(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a triangular domain with nodes at [0,1],
			 * [0,0] and [0,1], Where R_u is computed using the analitical
			 * solution proposed.
  			*/
			assert((node_idx > 0) && (node_idx < 4));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			return 6.0*cu_2*(r1 + (r2 - r1)*u + (r3 - r1)*v)*SIGMA_UR*
				jacobian(coordinates)*phi(u, v, node_idx);
		}
		
		precision_t integrand_analytical_v(precision_t u, precision_t v, 
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx, precision_t cv_2){
			/*
			 * returns the integrand r*R_v(C_u,C_v)*phi_j(r,z) for the domain 
			 * omega transformed to a triangular domain with nodes at [0,1],
			 * [0,0] and [0,1], Where R_v is computed using the analitical
			 * solution proposed.
  			*/
			assert((node_idx > 0) && (node_idx < 4));
			precision_t r1 = coordinates[this->nodes_[0]][0];
			precision_t r2 = coordinates[this->nodes_[1]][0];
			precision_t r3 = coordinates[this->nodes_[2]][0];
			return -6.0*cv_2*(r1 + (r2 - r1)*u + (r3 - r1)*v)*SIGMA_VR*
				jacobian(coordinates)*phi(u, v, node_idx);
		}

		precision_t diff_integrand_u(precision_t epsilon, precision_t eta, 
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				node_t node_idx, int coef_idx){
			/*
			 * Computes the value of the differenciated integrand of the u
			 * equations, node_idx represents the row to wich the integrand
			 * belongs, and coeff_idx the number of the coefficient that
			 * differentiates the equation
			*/
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
			/*
			 * Computes the value of the differenciated integrand of the v
			 * equations, node_idx represents the row to wich the integrand
			 * belongs, and coeff_idx the number of the coefficient that
			 * differentiates the equation
			*/
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

		int integrate_non_linear_term(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_vector* result_vector){
			/*
			 * Computes the value of the nonlinear function H at the value
			 * coefficients.
			*/
			precision_t result_u;
			precision_t result_v;
			int phi_idx = 1;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				result_u = integrand_u(1.0/3.0, 1.0/3.0, coefficients, coordinates, 
						phi_idx);

				result_v = integrand_v(1.0/3.0, 1.0/3.0, coefficients, coordinates, 
						phi_idx); 

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

		int integrate_analytical_resp(
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_vector* result_vector, precision_t cu_2, precision_t cv_2){
			/*
			 * Computes the integral of the proposed analytical solution on the
			 * element domain
			*/
			precision_t result_u;
			precision_t result_v;
			int phi_idx = 1;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				result_u = integrand_analytical_u(1.0/3.0, 1.0/3.0, coordinates, 
						phi_idx, cu_2);
				result_v = integrand_analytical_v(1.0/3.0, 1.0/3.0, coordinates, 
						phi_idx, cv_2);

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
			/*
			 * Updates result_matrix with the local contribution to the
			 * jacobian computed at coefficients, method used by the 
			 * nonliner solver
			*/
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
					result_uu = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx); 

					result_uv = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx + 3);

					result_vu = diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx);

					result_vv = diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx + 3);

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
	

		int update_sp_with_linearized_int_0(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_spmatrix* global_stiffness){
			/*
			 * Updates global_stiffness with the approximation of H obtained by
			 * dropping terms and generatin a linear H.
			*/
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			node_t n_3 = this->nodes_[2];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t r3 = coordinates[n_3][0];

			precision_t r_u1 = V_MU/K_MU;
			precision_t r_u2 = V_MU/K_MU;
			precision_t r_u3 = V_MU/K_MU;

			precision_t r_v1 = RESP_Q*V_MU/K_MU;
			precision_t r_v2 = RESP_Q*V_MU/K_MU;
			precision_t r_v3 = RESP_Q*V_MU/K_MU;

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
			
			return EXIT_SUCCESS;
		}
		
		int update_sp_with_linearized_int_2(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_spmatrix* global_stiffness){
			/*
			 * Updates global_stiffness with the linear approximation obtained
			 * by represention the H function in the phi basis
			*/
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
				global_vect_t* vector_f){
			/*
			 * Updates the f vector with the aproximation obtained by
			 * representing the nonlinear function in the phi basis
			*/
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
			
			gsl_vector_set(vector_f, n_1,
					gsl_vector_get(vector_f, n_1) 
					+ f_u_1);
			gsl_vector_set(vector_f, n_2,
					gsl_vector_get(vector_f, n_2) 
					+ f_u_2);
			gsl_vector_set(vector_f, n_3,
					gsl_vector_get(vector_f, n_3) 
					+ f_u_3);
			gsl_vector_set(vector_f, n_1 + NUM_NODES,
					gsl_vector_get(vector_f, n_1 + NUM_NODES) 
					- f_v_1);
			gsl_vector_set(vector_f, n_2 + NUM_NODES,
					gsl_vector_get(vector_f, n_2 + NUM_NODES) 
					- f_v_2);
			gsl_vector_set(vector_f, n_3 + NUM_NODES,
					gsl_vector_get(vector_f, n_3 + NUM_NODES) 
					- f_v_3);
			
			return EXIT_SUCCESS;	
		}

		int update_sp_with_jacobian(const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_spmatrix* result_matrix){
			/*
			 * Updates result_matrix with te contribution to the jacobian of the
			 * element
			*/
			precision_t result_uu = 0;
			precision_t result_uv = 0;
			precision_t result_vu = 0;
			precision_t result_vv = 0;
			size_t node_1;
			size_t node_2;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				for (int coeff_idx = 1; coeff_idx < 4; coeff_idx++){
				node_1 = this->nodes_[node_idx];	
				node_2 = this->nodes_[coeff_idx - 1];	
				result_uu = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients,
						coordinates, node_idx, coeff_idx); 

				result_uv = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients, 
						coordinates, node_idx, coeff_idx + 3);

				result_vu = diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
						coordinates, node_idx, coeff_idx);

				result_vv = diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
						coordinates, node_idx, coeff_idx + 3); 

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
				stiff_mat_t* global_stiffness){
			/*
			 * Updates the global_stiffness, using the local stiffness matrix
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
			precision_t k = (r1 + r2 + r3)/(6*jacobian(coordinates));
			precision_t s11 = (SIGMA_UR*pow(z3 - z2, 2) + 
					SIGMA_UZ*pow(r3 - r2, 2))*k;
			precision_t s12 = (SIGMA_UR*(z3 - z2)*(z1 - z3) +
					SIGMA_UZ*(r1 - r3)*(r3 - r2))*k;
			precision_t s13 = (SIGMA_UZ*(r1 - r2)*(r2 - r3) +
					SIGMA_UR*(z1 - z2)*(z2 - z3))*k;
			precision_t s22 = (SIGMA_UR*pow(z1 - z3, 2) + 
					SIGMA_UZ*pow(r1 - r3, 2))*k;
			precision_t s23 = -(SIGMA_UR*(z1 - z2)*(z1 - z3) +
					SIGMA_UZ*(r1 - r2)*(r1 - r3))*k;
			precision_t s33 = (SIGMA_UR*pow(z1 - z2, 2) + 
					SIGMA_UZ*pow(r1 - r2, 2))*k;
			
			gsl_spmatrix_set(global_stiffness, n_1, n_1, 
					gsl_spmatrix_get(global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2, 
					gsl_spmatrix_get(global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(global_stiffness, n_1, n_3, 
					gsl_spmatrix_get(global_stiffness, n_1, n_3) + s13);
			gsl_spmatrix_set(global_stiffness, n_2, n_1, 
					gsl_spmatrix_get(global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2, 
					gsl_spmatrix_get(global_stiffness, n_2, n_2) + s22);
			gsl_spmatrix_set(global_stiffness, n_2, n_3, 
					gsl_spmatrix_get(global_stiffness, n_2, n_3) + s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_1, 
					gsl_spmatrix_get(global_stiffness, n_3, n_1) + s13);
			gsl_spmatrix_set(global_stiffness, n_3, n_2, 
					gsl_spmatrix_get(global_stiffness, n_3, n_2) + s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_3, 
					gsl_spmatrix_get(global_stiffness, n_3, n_3) + s33);
			
			n_1 += NUM_NODES;
			n_2 += NUM_NODES;			
			n_3 += NUM_NODES;

			s11 = (SIGMA_VR*pow(z3 - z2, 2) + 
					SIGMA_VZ*pow(r3 - r2, 2))*k;
			s12 = (SIGMA_VR*(z3 - z2)*(z1 - z3) +
					SIGMA_VZ*(r1 - r3)*(r3 - r2))*k;
			s13 = (SIGMA_VZ*(r1 - r2)*(r2 - r3) +
					SIGMA_VR*(z1 - z2)*(z2 - z3))*k;
			s22 = (SIGMA_VR*pow(z1 - z3, 2) + 
					SIGMA_VZ*pow(r1 - r3, 2))*k;
			s23 = -(SIGMA_VR*(z1 - z2)*(z1 - z3) +
					SIGMA_VZ*(r1 - r2)*(r1 - r3))*k;
			s33 = (SIGMA_VR*pow(z1 - z2, 2) + 
					SIGMA_VZ*pow(r1 - r2, 2))*k;

			gsl_spmatrix_set(global_stiffness, n_1, n_1, 
					gsl_spmatrix_get(global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2, 
					gsl_spmatrix_get(global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(global_stiffness, n_1, n_3, 
					gsl_spmatrix_get(global_stiffness, n_1, n_3) + s13);
			gsl_spmatrix_set(global_stiffness, n_2, n_1, 
					gsl_spmatrix_get(global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2, 
					gsl_spmatrix_get(global_stiffness, n_2, n_2) + s22);
			gsl_spmatrix_set(global_stiffness, n_2, n_3, 
					gsl_spmatrix_get(global_stiffness, n_2, n_3) + s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_1, 
					gsl_spmatrix_get(global_stiffness, n_3, n_1) + s13);
			gsl_spmatrix_set(global_stiffness, n_3, n_2, 
					gsl_spmatrix_get(global_stiffness, n_3, n_2) + s23);
			gsl_spmatrix_set(global_stiffness, n_3, n_3, 
					gsl_spmatrix_get(global_stiffness, n_3, n_3) + s33);
			
			return EXIT_SUCCESS;
		}
		int update_stiffness_with_linearized_integral(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates,
				stiff_mat_t* global_stiffness){
			update_sp_with_jacobian(coefficients, 
					coordinates, &global_stiffness);
			return EXIT_SUCCESS;
		
		}
		
		int update_f_vector_with_linearized_integral(
				const gsl_vector* coefficients,
				const std::vector<std::vector<precision_t>>& coordinates, 
				global_vect_t* vector_f){
			/*
			 * Updates the f vector with the constant vector obtained by
			 * linearising the non linear function H using Taylor approximation
			*/
			integrate_non_linear_term(coefficients, coordinates, vector_f);
			
			precision_t result_uu = 0;
			precision_t result_uv = 0;
			precision_t result_vu = 0;
			precision_t result_vv = 0;
			size_t node_1;
			for (int node_idx = 0; node_idx < 3; node_idx++){
				for (int coeff_idx = 1; coeff_idx < 4; coeff_idx++){
					node_1 = this->nodes_[node_idx];	
					result_uu = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx);

					result_uv = diff_integrand_u(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx + 3);

					result_vu =	diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx);

					result_vv = diff_integrand_v(1.0/3.0, 1.0/3.0, coefficients, 
							coordinates, node_idx, coeff_idx + 3);

					gsl_vector_set(vector_f, node_1,
							gsl_vector_get(vector_f, node_1) 
							- (result_uu*C_U_AMB + result_uv*C_V_AMB));
					gsl_vector_set(vector_f, node_1 + NUM_NODES,
							gsl_vector_get(vector_f, node_1 + NUM_NODES) 
							+ (result_vu*C_U_AMB + result_vv*C_V_AMB));
					} 
			}
			return EXIT_SUCCESS;	
		}
};

template <typename P>
class ElementBoundary : public Element<P>{
	public:
		using typename Element<P>::precision_t;
		using typename Element<P>::node_t;
		using typename Element<P>::stiff_mat_t;
		using typename Element<P>::global_vect_t;
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

		inline precision_t length(
				const std::vector<std::vector<precision_t>>& coordinates){
			return sqrt(pow((coordinates[this->nodes_[0]][0] - 
						coordinates[this->nodes_[1]][0]), 2) +
				pow((coordinates[this->nodes_[0]][1] - 
							coordinates[this->nodes_[1]][1]), 2));
		}

		// Access
		const std::vector<node_t>& nodes(){
			return this->nodes_;
		}
		
		int update_stiffness_matrix(
				const std::vector<std::vector<precision_t>>& coordinates,
				stiff_mat_t* global_stiffness){
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

			gsl_spmatrix_set(global_stiffness, n_1, n_1, gsl_spmatrix_get(
						global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2, gsl_spmatrix_get(
						global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_1, gsl_spmatrix_get(
						global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2, gsl_spmatrix_get(
						global_stiffness, n_2, n_2) + s22);

			n_1 += NUM_NODES;
			n_2 += NUM_NODES;
			s11 = RHO_V*length(coordinates)*(r1/4  + r2/12);
			s12 = RHO_V*length(coordinates)*(r1 + r2)/12;
			s22 = RHO_V*length(coordinates)*(r2/4  + r1/12);

			gsl_spmatrix_set(global_stiffness, n_1, n_1, gsl_spmatrix_get(
						global_stiffness, n_1, n_1) + s11);
			gsl_spmatrix_set(global_stiffness, n_1, n_2, gsl_spmatrix_get(
						global_stiffness, n_1, n_2) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_1, gsl_spmatrix_get(
						global_stiffness, n_2, n_1) + s12);
			gsl_spmatrix_set(global_stiffness, n_2, n_2, gsl_spmatrix_get(
						global_stiffness, n_2, n_2) + s22);
			return EXIT_SUCCESS;
		}

		int update_vector_f(
				const std::vector<std::vector<precision_t>>& coordinates, 
				global_vect_t* vector_f){
			if (axis_flag){
				return EXIT_SUCCESS;
			}
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			gsl_vector_set(vector_f, n_1, gsl_vector_get(vector_f, n_1) - 
					RHO_U*length(coordinates)*C_U_AMB*(r1/3 + r2/6));
			gsl_vector_set(vector_f, n_2, gsl_vector_get(vector_f, n_2) - 
					RHO_U*length(coordinates)*C_U_AMB*(r2/3 + r1/6));

			n_1 += NUM_NODES;
			n_2 += NUM_NODES;
			gsl_vector_set(vector_f, n_1, gsl_vector_get(vector_f, n_1) - 
					RHO_V*length(coordinates)*C_V_AMB*(r1/3 + r2/6));
			gsl_vector_set(vector_f, n_2, gsl_vector_get(vector_f, n_2) - 
					RHO_V*length(coordinates)*C_V_AMB*(r2/3 + r1/6));
			return EXIT_SUCCESS;
		}

		int get_normal_dir_coords(
				const std::vector<std::vector<precision_t>>& coordinates,
				std::vector<precision_t>& res){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t z1 = coordinates[n_1][1];
			precision_t z2 = coordinates[n_2][1];
			precision_t slope = (z2 - z1)/(r2 - r1);
			res[0] = (r1 + r2)/2 + 0.5*length(coordinates)/
				sqrt((1 + pow(-1/slope),2));
			res[1] = (z1 + z2)/2 - 0.5*length(coordinates)/
				(slope*sqrt((1 + pow(-1/slope),2)));
			return EXIT_SUCCESS;
		}

		int get_midpoint_val(
				const std::vector<std::vector<precision_t>>& coordinates,
				gsl_vector* coefficients, std::vector<precision_t>& res){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			res[0] = (gsl_vector_get(coefficients, n_1) + 
					gsl_vector_get(coefficients, n_2))/2;
			res[1] = (gsl_vector_get(coefficients, n_1 + NUM_NODES) + 
					gsl_vector_get(coefficients, n_2 + NUM_NODES))/2;
			return EXIT_SUCCESS;
		}

		int get_grad(std::vector<precision_t>& coeff_1, 
				std::vector<precision_t>& coeff_2,
				const std::vector<std::vector<precision_t>>& coordinates, 
				std::vector<precision_t>& grad){
			node_t n_1 = this->nodes_[0];
			node_t n_2 = this->nodes_[1];
			precision_t r1 = coordinates[n_1][0];
			precision_t r2 = coordinates[n_2][0];
			precision_t z1 = coordinates[n_1][1];
			precision_t z2 = coordinates[n_2][1];
			precision_t slope = (z2 - z1)/(r2 - r1);
			grad[0] = (coeff_1[0] - coeff_2[0])/(0.5*length(coordinates)/
				sqrt((1 + pow(-1/slope),2)));
			grad[1] = (coeff_1[0] - coeff_2[0])/(-0.5*length(coordinates)/
				(slope*sqrt((1 + pow(-1/slope),2))));
			grad[2] = (coeff_1[1] - coeff_2[1])/(0.5*length(coordinates)/
				sqrt((1 + pow(-1/slope),2)));
			grad[3] = (coeff_1[1] - coeff_2[1])/(-0.5*length(coordinates)/
				(slope*sqrt((1 + pow(-1/slope),2))));
			return EXIT_SUCCESS;
		}
};
}
