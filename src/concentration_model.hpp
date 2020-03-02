#include <gsl/gsl_multiroots.h>

#include <vector>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>

#ifndef IMPORTER_HPP
#define IMPORTER_HPP
#include "importer.hpp"
#endif

#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include "element.hpp"
#endif

#ifndef PEAR_FUNCTORS_HPP
#define PEAR_FUNCTORS_HPP
#include "pear_functors.hpp"
#endif

namespace FEM_module{

template<typename P, typename I>
class NonLinearSystemFunctor;

template<typename P, typename I>
class ConcentrationModel{
	public:
		typedef P precision_t;
		typedef I node_t;
		typedef std::vector<std::vector<precision_t>> coord_mtx_t;
		typedef std::vector<FEM_module::ElementTriangular<P, I>> elem_list_t;
		typedef std::vector<FEM_module::ElementBoundary<P, I>> bound_list_t;
	private:
		node_t number_elements_;
		node_t number_boundaries_;
		node_t number_nodes_;
		std::vector<std::vector<precision_t>> coordinates_;
		elem_list_t elements_;
		bound_list_t boundaries_;
		gsl_vector* f_vector_;
		gsl_vector* coefficients_;
		gsl_vector* helper_;
		gsl_spmatrix* stiffness_matrix_;
	public:
		/*
		 *
		 *
		 *
		 *
		 *
		 *
		 * Very important question: can I free the importer after using it, I've
		 * try that before, and when the importer runs out of scope the default
		 * destructor gets called and i get an error, but I don't ant to wait
		 * until the importer gets out of scope, I can free the memory before.
		 * Also, if I pass the coordinates as a reference from the importer to
		 * the ConcentrationModel, and then destroy the Importer, it's the
		 * matrix still reachable from the ConcentrationModel or would I end up
		 * with a dangling pointer, should I pass the matrices as copies?.
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		*/
		// Constructor
		ConcentrationModel(Importer<precision_t, node_t>& importer, 
				std::vector<precision_t>& interior_point)
		: number_elements_()
		, number_boundaries_()
		, number_nodes_()
		, elements_(elem_list_t())
		, boundaries_(bound_list_t())
		, coordinates_()
		, f_vector_()
		, helper_()
		, coefficients_()
		// once we have an estimate of the number of nonzero elements on the
		// stiffness matrix it would be better to allocate using 
		// gsl_spmatrix_alloc_nzmax()
		, stiffness_matrix_()
		{
			number_nodes_ = importer.node_num();
			number_elements_ = importer.element_num();
			number_boundaries_ = importer.boundary_num();
			stiffness_matrix_ = gsl_spmatrix_alloc((size_t)(2*number_nodes_), 
					(size_t)(2*number_nodes_));
			gsl_spmatrix_set_zero(stiffness_matrix_);
			f_vector_ = gsl_vector_calloc((size_t)(2*number_nodes_));
			helper_ = gsl_vector_calloc((size_t)(2*number_nodes_));
			coefficients_ = gsl_vector_calloc((size_t)(2*number_nodes_));
			coordinates_ = importer.node_matrix();
			for (auto elem : importer.element_matrix()){
				FEM_module::ElementTriangular<precision_t, node_t> 
					temp_element(coordinates_, elem);
				elements_.push_back(temp_element);
			}
			for (auto elem : importer.boundary_matrix()){
				FEM_module::ElementBoundary<precision_t, node_t> 
					temp_boundary(coordinates_, elem, interior_point);
				boundaries_.push_back(temp_boundary);
			}	
		}
		// Destructor
		~ConcentrationModel(){
			gsl_spmatrix_free(stiffness_matrix_);
			gsl_vector_free(f_vector_);
			gsl_vector_free(helper_);
			gsl_vector_free(coefficients_);
		}
		// Access
		const gsl_spmatrix* stiffness_matrix(){
			return stiffness_matrix_;
		}

		const gsl_vector* f_vector(){
			return f_vector_;
		}

		const gsl_vector* coefficients(){
			return coefficients_;
		}

		gsl_vector* helper(){
			return helper_;
		}

		node_t number_nodes(){
			return number_nodes_;
		}
		// Setters
		
		int set_coefficients(const gsl_vector* source){
			gsl_vector_memcpy(coefficients_, source);
			return EXIT_SUCCESS;
		}

		// Functionality
		int generate_stiffness_matrix(){
			for (auto elem : elements_){
				elem.update_stiffness_matrix(coordinates_, 
						*stiffness_matrix_);
			}
			for (auto bound : boundaries_){
				bound.update_stiffness_matrix(coordinates_, 
						*stiffness_matrix_);
			}
			return EXIT_SUCCESS;
		}
		
		int generate_f_vector(){
			for (auto bound : boundaries_){
				bound.update_vector_f(coordinates_, *f_vector_);
			}
			return EXIT_SUCCESS;
		}

		int solve_linear_model(){
			gsl_spmatrix* stiff_mat_cc;
			const precision_t tolerance = 1.0e-10;
			const size_t max_iter = 100000;
			const gsl_splinalg_itersolve_type* itersolve_type = 
				gsl_splinalg_itersolve_gmres;
			gsl_splinalg_itersolve *work = 
				gsl_splinalg_itersolve_alloc(itersolve_type, 
						2*number_nodes_, 0);
			size_t iter = 0;
			precision_t residual;
			int status;

			gsl_vector_set_zero(coefficients_);
			
			stiff_mat_cc = gsl_spmatrix_ccs(stiffness_matrix_);
			do{
				status = gsl_splinalg_itersolve_iterate(stiff_mat_cc, f_vector_,
						tolerance, coefficients_, work);
				std::cout<<coefficients()<<std::endl;
				residual = gsl_splinalg_itersolve_normr(work);
				if (status == GSL_SUCCESS)
					fprintf(stderr, "Converged\n");
			}
			while (status == GSL_CONTINUE && ++iter < max_iter);

			gsl_splinalg_itersolve_free(work);
			gsl_spmatrix_free(stiff_mat_cc);
			return EXIT_SUCCESS;
		}

		int	get_integral_vector(gsl_vector* integral_vect, 
				size_t cuad_points){
			gsl_vector_set_all(integral_vect, 0.0);
			for (auto elem : elements_){
				elem.integrate_non_linear_term(coefficients_, coordinates_,
						cuad_points, integral_vect);
			}
			return EXIT_SUCCESS;
		}

		int solve_nonlinear_model(){
			int condition = 0;
			const gsl_multiroot_fsolver_type* nonlinear_solver_type;
			nonlinear_solver_type = gsl_multiroot_fsolver_dnewton;
			gsl_multiroot_fsolver* nonlinear_solver = 
				gsl_multiroot_fsolver_alloc(nonlinear_solver_type,
						(size_t)(2*number_nodes_));
			// We need the initial guess for the solver, for now I'm just
			// passing 0
			gsl_vector_set_all(coefficients_, 0);
			
			FEM_module::NonLinearSystemFunctor<precision_t, node_t> 
				nls_functor(*this, 5);
			
			gsl_multiroot_function nls_function;
			nls_function.f = &FEM_module::non_linear_function;
			nls_function.n = 2*number_nodes_;
			nls_function.params = &nls_functor;

			gsl_multiroot_fsolver_set(nonlinear_solver, &nls_function, 
					coefficients_);
			do {
				gsl_multiroot_fsolver_iterate(nonlinear_solver);
				condition = gsl_multiroot_test_residual(
						gsl_multiroot_fsolver_f(nonlinear_solver), 1e-10);
			} while(condition != GSL_SUCCESS);
			gsl_multiroot_fsolver_free(nonlinear_solver);	
			return EXIT_SUCCESS;
		}

	// Contructor{
	// Create an Importer
	// Read the mesh
	// Initialize model using the Importer obeject
	// Create global stiffness matrix
	// Vector f
	// function handle for the non-linear part( taken as input)
	// }
	//
	// Functionality
	// 1_ compute initial guess{
	// creates a matrix 
	// create linear system
	// solve it
	// }
	//
	// 2_ solve(){
	// solve the system using guess and the non linear solver
	// }
};

std::ostream& operator<<(std::ostream& os, const gsl_spmatrix* sp_mat){
	for (size_t i = 0; i < sp_mat->size1; i++){
		for (size_t j = 0; j < sp_mat->size2; j++){
			os<<gsl_spmatrix_get(sp_mat, i, j)<<" ";
		}
		os<<std::endl;
	}
    return os;
}

std::ostream& operator<<(std::ostream& os, const gsl_vector* vect){
	for (size_t i = 0; i < vect->size; i++){
		os<<gsl_vector_get(vect, i)<<" ";
	}
	os<<std::endl;
    return os;
}

template <typename P, typename I>
std::ostream& operator<<(std::ostream& os, ConcentrationModel<P, I>& model){
    os<<"Stiffness Matrix"<<std::endl;
	os<<model.stiffness_matrix()<<std::endl;
    os<<"Vector f"<<std::endl;
	os<<model.f_vector()<<std::endl;
    os<<"Concentration Coefficients"<<std::endl;
	os<<model.coefficients()<<std::endl;
    return os;
}
}
