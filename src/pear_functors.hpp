#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_integration.h>

#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include "element.hpp"
#endif

namespace FEM_module{

template<typename P, typename I>
class ConcentrationModel;

template<typename P, typename I>
class NonLinearSystemFunctor{
	public:
		typedef P precision_t;
		typedef I node_t;
	private:
		FEM_module::ConcentrationModel<precision_t, node_t>& model_;
		size_t quad_points_;
	public:
		NonLinearSystemFunctor(FEM_module::ConcentrationModel<precision_t, 
				node_t>& model, size_t quad_points)
		: model_(model)
		, quad_points_(quad_points)
		{}
		
		int operator()(const gsl_vector* x, void* params,
				gsl_vector* f){
			model_.set_coefficients(x);
			gsl_vector_set_all(f, 0.0);
			gsl_spblas_dgemv(CblasNoTrans, 1.0, model_.stiffness_matrix(), 
					x, 0.0, f);
			gsl_vector_add(f, model_.f_vector());
			model_.get_integral_vector(model_.helper(), quad_points_);
			gsl_vector_add(f, model_.helper());
			return EXIT_SUCCESS;
		}
};

extern "C" int non_linear_function(const gsl_vector* x, void *param, 
		gsl_vector* f){
	NonLinearSystemFunctor<double, int> *my_functor =
		(NonLinearSystemFunctor<double, int> *)param;
	return (*my_functor)(x, NULL, f);
}

std::string vector_to_string(const gsl_vector* vector_to_print){
	std::stringstream str_stream;
	str_stream<<"("<<vector_to_print->size<<")"<<"[";
	for (size_t i = 0; i < vector_to_print->size; i++){
		str_stream<<gsl_vector_get(vector_to_print, i)<<" ";
	}
	str_stream<<"]"<<std::endl;
	return str_stream.str();
}

//template<typename P, typename I>
//class FirstIntegrandUFunctor{
//	public:
//		typedef P precision_t;
//		typedef I node_t;
//	private:
//		FEM_module::ElementTriangular<precision_t, node_t>& element_;
//		precision_t u_;
//		const gsl_vector* coefficients_;
//		const std::vector<std::vector<precision_t>>& coordinates_;
//		node_t node_idx_;
//	public:
//		IntegrandUFunctor(FEM_module::ElementTriangular<precision_t, 
//				node_t>& element, const gsl_vector* coefficients,
//				const std::vector<std::vector<precision_t>>& coordinates,
//				node_t node_idx)
//		: element_(element)
//		, u_()
//		, coefficients_(coefficients)
//		, coordinates_(coordinates)
//		, node_idx_(node_idx)
//		{}
//		
//		int set_u(precision_t u){
//			u_ = u;
//		}
//
//		precision_t operator()(double x, void* params){
//			return *element_.integrand_u(u_, x, coefficients_, coordinates_, 
//					node_idx_);
//		}
//};
//
//template<typename P, typename I>
//class SecondIntegrandUFunctor{
//	public:
//		typedef P precision_t;
//		typedef I node_t;
//	private:
//		FirstIntegrandUFunctor<precision_t, node_t> first_integrand_functor_;
//		size_t points_;
//	public:
//		IntegrandUFunctor(FEM_module::ElementTriangular<precision_t, 
//				node_t>& element, const gsl_vector* coefficients,
//				const std::vector<std::vector<precision_t>>& coordinates,
//				node_t node_idx, size_t points)
//		: first_integrand_functor(element, coefficients, coordinates, 
//				node_idx)
//		, points_(points)
//		{}
//
//		precision_t operator()(double x, void* params){
//			precision_t u;
//			precision_t w;
//			gsl_integration_glfixed_table* table = 
//				gsl_integration_glfixed_table_alloc(n);
//			gsl
//			first_integrand_functor_.set_u(  ...  )
//			return gsl_integration_glfixed(first_integrand_functor_,
//					0, 1, table);
//		}
//};
}
