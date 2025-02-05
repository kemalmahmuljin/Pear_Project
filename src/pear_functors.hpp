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

template<typename P>
class ConcentrationModel;

template<typename P>
class NonLinearSystemFunctor{
	public:
		typedef P precision_t;
		typedef size_t node_t;
	private:
		FEM_module::ConcentrationModel<precision_t>& model_;
	public:
		NonLinearSystemFunctor(
				FEM_module::ConcentrationModel<precision_t>& model)
		: model_(model)
		{}
		
		int operator()(const gsl_vector* x, void* params,
				gsl_vector* f){
			model_.set_coefficients(x);
			gsl_vector_set_all(f, 0.0);
			gsl_spblas_dgemv(CblasNoTrans, 1.0, model_.stiffness_matrix(), 
					x, 0.0, f);
			gsl_vector_add(f, model_.f_vector());
			model_.get_integral_vector();
			gsl_vector_add(f, model_.helper());
			return EXIT_SUCCESS;
		}
};

template<typename P>
class JacobianFunctor{
	public:
		typedef P precision_t;
		typedef size_t node_t;
	private:
		FEM_module::ConcentrationModel<precision_t>& model_;
	public:
		JacobianFunctor(FEM_module::ConcentrationModel<precision_t>& model)
		: model_(model)
		{}
		
		int operator()(const gsl_vector* x, void* params,
				gsl_matrix* j){
			model_.set_coefficients(x);
			gsl_spmatrix_sp2d(j, model_.stiffness_matrix());
			model_.update_matrix_with_jacobian(j);
			return EXIT_SUCCESS;
		}
};

template<typename P>
class FiniteDifferenceFunctor{
	public:
		typedef P precision_t;
		typedef size_t node_t;
	private:
		FEM_module::ConcentrationModel<precision_t>& model_;
		FEM_module::NonLinearSystemFunctor<precision_t>& funct_;
		precision_t epsilon_;
		gsl_vector* delta_;
		gsl_vector* function_val_;
		gsl_vector* function_val_delta_;
		gsl_matrix* stiff_;
	public:
		FiniteDifferenceFunctor(
				FEM_module::ConcentrationModel<precision_t>& model, 
				FEM_module::NonLinearSystemFunctor<precision_t>& funct, 
				precision_t epsilon)
		: model_(model)
		, funct_(funct)
		, epsilon_(epsilon){
			delta_ = gsl_vector_alloc(2*model_.number_nodes());
			function_val_ = gsl_vector_alloc(2*model_.number_nodes());
			function_val_delta_ = gsl_vector_alloc(2*model_.number_nodes());
			stiff_ = gsl_matrix_alloc(2*model_.number_nodes(), 
					2*model_.number_nodes());
			gsl_spmatrix_sp2d(stiff_, model_.stiffness_matrix());
		}
		~FiniteDifferenceFunctor(){
			gsl_vector_free(delta_);
			gsl_vector_free(function_val_);
			gsl_vector_free(function_val_delta_);
			gsl_matrix_free(stiff_);
		}

		int operator()(const gsl_vector* x, void* params,
				gsl_matrix* j){
			model_.set_coefficients(x);
			gsl_vector_memcpy(delta_, x);
			funct_(x, NULL, function_val_);
			for (size_t idx = 0; idx < j->size1; idx++){
				gsl_vector_set(delta_, idx, gsl_vector_get(delta_, idx) + 
						epsilon_);
				funct_(delta_, NULL, function_val_delta_);
				gsl_vector_sub(function_val_delta_, function_val_);
				gsl_vector_scale(function_val_delta_, 1.0/epsilon_);
				gsl_matrix_set_col(j, idx, function_val_delta_);
				gsl_vector_set(delta_, idx, gsl_vector_get(delta_, idx) - 
						epsilon_);
			}
			gsl_matrix_add(j, stiff_);
			return EXIT_SUCCESS;
		}

		int change_epsilon(precision_t new_epsilon){
			epsilon_ = new_epsilon;
			return EXIT_SUCCESS;	
		}
};

struct solver_params {
	FEM_module::NonLinearSystemFunctor<double>* func;
	FEM_module::JacobianFunctor<double>* jac;
};

extern "C" int non_linear_function(const gsl_vector* x, void *param, 
		gsl_vector* f){
	struct FEM_module::solver_params* parameters = 
		(struct FEM_module::solver_params*)param;
	NonLinearSystemFunctor<double> *my_functor =
		(NonLinearSystemFunctor<double> *)parameters->func;
	return (*my_functor)(x, NULL, f);
}

extern "C" int jacobian_function(const gsl_vector* x, void *param, 
		gsl_matrix* j){
	struct FEM_module::solver_params* parameters = 
		(struct FEM_module::solver_params*)param;
	JacobianFunctor<double> *my_functor =
		(JacobianFunctor<double>*)parameters->jac;
	return (*my_functor)(x, NULL, j);
}

extern "C" int fdf_function(const gsl_vector* x, void *param,
		gsl_vector* f, gsl_matrix* j){
	struct FEM_module::solver_params* parameters = 
		(struct FEM_module::solver_params*)param;
	NonLinearSystemFunctor<double>* my_functor1 =
		(NonLinearSystemFunctor<double>*)parameters->func;
	JacobianFunctor<double>* my_functor2 =
		(JacobianFunctor<double> *)parameters->jac;
	(*my_functor1)(x, NULL, f);
	(*my_functor2)(x, NULL, j);
	return EXIT_SUCCESS;
}

std::string vector_to_string(const gsl_vector* vector_to_print){
	std::stringstream str_stream;
	str_stream<<"("<<vector_to_print->size<<")"<<std::endl;
	for (size_t i = 0; i < vector_to_print->size; i++){
		str_stream<<gsl_vector_get(vector_to_print, i)<<" ";
	}
	str_stream<<std::endl;
	return str_stream.str();
}

std::string spmatrix_to_string(const gsl_spmatrix* matrix_to_print){
	std::stringstream str_stream;
	str_stream<<"("<<matrix_to_print->size1<<", "<<matrix_to_print->size2<<
		")"<<std::endl;
	for (size_t i = 0; i < matrix_to_print->size1; i++){
		for (size_t j = 0; j < matrix_to_print->size2; j++){
			str_stream<<gsl_spmatrix_get(matrix_to_print, i, j)<<" ";
		}
		str_stream<<std::endl;
	} 
	str_stream<<std::endl;
	return str_stream.str();
}

std::string matrix_to_string(const gsl_matrix* matrix_to_print){
	std::stringstream str_stream;
	str_stream<<"("<<matrix_to_print->size1<<", "<<matrix_to_print->size2<<
		")"<<std::endl;
	for (size_t i = 0; i < matrix_to_print->size1; i++){
		for (size_t j = 0; j < matrix_to_print->size2; j++){
			str_stream<<gsl_matrix_get(matrix_to_print, i, j)<<" ";
		}
		str_stream<<std::endl;
	} 
	str_stream<<std::endl;
	return str_stream.str();
}

int write_matrix_to_file(const gsl_spmatrix* mtrx, std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::out);
	myfile<<spmatrix_to_string(mtrx);
	myfile.close();	
	return EXIT_SUCCESS;
}

int write_den_matrix_to_file(const gsl_matrix* mtrx, std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::out);
	myfile<<matrix_to_string(mtrx);
	myfile.close();	
	return EXIT_SUCCESS;
}

int write_vector_to_file(const gsl_vector* vect, std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::out);
	myfile<<vector_to_string(vect);
	myfile.close();	
	return EXIT_SUCCESS;
}

double norm_1(gsl_vector* vect){
	double result = 0;
	for (size_t idx = 0; idx < vect->size; idx++){
		result += abs(gsl_vector_get(vect, idx));
	}
	return result;
}

double norm_2(gsl_vector* vect){
	double result = 0;
	for (size_t idx = 0; idx < vect->size; idx++){
		result += pow(gsl_vector_get(vect, idx), 2);
	}
	return result;
}

}
