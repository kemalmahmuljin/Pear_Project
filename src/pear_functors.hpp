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
			model_.get_integral_vector(quad_points_);
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

std::string matrix_to_string(const gsl_spmatrix* matrix_to_print){
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

int write_matrix_to_file(const gsl_spmatrix* mtrx, std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::out);
	myfile<<matrix_to_string(mtrx);
	myfile.close();	
}

int write_vector_to_file(const gsl_vector* vect, std::string filename){
	std::ofstream myfile;
	myfile.open(filename, std::ios::out);
	myfile<<vector_to_string(vect);
	myfile.close();	
}
}
