#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "Importer.hpp"


using namespace std;

namespace FEM_module{

template <typename P, typename I>
class ElementTriangular{
	public:
		typedef P precision_t;
		typedef I node_t;
	protected:
		vector<node_t> nodes_;
	public:
		ElementTriangular(const vector<vector<precision_t>>& coordinates, 
				const vector<node_t>& nodes)
		: nodes_(vector<node_t>())
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
		}

		int get_global_coords(const vector<precision_t>& coordinates, 
				vector<precision_t>& output_coords){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 2; j++){
					output_coords[i][j] = coordinates[i][j];
				}
			}
			return EXIT_SUCCESS;
		}
		
		// Access
		const vector<node_t>& nodes(){
			return nodes_;
		}

	private:
		bool get_rotation(const vector<vector<precision_t>>& coordinates, 
				const vector<node_t>& nodes){
		/* Function use to get the orientation of the nodes, in order to order
		them always counter-clockwise.
		Output:
			True: Counter-Clockwise order.
			False: Clockwise order.
		*/
		precision_t slope_1 = 
			(coordinates[nodes[1]][2] - coordinates[nodes[0]][2])/
			(coordinates[nodes[1]][1] - coordinates[nodes[0]][1]);
		precision_t slope_2 = 
			(coordinates[nodes[2]][2] - coordinates[nodes[1]][2])/
			(coordinates[nodes[2]][1] - coordinates[nodes[1]][1]);
		return (slope_1 < slope_2);
		}
};

template <typename P, typename I>
class ElementBoundary : public ElementTriangular<P, I>{
	public:
		using typename ElementTriangular<P, I>::precision_t;
		using typename ElementTriangular<P, I>::node_t;
	private:
		vector<node_t> boundary_nodes;
	public:
		ElementBoundary(
				const vector<vector<precision_t>>& coordinates, 
				const vector<node_t>& nodes)
		: ElementTriangular<precision_t, node_t>(coordinates, nodes)
		{}
};
}

int test2(){
	FEM_module::ImporterMsh<double, long> mesh_importer("pear.msh");
	mesh_importer.process_file();
	const vector<vector<double>> coords = mesh_importer.node_matrix();
	vector<long> nodes_from_element;
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
	test2();
	return EXIT_SUCCESS;
}
