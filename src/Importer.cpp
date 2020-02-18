#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <iostream>
#include <assert.h>


using namespace std;

namespace FEM_module{

std::vector<std::string> split(const std::string& input, char delimiter);
template <typename P, typename I>
class Importer{
	public:
		typedef P precision_t;
		typedef I mesh_size_t;
		typedef std::vector<vector<precision_t>> coord_vector_t;
		typedef std::vector<vector<mesh_size_t>> element_vector_t;
	protected:
		coord_vector_t node_matrix_;
		element_vector_t element_matrix_;
		element_vector_t boundary_matrix_;
		std::string file_path_;
		bool initialized_flag_;
	public:
		Importer()
		: node_matrix_(coord_vector_t())
		, element_matrix_(element_vector_t())
		, boundary_matrix_(element_vector_t())
		, initialized_flag_(false)
		{}
		Importer(std::string file_path)
		: file_path_(file_path)
		, node_matrix_(coord_vector_t())
		, element_matrix_(element_vector_t())
		, boundary_matrix_(element_vector_t())
		, initialized_flag_(false)
		{}
		
		// Interface declaration	
		
		// Access
		virtual const coord_vector_t node_matrix() = 0;
		virtual const element_vector_t element_matrix() = 0;
		virtual const element_vector_t boundary_matrix() = 0;
		
		// Functionality
		virtual int process_file() = 0;
};

template <typename P, typename I>
class ImporterMsh : public Importer<P, I>{
	public:
		using typename Importer<P, I>::precision_t;
		using typename Importer<P, I>::mesh_size_t;
		using typename Importer<P, I>::coord_vector_t;
		using typename Importer<P, I>::element_vector_t;
	public:
		ImporterMsh()
		: Importer<P, I>()
		{}
		ImporterMsh(std::string file_path)
		: Importer<P, I>(file_path)
		{}

		const coord_vector_t node_matrix(){
			return this->node_matrix_;
		}
		const element_vector_t element_matrix(){
			return this->element_matrix_;
		}
		const element_vector_t boundary_matrix(){
			return this->boundary_matrix_;
		}

		int get_node_matrix(ifstream& file_stream, std::string& line){
			mesh_size_t node_num = 0;
			mesh_size_t node_blocks = 0;
			mesh_size_t nodes_in_block = 0;
			mesh_size_t skiped_lines = 0;
			precision_t help_val;


			stringstream str_to_num;
			std::string data;
			std::vector<std::string> line_data;

			std::getline(file_stream, line);
			line_data = split(line, ' ');
			str_to_num << line_data[0];
			str_to_num >> node_blocks;
			str_to_num.clear();
			str_to_num <<line_data[1];
			str_to_num >> node_num;
			str_to_num.clear();
			
			// Reads each block of nodes, skips the node tags, because there are
			// ordered
			for (mesh_size_t block = 1; block <= node_blocks; block++){
				std::getline(file_stream, line);
				line_data = split(line, ' ');
				nodes_in_block = stol(line_data[3]);
				// skipping node tags
				skiped_lines = 0;
				while (skiped_lines < nodes_in_block){
					file_stream.ignore(std::numeric_limits<streamsize>::max(), 
							'\n');
					skiped_lines++;
				}
				for (mesh_size_t i = 0; i < nodes_in_block; i++){
					vector<precision_t> point_data;
					std::getline(file_stream, line);
					line_data = split(line, ' ');
					for (int i = 0; i < 3; i++){
						str_to_num << line_data[i];
						str_to_num >> help_val;
						str_to_num.clear();
						point_data.push_back(help_val);
					} 
					this->node_matrix_.push_back(point_data);
				}
			}
			std::getline(file_stream, line);
			assert(line == "$EndNodes");
			return EXIT_SUCCESS;
		}

		int get_elements(ifstream& file_stream, std::string& line){
			mesh_size_t elem_num = 0;
			mesh_size_t elem_blocks = 0;
			mesh_size_t elems_in_block = 0;
			
			mesh_size_t help_integer = 1;
			int elems_type = 0;
			std::string data;
			std::vector<std::string> line_data;
			stringstream str_to_num;

			std::getline(file_stream, line);
			assert(line == "$Elements");
			std::getline(file_stream, line);
			line_data = split(line, ' ');
			
			str_to_num << line_data[0];
			str_to_num >> elem_blocks;
			str_to_num.clear();
			str_to_num <<line_data[1];
			str_to_num >> elem_num;
			str_to_num.clear();
			
			// Reads each block of nodes, skips the node tags, because there are
			// ordered
			for (mesh_size_t block = 1; block <= elem_blocks; block++){
				std::getline(file_stream, line);
				line_data = split(line, ' ');
				elems_in_block = stol(line_data[3]);
				elems_type = stoi(line_data[2]);
				for (mesh_size_t i = 0; i < elems_in_block; i++){
					vector<mesh_size_t> element_data;
					std::getline(file_stream, line);
					line_data = split(line, ' ');
					if (elems_type == 1){
						for (int i = 1; i < 3; i++){

							element_data.push_back(stol(line_data[i]));
						} 
						this->boundary_matrix_.push_back(element_data);
					}
					else if(elems_type == 2){
						for (int i = 0; i < 3; i++){
							str_to_num << line_data[i];
							str_to_num >> help_integer;
							str_to_num.clear();
							element_data.push_back(help_integer);
						} 
						this->element_matrix_.push_back(element_data);
					}
				}
			}
			std::getline(file_stream, line);
			assert(line == "$EndElements");
			return EXIT_SUCCESS;
		}

		int process_file(){
			/* Function used to extract node, boundary, and elements data from a
			 * file
			*/
			std::string line;
			ifstream file_stream;
			
			file_stream.open(this->file_path_);
			std::getline(file_stream, line);
			// Starts rocessing nodes, gets the number of nodes.
			while (!(line == "$Nodes")){
				std::getline(file_stream, line);
			}

			get_node_matrix(file_stream, line);
			get_elements(file_stream, line);
			this->initialized_flag_ = true;
			return EXIT_SUCCESS;
		}
};

std::vector<std::string> split(const std::string& input, char delimiter){
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(input);
	while (std::getline(tokenStream, token, delimiter)){
		tokens.push_back(token);
	}
	return tokens;
}

template <typename P, typename I>
std::ostream& operator<<(std::ostream& os, ImporterMsh<P, I> importer) {
    os<<"Node Num - x - y"<<std::endl;
	I count = 1;
    for (auto elem : importer.node_matrix()) {
		os<<"Node "<<count<<" - "<<elem[0]<<" - "<<elem[1]<<std::endl;
		count++;
    }
    os<<std::endl<<"Element Num - node 1 - node 2 - node 3"<<std::endl;
	count = 1;
    for (auto elem : importer.element_matrix()) {
		os<<"Element "<<count<<" - "<<elem[0]<<" - "<<elem[1]<<" - "<<
			elem[2]<<std::endl;
		count++;
    }
    os<<std::endl<<"Boundary Num - node 1 - node 2"<<std::endl;
	count = 1;
    for (auto elem : importer.boundary_matrix()) {
		os<<"Boundary "<<count<<" - "<<elem[0]<<" - "<<elem[1]<<std::endl;
		count++;
    }
    return os;
}

}

int test1(){
	/* Test to check the correct reading of a file */
	FEM_module::ImporterMsh<double, long> mesh_importer("../Input/pear.msh");
	mesh_importer.process_file();
	std::cout<<mesh_importer<<std::endl;
	return 0;
}

int main(){
	test1();
}
