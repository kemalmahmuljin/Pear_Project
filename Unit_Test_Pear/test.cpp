#include "pch.h"
#include "../src/importer.hpp"


TEST(TestImporter, NodeCoordinates) {
  FEM_module::ImporterMsh<double, long> mesh_importer("../input/pear.msh");
  mesh_importer.process_file();
  const std::vector<std::vector<double>> nodes = mesh_importer.node_matrix();
  for (int i = 0; i < nodes.size(); ++i) {
      // Checking if each node has 3 inputs: x,y,z coordinates
      ASSERT_EQ(nodes[i].size(), 3) << "Vector does not have size 3 at index" << i;
  }
}

TEST(TestImporter, Boundary ) {

}