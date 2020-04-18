clear all; close all;
%% Creating mesh using Matlab
% Import the surface which we mesh (pear)
model = createpde(1);
importGeometry(model,'../Input/fine_pear_matlab.stl');
% Mesh the pear with the matlab mesher
generateMesh(model, 'GeometricOrder', 'linear','Hmax', 2.5) % hmax determines max triangle length.

B = boundary(model.Mesh.Nodes', 1);
%% To show that B only contains boundary nodes
% figure;
% pdeplot(model)
% hold on
% for i=1:size(B,1)
%     plot(model.Mesh.Nodes(1,B(i)), model.Mesh.Nodes(2,B(i)), 'y*', 'MarkerSize', 15)
%     hold on
% end

%% Write information to code compatible files
boundariess = zeros(size(B,1) - 1,3);

for i= 1:(size(B,1) - 1)
   boundariess(i,1) = B(i + 1);
   boundariess(i,2) = B(i);
   if( ~(model.Mesh.Nodes(1,B(i))==0 && model.Mesh.Nodes(1,B(i + 1))==0) ) % splitting by type of boundary
       boundariess(i,3) = 1;
   end
end

%% Write text files in output folder
delete('../output/M_boundaries');
dlmwrite('../output/M_boundaries',boundariess -1,'delimiter',' ','precision',12,'-append');
delete('../output/M_coords');
dlmwrite('../output/M_coords',model.Mesh.Nodes'/1000,'delimiter',' ','precision',12,'-append');
delete('../output/M_elements');
dlmwrite('../output/M_elements',( model.Mesh.Elements - 1)','delimiter',' ','precision',12,'-append');
