clear all; close all;
%% Creating mesh using Matlab
% Import the surface which we mesh (pear)
model = createpde(1);
importGeometry(model,'../Input/fine_pear_matlab.stl');
% Mesh the pear with the matlab mesher
generateMesh(model, 'GeometricOrder', 'linear','Hmax', 2.5) % hmax determines max triangle length.

B = boundary(model.Mesh.Nodes', 0.75); %  parameter needs to be changed together with HMAX to get all the boundaries
%% To show that B only contains boundary nodes
% figure;
% pdeplot(model)
% hold on
% for i=1:size(B,1)
%     plot(model.Mesh.Nodes(1,B(i)), model.Mesh.Nodes(2,B(i)), 'y*', 'MarkerSize', 15)
%     hold on
% end

%% Write information to code compatible files
boundariess = ones(size(B,1) - 1,3);

for i= 1:(size(B,1) - 1)
   boundariess(i,1) = B(i + 1);
   boundariess(i,2) = B(i);
   if( ~(model.Mesh.Nodes(1,B(i))==0 && model.Mesh.Nodes(1,B(i + 1))==0) ) % splitting by type of boundary
       boundariess(i,3) = 0;
   end
end


%% Create matrix once
n = model.Mesh.Elements;
s = model.Mesh.Nodes;

K = zeros(size(s,2), size(s,2));
for e=1:size(n,2)
    r1 = n(1,e);
    r2 = n(2,e);
    r3 = n(3,e);
    K(r1,r1) = 1;
    K(r1,r2) = 1;
    K(r1,r3) = 1;
    K(r2,r1) = 1;
    K(r2,r2) = 1;
    K(r2,r3) = 1;
    K(r3,r1) = 1;
    K(r3,r2) = 1;
    K(r3,r3) = 1;
end
% spy(K)
tK = symrcm(K);
% Kn = K(tK,tK);
% spy(Kn);


%% All nodes need to be reordered
% start with boundaries
tempB = zeros(size(boundariess,1), size(boundariess,2));
for i=1:size(boundariess,1)
    r1 = boundariess(i,1);
    r2 = boundariess(i,2);
    tempB(i,1) = tK(r1);
    tempB(i,2) = tK(r2);
end
boundariess(:,1) = tempB(:,1);
boundariess(:,2) = tempB(:,2);
boundariess(:,1:2) = boundariess(:,1:2) - 1;

% reorder coordinates
Nodes = zeros(size(s,1),size(s,2));
for i=1:size(s,2)
    check = tK(i);
    Nodes(:,check) = s(:,i);
end

% check if correctly reodered
% figure;
% pdeplot(model)
% hold on
% for i=1:size(boundariess,1)
%     plot(Nodes(1,boundariess(i,1)), Nodes(2,boundariess(i,1)), 'y*', 'MarkerSize', 15)
%     hold on
% end

% reordering elements
tempE = zeros(size(n,1),size(n,2));
for i=1:size(n,2)
   r1 = n(1,i);
   r2 = n(2,i);
   r3 = n(3,i);
   c1 = tK(r1);
   c2 = tK(r2);
   c3 = tK(r3);
   tempE(1,i) = tK(r1);
   tempE(2,i) = tK(r2);
   tempE(3,i) = tK(r3);
end
%% Write text files in output folder
delete('../output/M_boundaries');
dlmwrite('../output/M_boundaries',boundariess,'delimiter',' ','precision',12,'-append');
delete('../output/M_coords');
dlmwrite('../output/M_coords',Nodes'/1000,'delimiter',' ','precision',12,'-append');
delete('../output/M_elements');
dlmwrite('../output/M_elements',( tempE - 1)','delimiter',' ','precision',12,'-append');
