clear all; close all;
%% Creating mesh using Matlab
% Import the surface which we mesh (pear)
model = createpde(1);
importGeometry(model,'../Input/fine_pear_matlab.stl');
% Mesh the pear with the matlab mesher
generateMesh(model, 'GeometricOrder', 'linear','Hmax', 5) % hmax determines max triangle length.


%% Create matrix once to be able to transform it at the end
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

tK = symrcm(K);

%% Transforming Nodes
tempE = zeros(size(n,1),size(n,2));
for i=1:size(n,2)
   tempE(1,i) = find(tK == n(1,i)); 
   tempE(2,i) = find(tK == n(2,i)); 
   tempE(3,i) = find(tK == n(3,i)); 
end
s1 = s(1,tK);
s2 = s(2,tK);
Nodes = [s1;s2];

%% Determine boundaries (splitting coefficient REALLY IMPORTANT)
B = boundary(Nodes', 1); %  parameter needs to be changed together with HMAX to get all the boundaries

%% Write information to code compatible files
boundariess = ones(size(B,1) - 1,3); % to separate between 2 boundary types
for i= 1:(size(B,1) - 1)
   boundariess(i,1) = B(i + 1);
   boundariess(i,2) = B(i);
   if( ~(Nodes(1,B(i))==0 && Nodes(1,B(i + 1))==0) ) % splitting by type of boundary
       boundariess(i,3) = 0;
   end
end

%% To show that B only contains boundary nodes
figure;
pdeplot(model)
hold on
for i=1:(size(B,1) - 1)
    plot(Nodes(1,B(i)), Nodes(2,B(i)), 'y*', 'MarkerSize', 15)
    hold on
end
hold off;
%% Plot to check if sparsity is achieved
%%{
K = zeros(size(s,2), size(s,2));
for e=1:size(tempE,2)
    r1 = tempE(1,e);
    r2 = tempE(2,e);
    r3 = tempE(3,e);
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
figure();
spy(K)
%}
%% Write text files in output folder
boundariess(:,1:2) = boundariess(:,1:2) - 1;
delete('../output/M_boundaries');
dlmwrite('../output/M_boundaries',boundariess,'delimiter',' ','precision',12,'-append');
delete('../output/M_coords');
dlmwrite('../output/M_coords',Nodes'/1000,'delimiter',' ','precision',12,'-append');
delete('../output/M_elements');
dlmwrite('../output/M_elements',( tempE - 1)','delimiter',' ','precision',12,'-append');
