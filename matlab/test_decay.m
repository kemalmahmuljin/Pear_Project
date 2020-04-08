%% Importing everything
boundaries = dlmread('../output/boundaries', ' ', 1, 0);
boundaries = boundaries(:, 1:3);
elements = dlmread('../output/elements', ' ', 1, 0);
elements = elements(:, 1:3);
coords = dlmread('../output/coords', ' ', 1, 0);
coords = coords(:, 1:2);

stiffness = dlmread('../output/stiff', ' ', 1, 0);
stiffness = stiffness(:, 1:end - 1);
stiffness_lin = dlmread('../output/stiff_lin', ' ', 1, 0);
stiffness_lin = stiffness_lin(:, 1:end - 1);
fvector = dlmread('../output/f_vector', ' ', 1, 0);
fvector = fvector(:, 1:end - 1)';
fvector_lin = dlmread('../output/f_vector_lin', ' ', 1, 0);
fvector_lin = fvector_lin(:, 1:end - 1)';

x0 = dlmread('../output/initial_coeff', ' ', 1, 0);
x0= x0(:, 1:end - 1)';

%% function decaying with inverse square
%{
result_coeff = ones(2*size(coords,1),1)';

for i=1:size(coords,1)
    r = coords(i,1);
    z = coords(i,2) - 0.04;
    result_coeff(i) = (0.001 + r^2 + z^2)^(-1); % 0.001 to not go to inf
    result_coeff(i + size(coords,1)) = (0.001 + r^2 + z^2)^(-1); % 0.001 to not go to inf
    %disp(result_coeff(i));
end
f_calc = stiffness*result_coeff';
f_calc_lin = stiffness_lin*result_coeff';

dlmwrite('output/f_calc',4,'delimiter',' ','precision',12,'-append');
dlmwrite('output/f_calc',f_calc','delimiter',' ','precision',12,'-append');
dlmwrite('output/f_calc_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('output/f_calc_lin',f_calc_lin','delimiter',' ','precision',12,'-append');
dlmwrite('output/calculated_coeff',4,'delimiter',' ','precision',12,'-append');
dlmwrite('output/calculated_coeff',result_coeff,'delimiter',' ','precision',12,'-append');
%}


%% function linearly increasing C in r
%%{
result_coeff = ones(2*size(coords,1),1)';

for i=1:size(coords,1)
    r = coords(i,1);
    result_coeff(i) = 150*r;
    result_coeff(i + size(coords,1)) = 0.2*r;
end
%dlmwrite('../output/calculated_coeff',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/calculated_coeff',result_coeff,'delimiter',' ','precision',12);
%}
