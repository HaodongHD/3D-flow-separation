% the velcocity field
function  f123 = velocity_field_curved_ridge(t,x123_arr)


% the velocity field has to transform a vector to a vector, which means
% that x123_arr is a array that contains x1,x2,x3.

% Here the input x123_arr must be a column vector

% Here we use the velocity field in the saddle type seperation with b=2,
% and the rest of the parameters = 1

% number of points used to evaluate velocity field
num_points = length(x123_arr)/3;

% extract the x1 x2 x3 positions
x1 = x123_arr(1:num_points);
x2 = x123_arr(num_points+1:2*num_points);
x3 = x123_arr(2*num_points+1:3*num_points);

% calculate the component of the velocity field
f1 = 0*x1;
f2 = 0*x2;
epsilon = 0.2;
%f3 = 0.1./(x1.^2 +epsilon)+0.1*(0.05*x2.^2+1);
f3 = 0.1*(0.1*x2.^2+1).*x3.^2./(x1.^2 +epsilon);
f123 = [f1;f2;f3];

end


