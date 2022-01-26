function [x,y,z] = sigma_surface_plane(u,v,h)

num_parallel_surface = length(h);

x = zeros([size(u),num_parallel_surface]);
y = zeros([size(u),num_parallel_surface]);
z = zeros([size(u),num_parallel_surface]);

% coordinates of the parallel surface in 3D
for i = 1:num_parallel_surface
    
    x(:,:,i) = u;
    
    y(:,:,i) = v;
    
    z(:,:,i) = h(i)*ones(size(u));
end

end