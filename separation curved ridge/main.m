%% calculate the large principle curvature



clear all; clc
delete(gcp('nocreate'))
p = parpool(16);

currentFolder = pwd;

sheets_choice = 1;
warning('off','all')

% generate the flow map
t = 3;
tspan = [0:0.1:t];

%% generate surface plot

% grid separation for the 3d meshgrid
grid_separation = 0.005;

% choose the grid point for(u,v)
u_list = [-2:grid_separation:2];
v_list = [-2:grid_separation:2];
% we only need u v to characterize a surface
% however, in order to calculate the del F, we need grid points in 3d, so
% we choose to label each parallel surface with h_list
h_list = [0.1,0.2,0.3];
h_list = [0.01,h_list];

% number of parallel surfaces
num_parallel_surface = length(h_list);

% notice that u_grid = actual x position = column in the matrix (from left to right)
% notice that v_grid = actual y position = row in the matrix (from top to bottom)
[u_2d_grid,v_2d_grid] = meshgrid(u_list,v_list);


[x1_0,x2_0,x3_0] = sigma_surface_plane(u_2d_grid,v_2d_grid,h_list);


% 3d grid needed to calculate del F
% first we need to calculate the flow fields for grid point in order to
% calculate del F
min_x1 = min(x1_0(:))-5*grid_separation;
max_x1 = max(x1_0(:))+5*grid_separation;
min_x2 = min(x2_0(:))-5*grid_separation;
max_x2 = max(x2_0(:))+5*grid_separation;
min_x3 = min(x3_0(:))-5*grid_separation;
max_x3 = max(x3_0(:))+5*grid_separation;

[F1_0,F2_0,F3_0] = meshgrid([min_x1:grid_separation:max_x1],...
    [min_x2:grid_separation:max_x2],[min_x3:grid_separation:max_x3]);


% plot the surface using surf
figure()
for i = 1:num_parallel_surface
    surf(x1_0(:,:,i),x2_0(:,:,i),x3_0(:,:,i));
    hold on;
end
% checking if the grid points cover the surface
%scatter3(F1_0(:),F2_0(:),F3_0(:))
xlabel('x')
ylabel('y')
zlabel('z')
shading interp
title('surface with radom color using surf');


%% calculate tangent vectors sigma_u and sigma_v

% tagent vector at t=0
% since we know the surface mapping, we can calculate the tangent vectors
% for the plane
sigma_u_1 = ones(size(x1_0));
sigma_u_2 = zeros(size(x1_0));
sigma_u_3 = zeros(size(x1_0));

sigma_v_1 = zeros(size(x1_0));
sigma_v_2 = ones(size(x1_0));
sigma_v_3 = zeros(size(x1_0));


%% second order derivative of sigma

% since we know the surface mapping

% sigma_uu
sigma_uu_1 = zeros(size(x1_0));
sigma_uu_2 = zeros(size(x1_0));
sigma_uu_3 = zeros(size(x1_0));

% sigma_uv
sigma_uv_1 = zeros(size(x1_0));
sigma_uv_2 = zeros(size(x1_0));
sigma_uv_3 = zeros(size(x1_0));

% sigma_vv
sigma_vv_1 = zeros(size(x1_0));
sigma_vv_2 = zeros(size(x1_0));
sigma_vv_3 = zeros(size(x1_0));


%% prepare for parallel computer in ode45

% assign points to multiple cores to speed up ode45 calculation

% number of points in the grid
num_points=length(F1_0(:));

% max number of cores you want to assign the task to
MaxNcores = 16;
cpu_num = min(MaxNcores,num_points);

% id of the boundary of points assigned to each core
id = ceil( linspace(0,num_points,cpu_num+1) );

% % flow map of the grid

%optionsT= odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep',1e-2);

% calculate the flow map
tic
[t_series, F1_tSeries_grid,F2_tSeries_grid,F3_tSeries_grid] = flow_map_curved_ridge(tspan,F1_0,F2_0,F3_0,cpu_num,id);
toc



%% prepare for calculate eigenvalues

small_eigenvalues_t = cell(length(t_series),1);

large_eigenvalues_t = cell(length(t_series),1);

gaussian_curvature_t = cell(length(t_series),1);

% eigenvector in tangent space
zeta_2_u = cell(length(t_series),1);
zeta_2_v = cell(length(t_series),1);

x1_all_t = cell(length(t_series),1);
x2_all_t = cell(length(t_series),1);
x3_all_t = cell(length(t_series),1);

for kk = 1:length(t_series)
    %% choose the last time point to analyze first
    
    % get the flow map at the last time point t = 0.5
    % F1_t = F1_tSeries_grid{end};
    % F2_t = F2_tSeries_grid{end};
    % F3_t = F3_tSeries_grid{end};
    
    F1_t = F1_tSeries_grid{kk};
    F2_t = F2_tSeries_grid{kk};
    F3_t = F3_tSeries_grid{kk};
    
    % get the flow map at the last time point for the surface
    % x1_t_ode = x1_tSeries_ode{end};
    % x2_t_ode = x2_tSeries_ode{end};
    % x3_t_ode = x3_tSeries_ode{end};
    
    %% interpolation of the flow map and flow map of the surface
    
    % From above, we have calculated the flow map on the grid points. Next, we
    % need to calculate the gridded interpolation to extract the flow map for
    % the surface
    
    
    % use interpolation to calculate the flow map of the surfaces
    x1_t_interp = interp3(F1_0,F2_0,F3_0,F1_t,x1_0,x2_0,x3_0);
    x2_t_interp = interp3(F1_0,F2_0,F3_0,F2_t,x1_0,x2_0,x3_0);
    x3_t_interp = interp3(F1_0,F2_0,F3_0,F3_t,x1_0,x2_0,x3_0);
    
    
    
    %% calculate the gradient of the flow map DelF
    
    [F1_x,F1_y,F1_z] = gradient(F1_t,grid_separation);
    [F2_x,F2_y,F2_z] = gradient(F2_t,grid_separation);
    [F3_x,F3_y,F3_z] = gradient(F3_t,grid_separation);
    
    % calculate the interpolated delF
    DelF1_x = interp3(F1_0,F2_0,F3_0,F1_x,x1_0,x2_0,x3_0);
    DelF1_y = interp3(F1_0,F2_0,F3_0,F1_y,x1_0,x2_0,x3_0);
    DelF1_z = interp3(F1_0,F2_0,F3_0,F1_z,x1_0,x2_0,x3_0);
    
    DelF2_x = interp3(F1_0,F2_0,F3_0,F2_x,x1_0,x2_0,x3_0);
    DelF2_y = interp3(F1_0,F2_0,F3_0,F2_y,x1_0,x2_0,x3_0);
    DelF2_z = interp3(F1_0,F2_0,F3_0,F2_z,x1_0,x2_0,x3_0);
    
    DelF3_x = interp3(F1_0,F2_0,F3_0,F3_x,x1_0,x2_0,x3_0);
    DelF3_y = interp3(F1_0,F2_0,F3_0,F3_y,x1_0,x2_0,x3_0);
    DelF3_z = interp3(F1_0,F2_0,F3_0,F3_z,x1_0,x2_0,x3_0);
    
    
    
    %% tangent vector at Fu Fv
    
    Fu_1 = DelF1_x.*sigma_u_1 + DelF1_y.*sigma_u_2 + DelF1_z.*sigma_u_3;
    Fu_2 = DelF2_x.*sigma_u_1 + DelF2_y.*sigma_u_2 + DelF2_z.*sigma_u_3;
    Fu_3 = DelF3_x.*sigma_u_1 + DelF3_y.*sigma_u_2 + DelF3_z.*sigma_u_3;
    
    Fv_1 = DelF1_x.*sigma_v_1 + DelF1_y.*sigma_v_2 + DelF1_z.*sigma_v_3;
    Fv_2 = DelF2_x.*sigma_v_1 + DelF2_y.*sigma_v_2 + DelF2_z.*sigma_v_3;
    Fv_3 = DelF3_x.*sigma_v_1 + DelF3_y.*sigma_v_2 + DelF3_z.*sigma_v_3;
    
    
    
    %% calculate the first fundamental form
    
    E = Fu_1.^2 + Fu_2.^2 + Fu_3.^2;
    F = Fu_1.*Fv_1+Fu_2.*Fv_2+Fu_3.*Fv_3;
    G = Fv_1.^2 + Fv_2.^2 + Fv_3.^2;
    J = sqrt( E.*G-F.^2 );
    
    
    
    %% calculate the normal vector
    
    N1 = -(Fu_2.*Fv_3 - Fu_3.*Fv_2)./J;
    N2 = +(Fu_1.*Fv_3 - Fu_3.*Fv_1)./J;
    N3 = -(Fu_1.*Fv_2 - Fu_2.*Fv_1)./J;
    
    check_normal_unit = N1.^2 + N2.^2 + N3.^3;
    
    
    %% calculate the Del^2F
    
    [F1_xx,F1_xy,F1_xz] = gradient(F1_x,grid_separation);
    [F1_yx,F1_yy,F1_yz] = gradient(F1_y,grid_separation);
    [F1_zx,F1_zy,F1_zz] = gradient(F1_z,grid_separation);
    
    [F2_xx,F2_xy,F2_xz] = gradient(F2_x,grid_separation);
    [F2_yx,F2_yy,F2_yz] = gradient(F2_y,grid_separation);
    [F2_zx,F2_zy,F2_zz] = gradient(F2_z,grid_separation);
    
    [F3_xx,F3_xy,F3_xz] = gradient(F3_x,grid_separation);
    [F3_yx,F3_yy,F3_yz] = gradient(F3_y,grid_separation);
    [F3_zx,F3_zy,F3_zz] = gradient(F3_z,grid_separation);
    
    % Del^2 F_ij1
    F1_xx_interp = interp3(F1_0,F2_0,F3_0,F1_xx,x1_0,x2_0,x3_0);
    F1_yx_interp = interp3(F1_0,F2_0,F3_0,F1_yx,x1_0,x2_0,x3_0);
    F1_zx_interp = interp3(F1_0,F2_0,F3_0,F1_zx,x1_0,x2_0,x3_0);
    
    F2_xx_interp = interp3(F1_0,F2_0,F3_0,F2_xx,x1_0,x2_0,x3_0);
    F2_yx_interp = interp3(F1_0,F2_0,F3_0,F2_yx,x1_0,x2_0,x3_0);
    F2_zx_interp = interp3(F1_0,F2_0,F3_0,F2_zx,x1_0,x2_0,x3_0);
    
    F3_xx_interp = interp3(F1_0,F2_0,F3_0,F3_xx,x1_0,x2_0,x3_0);
    F3_yx_interp = interp3(F1_0,F2_0,F3_0,F3_yx,x1_0,x2_0,x3_0);
    F3_zx_interp = interp3(F1_0,F2_0,F3_0,F3_zx,x1_0,x2_0,x3_0);
    
    % Del^2 F_ij2
    F1_xy_interp = interp3(F1_0,F2_0,F3_0,F1_xy,x1_0,x2_0,x3_0);
    F1_yy_interp = interp3(F1_0,F2_0,F3_0,F1_yy,x1_0,x2_0,x3_0);
    F1_zy_interp = interp3(F1_0,F2_0,F3_0,F1_zy,x1_0,x2_0,x3_0);
    
    F2_xy_interp = interp3(F1_0,F2_0,F3_0,F2_xy,x1_0,x2_0,x3_0);
    F2_yy_interp = interp3(F1_0,F2_0,F3_0,F2_yy,x1_0,x2_0,x3_0);
    F2_zy_interp = interp3(F1_0,F2_0,F3_0,F2_zy,x1_0,x2_0,x3_0);
    
    F3_xy_interp = interp3(F1_0,F2_0,F3_0,F3_xy,x1_0,x2_0,x3_0);
    F3_yy_interp = interp3(F1_0,F2_0,F3_0,F3_yy,x1_0,x2_0,x3_0);
    F3_zy_interp = interp3(F1_0,F2_0,F3_0,F3_zy,x1_0,x2_0,x3_0);
    
    % Del^2 F_ij3
    
    F1_xz_interp = interp3(F1_0,F2_0,F3_0,F1_xz,x1_0,x2_0,x3_0);
    F1_yz_interp = interp3(F1_0,F2_0,F3_0,F1_yz,x1_0,x2_0,x3_0);
    F1_zz_interp = interp3(F1_0,F2_0,F3_0,F1_zz,x1_0,x2_0,x3_0);
    
    F2_xz_interp = interp3(F1_0,F2_0,F3_0,F2_xz,x1_0,x2_0,x3_0);
    F2_yz_interp = interp3(F1_0,F2_0,F3_0,F2_yz,x1_0,x2_0,x3_0);
    F2_zz_interp = interp3(F1_0,F2_0,F3_0,F2_zz,x1_0,x2_0,x3_0);
    
    F3_xz_interp = interp3(F1_0,F2_0,F3_0,F3_xz,x1_0,x2_0,x3_0);
    F3_yz_interp = interp3(F1_0,F2_0,F3_0,F3_yz,x1_0,x2_0,x3_0);
    F3_zz_interp = interp3(F1_0,F2_0,F3_0,F3_zz,x1_0,x2_0,x3_0);
    
    
    %% calculate second derivative in F
    
    % del2 sigma_u
    d2su_11 = sigma_u_1.*F1_xx_interp + sigma_u_2.*F1_xy_interp + sigma_u_3.*F1_xz_interp;
    d2su_12 = sigma_u_1.*F1_yx_interp + sigma_u_2.*F1_yy_interp + sigma_u_3.*F1_yz_interp;
    d2su_13 = sigma_u_1.*F1_zx_interp + sigma_u_2.*F1_zy_interp + sigma_u_3.*F1_zz_interp;
    
    d2su_21 = sigma_u_1.*F2_xx_interp + sigma_u_2.*F2_xy_interp + sigma_u_3.*F2_xz_interp;
    d2su_22 = sigma_u_1.*F2_yx_interp + sigma_u_2.*F2_yy_interp + sigma_u_3.*F2_yz_interp;
    d2su_23 = sigma_u_1.*F2_zx_interp + sigma_u_2.*F2_zy_interp + sigma_u_3.*F2_zz_interp;
    
    d2su_31 = sigma_u_1.*F3_xx_interp + sigma_u_2.*F3_xy_interp + sigma_u_3.*F3_xz_interp;
    d2su_32 = sigma_u_1.*F3_yx_interp + sigma_u_2.*F3_yy_interp + sigma_u_3.*F3_yz_interp;
    d2su_33 = sigma_u_1.*F3_zx_interp + sigma_u_2.*F3_zy_interp + sigma_u_3.*F3_zz_interp;
    
    % del2 sigma_v
    d2sv_11 = sigma_v_1.*F1_xx_interp + sigma_v_2.*F1_xy_interp + sigma_v_3.*F1_xz_interp;
    d2sv_12 = sigma_v_1.*F1_yx_interp + sigma_v_2.*F1_yy_interp + sigma_v_3.*F1_yz_interp;
    d2sv_13 = sigma_v_1.*F1_zx_interp + sigma_v_2.*F1_zy_interp + sigma_v_3.*F1_zz_interp;
    
    d2sv_21 = sigma_v_1.*F2_xx_interp + sigma_v_2.*F2_xy_interp + sigma_v_3.*F2_xz_interp;
    d2sv_22 = sigma_v_1.*F2_yx_interp + sigma_v_2.*F2_yy_interp + sigma_v_3.*F2_yz_interp;
    d2sv_23 = sigma_v_1.*F2_zx_interp + sigma_v_2.*F2_zy_interp + sigma_v_3.*F2_zz_interp;
    
    d2sv_31 = sigma_v_1.*F3_xx_interp + sigma_v_2.*F3_xy_interp + sigma_v_3.*F3_xz_interp;
    d2sv_32 = sigma_v_1.*F3_yx_interp + sigma_v_2.*F3_yy_interp + sigma_v_3.*F3_yz_interp;
    d2sv_33 = sigma_v_1.*F3_zx_interp + sigma_v_2.*F3_zy_interp + sigma_v_3.*F3_zz_interp;
    
    % Fuu
    Fuu_1 = sigma_u_1.*d2su_11 + sigma_u_2.*d2su_12 + sigma_u_3.*d2su_13;
    Fuu_1 = Fuu_1 + sigma_uu_1.*DelF1_x + sigma_uu_2.*DelF1_y + sigma_uu_3.*DelF1_z;
    
    Fuu_2 = sigma_u_1.*d2su_21 + sigma_u_2.*d2su_22 + sigma_u_3.*d2su_23;
    Fuu_2 = Fuu_2 + sigma_uu_1.*DelF2_x + sigma_uu_2.*DelF2_y + sigma_uu_3.*DelF2_z;
    
    Fuu_3 = sigma_u_1.*d2su_31 + sigma_u_2.*d2su_32 + sigma_u_3.*d2su_33;
    Fuu_3 = Fuu_3 + sigma_uu_1.*DelF3_x + sigma_uu_2.*DelF3_y + sigma_uu_3.*DelF3_z;
    
    % Fuv
    Fuv_1 = sigma_v_1.*d2su_11 + sigma_v_2.*d2su_12 + sigma_v_3.*d2su_13;
    Fuv_1 = Fuv_1 + sigma_uv_1.*DelF1_x + sigma_uv_2.*DelF1_y + sigma_uv_3.*DelF1_z;
    
    Fuv_2 = sigma_v_1.*d2su_21 + sigma_v_2.*d2su_22 + sigma_v_3.*d2su_23;
    Fuv_2 = Fuv_2 + sigma_uv_1.*DelF2_x + sigma_uv_2.*DelF2_y + sigma_uv_3.*DelF2_z;
    
    Fuv_3 = sigma_v_1.*d2su_31 + sigma_v_2.*d2su_32 + sigma_v_3.*d2su_33;
    Fuv_3 = Fuv_3 + sigma_uv_1.*DelF3_x + sigma_uv_2.*DelF3_y + sigma_uv_3.*DelF3_z;
    
    % Fvv
    Fvv_1 = sigma_v_1.*d2sv_11 + sigma_v_2.*d2sv_12 + sigma_v_3.*d2sv_13;
    Fvv_1 = Fvv_1 + sigma_vv_1.*DelF1_x + sigma_vv_2.*DelF1_y + sigma_vv_3.*DelF1_z;
    
    Fvv_2 = sigma_v_1.*d2sv_21 + sigma_v_2.*d2sv_22 + sigma_v_3.*d2sv_23;
    Fvv_2 = Fvv_2 + sigma_vv_1.*DelF2_x + sigma_vv_2.*DelF2_y + sigma_vv_3.*DelF2_z;
    
    Fvv_3 = sigma_v_1.*d2sv_31 + sigma_v_2.*d2sv_32 + sigma_v_3.*d2sv_33;
    Fvv_3 = Fvv_3 + sigma_vv_1.*DelF3_x + sigma_vv_2.*DelF3_y + sigma_vv_3.*DelF3_z;
    
    
    
    %% calculate the second fundamental form
    
    L = N1.*Fuu_1 + N2.*Fuu_2 + N3.*Fuu_3;
    
    M = N1.*Fuv_1 + N2.*Fuv_2 + N3.*Fuv_3;
    
    N = N1.*Fvv_1 + N2.*Fvv_2 + N3.*Fvv_3;
    
    
    det_F2 = L.*N-M.^2;
    
    
    %% calculate the weigarden map
    
    % deformation square
    J_2 = E.*G-F.^2;
    
    W11 = (L.*G-M.*F)./J_2;
    W12 = (M.*G-N.*F)./J_2;
    W21 = (M.*E-L.*F)./J_2;
    W22 = (N.*E-M.*F)./J_2;
    
    % turn W into matrix form for each u,v
    W_uv = cell(size(x1_0));
    
    % large eigenvalue
    large_eig = zeros(size(x1_0));
    % small eigenvalue
    small_eig = zeros(size(x1_0));
    
    % eigenvector corresponding to large eigenvalue in the tangent vector basis
    eig_vec_large_sigmaU = zeros(size(x1_0));
    eig_vec_large_sigmaV = zeros(size(x1_0));
    
    % eigenvector corresponding to small eigenvalue in the tangent vector basis
    eig_vec_small_sigmaU = zeros(size(x1_0));
    eig_vec_small_sigmaV = zeros(size(x1_0));
    
    length_u = length(u_list);
    length_v = length(v_list);
    
    parfor i = 1:length_v
        for j = 1:length_u
            for k = 1:num_parallel_surface
                
                W_uv{i,j,k} = [W11(i,j,k),W12(i,j,k);W21(i,j,k),W22(i,j,k)];
                
                if sum(isnan(W_uv{i,j,k}(:))) == 0
                    [eigen_vec_matrix,eig_uv_matrix] = eig(W_uv{i,j,k});
                    eig_uv_list = diag(eig_uv_matrix);
                    eig_uv_list_sort = sort((eig_uv_list));
                    large_eig(i,j,k) = eig_uv_list_sort(2);
                    
                    small_eig(i,j,k) = eig_uv_list_sort(1);
                    
                    
                    if eig_uv_list_sort(2) == (eig_uv_list(2))
                        % eigenvector corresponding to large eigenvalue
                        eig_vec_large_sigmaU(i,j,k) = eigen_vec_matrix(1,2);
                        eig_vec_large_sigmaV(i,j,k) = eigen_vec_matrix(2,2);
                        % eigenvector corresponding to small eigenvalue
                        eig_vec_small_sigmaU(i,j,k) = eigen_vec_matrix(1,1);
                        eig_vec_small_sigmaV(i,j,k) = eigen_vec_matrix(2,1);
                        
                    else
                        % eigenvector corresponding to large eigenvalue
                        eig_vec_large_sigmaU(i,j,k) = eigen_vec_matrix(1,1);
                        eig_vec_large_sigmaV(i,j,k) = eigen_vec_matrix(2,1);
                        % eigenvector corresponding to small eigenvalue
                        eig_vec_small_sigmaU(i,j,k) = eigen_vec_matrix(1,2);
                        eig_vec_small_sigmaV(i,j,k) = eigen_vec_matrix(2,2);
                    end
                    
                else
                    disp(strcat('i = ',num2str(i),',j = ',num2str(j)));
                    
                end
            end
        end
    end
    
    
    
    
    % eigenvector corresponding to large eigenvalue
    eig_vec_large_1 = eig_vec_large_sigmaU.*Fu_1+ eig_vec_large_sigmaV.*Fv_1;
    eig_vec_large_2 = eig_vec_large_sigmaU.*Fu_2+ eig_vec_large_sigmaV.*Fv_2;
    eig_vec_large_3 = eig_vec_large_sigmaU.*Fu_3+ eig_vec_large_sigmaV.*Fv_3;
    % eigenvector corresponding to small eigenvalue
    eig_vec_small_1 = eig_vec_small_sigmaU.*Fu_1+ eig_vec_small_sigmaV.*Fv_1;
    eig_vec_small_2 = eig_vec_small_sigmaU.*Fu_2+ eig_vec_small_sigmaV.*Fv_2;
    eig_vec_small_3 = eig_vec_small_sigmaU.*Fu_3+ eig_vec_small_sigmaV.*Fv_3;
    
    % eigenvector in the tangent space
    zeta_2_u{kk} = eig_vec_large_sigmaU;
    zeta_2_v{kk} = eig_vec_large_sigmaV;
    %% save the eigenvalues
    
    x1_all_t{kk} = x1_t_interp;
    x2_all_t{kk} = x2_t_interp;
    x3_all_t{kk} = x3_t_interp;
    
    small_eigenvalues_t{kk} = small_eig;
    large_eigenvalues_t{kk} = large_eig;
    %===============================================================================================
    % calculate the gaussian curvature
    
    gaussian_curvature_t{kk} = large_eig.*small_eig;
end



%% backbone from equation
tic

disp('extracting the backbone......')

% store the selection of points used as backbone
logical_condition_surf_sepa = cell(num_parallel_surface,1);

% choose the last time point to extract the backbone
kk = length(t_series);


for i = 1:num_parallel_surface
    
    % eigenvalues on each surface at the end time points
    large_eig = large_eigenvalues_t{kk}(:,:,i);
    small_eig = small_eigenvalues_t{kk}(:,:,i);
    
    % eigenvector zeta_2 in tangent space
    eig_vec_large_sigmaU = zeta_2_u{kk}(:,:,i);
    eig_vec_large_sigmaV = zeta_2_v{kk}(:,:,i);
    
    % gradient of the large eigenvalue
    [k2_u,k2_v] = gradient(large_eig,grid_separation);
    
    % grad^2 of the large eigenvalue
    [k2_uu,k2_uv] = gradient(k2_u,grid_separation);
    [k2_vu,k2_vv] = gradient(k2_v,grid_separation);
    
    % first order directional derivative
    D_zeta = k2_u.*eig_vec_large_sigmaU + k2_v.*eig_vec_large_sigmaV;
    
    % get the contours of zero
    cc = contour(x1_0(:,:,i),x2_0(:,:,i),D_zeta,[0,0]);
    
    % get the positions of the contour
    ss = getcontourlines(cc);
    x_D = ss(2).x;
    y_D = ss(2).y;
    
    % round to get the logical array
    x_D_round = round(x_D,2);
    y_D_round = round(y_D,2);
    x1_0_round = round(x1_0(:,:,i),2);
    x2_0_round = round(x2_0(:,:,i),2);
    
    [logical_xy_arr,~] = ismember([x1_0_round(:),x2_0_round(:)],[x_D_round', y_D_round'],'rows');
    
    condition_2 = reshape(logical_xy_arr,size(x1_0_round));
    
    % second order directional derivative
    D2_zeta = k2_uu.*eig_vec_large_sigmaU.^2 + k2_vv.*eig_vec_large_sigmaV.^2 ...
        + 2*k2_uv.*eig_vec_large_sigmaU.*eig_vec_large_sigmaV;
    
    % condition for the backbone
    condition_1 = (large_eig>0);
    
    condition_3 = (D2_zeta<0);
    overall_condition = logical(condition_1.*condition_2.*condition_3);
    
    logical_condition_surf_sepa{i} = overall_condition;
    
end

%%
figure()
for i = 1:num_parallel_surface
    subplot(1,4,i)
    x = x1_0(:,:,i);
    y = x2_0(:,:,i);
    scatter(x(logical_condition_surf_sepa{i}),y(logical_condition_surf_sepa{i}),2,'filled')
end

%% advect the backbone

sepa_t_x = cell(length(t_series),num_parallel_surface);
sepa_t_y = cell(length(t_series),num_parallel_surface);
sepa_t_z = cell(length(t_series),num_parallel_surface);

parfor ct = 1:length(t_series)
    x1_t_interp = x1_all_t{ct};
    x2_t_interp = x2_all_t{ct};
    x3_t_interp = x3_all_t{ct};
    
    
    for i = 1: num_parallel_surface
        
        x_backbone = [];
        y_backbone = [];
        z_backbone = [];
        
        overall_condition_curr_surf = logical_condition_surf_sepa{i};
        
        x_curr_surf = x1_t_interp(:,:,i);
        y_curr_surf = x2_t_interp(:,:,i);
        z_curr_surf = x3_t_interp(:,:,i);
        
        x_backbone = [x_backbone;x_curr_surf(overall_condition_curr_surf)];
        y_backbone = [y_backbone;y_curr_surf(overall_condition_curr_surf)];
        z_backbone = [z_backbone;z_curr_surf(overall_condition_curr_surf)];
        
        sepa_t_x{ct,i} = x_backbone;
        sepa_t_y{ct,i} = y_backbone;
        sepa_t_z{ct,i} = z_backbone;
    end
    
end

%% find the unique points

round_digit = 2;

sepa_t_x_unique_round = cell(size(sepa_t_x));
sepa_t_y_unique_round = cell(size(sepa_t_y));
sepa_t_z_unique_round = cell(size(sepa_t_z));

parfor ct = 1:length(t_series)
    for i = 1: num_parallel_surface
        x_t_i = sepa_t_x{ct,i};
        y_t_i = sepa_t_y{ct,i};
        z_t_i = sepa_t_z{ct,i};
        
        % round the number to find unique values
        x_t_i_round = round(x_t_i,round_digit);
        y_t_i_round = round(y_t_i,round_digit);
        z_t_i_round = round(z_t_i,round_digit);
        
        xyz_t_i_round = [x_t_i_round,y_t_i_round,z_t_i_round];
        
        
        xyz_t_i_unique = unique(xyz_t_i_round,'rows');
        
        sepa_t_x_unique_round{ct,i} = xyz_t_i_unique(:,1);
        sepa_t_y_unique_round{ct,i} = xyz_t_i_unique(:,2);
        sepa_t_z_unique_round{ct,i} = xyz_t_i_unique(:,3);
    end
end

%% connecting the backbone

disp('just backbone......')

folder_just_backbone = strcat(currentFolder,'\backbone',num2str(t));

if ~exist(folder_just_backbone, 'dir')
    mkdir(folder_just_backbone)
end

% choose first and last time point to plot
plot_t_idx_arr = [1,length(t_series)];

f1 = figure('units','normalized','outerposition',[0 0 1 1]);

for plot_idx = 1:2
    
    subplot(1,2,plot_idx)
    
    % current time point
    ct = plot_t_idx_arr(plot_idx);
    
    % flow map of the surface
    x1_t_interp = x1_all_t{ct};
    x2_t_interp = x2_all_t{ct};
    x3_t_interp = x3_all_t{ct};
    
    % real time
    t = t_series(ct);
    
    % eigenvalues at the current surface at time t
    large_eig = large_eigenvalues_t{ct};
    
    % plot the backbone and large principle curvature
    hold on
    for i = 1:sheets_choice:num_parallel_surface
       if i ~= 1
            surf( x1_t_interp(:,:,i),x2_t_interp(:,:,i),x3_t_interp(:,:,i),large_eig(:,:,i),'FaceAlpha',0.5);
            %scatter3(sepa_t_x{ct,i},sepa_t_y{ct,i},sepa_t_z{ct,i},2,'filled','r');
       end
    end
    grid off
    shading interp
    
    scatter3(sepa_t_x{ct,1},sepa_t_y{ct,1},sepa_t_z{ct,1},10,'filled','k');
    
    %===============================================================================================
    
    % connecting the backbone to surface
    
    
    % number of points in each surface
    numPts_surfs = zeros(1,num_parallel_surface);
    count = 0;
    for surf_idx = 1:num_parallel_surface
        count = count+length(sepa_t_x_unique_round{ct,surf_idx});
        numPts_surfs(surf_idx) = count;
    end
    numPts_surfs = [0,numPts_surfs];
    % backbone from all surfaces
    x_spline = [];
    y_spline = [];
    z_spline = [];
    for i = 1: num_parallel_surface
        x_spline = [x_spline;sepa_t_x_unique_round{ct,i}];
        y_spline = [y_spline;sepa_t_y_unique_round{ct,i}];
        z_spline = [z_spline;sepa_t_z_unique_round{ct,i}];
    end
    
    % get ride of triangles made of points entirely on the same layer
    tri = delaunay(y_spline,z_spline);
    
    % check which edge contains elements all in the same layer
    logical_all = ones(size(tri,1),1);
    
    for surf_idx = 1 : num_parallel_surface
        
        % idx of pts in the current layer
        idx_curr = [numPts_surfs(surf_idx)+1:numPts_surfs(surf_idx+1)];
        
        for line = 1:size(tri,1)
            
            if sum(ismember(tri(line,:),idx_curr))==3
                logical_all(line) = 0;
            end
            
        end
        
    end
    
    tri_new = tri(logical(logical_all),:);
    
    h = trisurf(tri_new, x_spline, y_spline,z_spline,'FaceColor','r','LineStyle','none','FaceAlpha',0.5);
    %===============================================================================================
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %zlim([0 max(x3_all_t{end}(:))])
    title(strcat('(a)','$\gamma =k_2$ and separation backbone at t = ',num2str(t)),'Interpreter','latex')
    colorbar
    box on
    view(30,40)
     
    
end
set(findall(gcf,'-property','FontSize'),'FontSize',24) 
savename='lambda_2_analysis_2plot';
%saveas(f1,strcat(folder_just_backbone,'\',savename,'.png'));


%% connecting the backbone


disp('4 plot movie with backbone......')

folder_4_plot_backbone = strcat(currentFolder,'\4 plot with backbone',num2str(t));

if ~exist(folder_4_plot_backbone, 'dir')
    mkdir(folder_4_plot_backbone)
end

for ct = 1:length(t_series)
    
    % flow map of the surface
    x1_t_interp = x1_all_t{ct};
    x2_t_interp = x2_all_t{ct};
    x3_t_interp = x3_all_t{ct};
    
    % real time
    t = t_series(ct);
    
    % eigenvalues at the current surface at time t
    large_eig = large_eigenvalues_t{ct};
    small_eig = small_eigenvalues_t{ct};
    
    gaussian_curvature = gaussian_curvature_t{ct};
    
    f1=figure('visible','off','units','normalized','outerposition',[0 0 1 1]);
    
    % plot the surface at time t
    subplot(2,2,1)
    hold on;
    for i = 1:sheets_choice:num_parallel_surface
        if i~=1
            surf(x1_t_interp(:,:,i),x2_t_interp(:,:,i),x3_t_interp(:,:,i),'FaceAlpha',0.7);
        end
    end
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %zlim([0 max(x3_all_t{end}(:))])
    title(strcat('Advected surface at t = ',num2str(t)),'Interpreter','latex')
    shading interp
    colorbar
    box on
    view(30,40)
    
    % plot the large principle curvature
    subplot(2,2,2)
    mean_curvature = (large_eig + small_eig)/2;
    hold on
    for i = 1:sheets_choice:num_parallel_surface
        if i~=1
            surf( x1_t_interp(:,:,i),x2_t_interp(:,:,i),x3_t_interp(:,:,i),mean_curvature(:,:,i),'FaceAlpha',0.7 );
        end
    end
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %zlim([0 max(x3_all_t{end}(:))])
    title(strcat('Mean curvature on the surface at t = ',num2str(t)),'Interpreter','latex')
    shading interp
    colorbar
    box on
    view(30,40)
    
    % plot the gaussian curvature
    subplot(2,2,3)
    hold on
    for i = 1:sheets_choice:num_parallel_surface
        if i ~=1
            surf( x1_t_interp(:,:,i),x2_t_interp(:,:,i),x3_t_interp(:,:,i),gaussian_curvature(:,:,i) ,'FaceAlpha',0.7);
        end
    end
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %zlim([0 max(x3_all_t{end}(:))])
    title(strcat('gaussian curvature on the surface at t = ',num2str(t)),'Interpreter','latex')
    shading interp
    colorbar
    box on
    view(30,40)
    %caxis([-1 100])
    
    % plot the separation line and large principle curvature
    % plot the separation line and large principle curvature
    subplot(2,2,4)
    hold on
    for i = 1:sheets_choice:num_parallel_surface
        if i~=1
            surf( x1_t_interp(:,:,i),x2_t_interp(:,:,i),x3_t_interp(:,:,i),large_eig(:,:,i),'FaceAlpha',0.7);
            %scatter3(sepa_t_x{ct,i},sepa_t_y{ct,i},sepa_t_z{ct,i},2,'filled','r');
        end
    end
    grid off
    shading interp
    scatter3(sepa_t_x{ct,1},sepa_t_y{ct,1},sepa_t_z{ct,1},10,'filled','k');
    %===============================================================================================
    
    % connecting the backbone to surface
    
    % backbone at the top surface
    x_spline_top = sepa_t_x{ct,num_parallel_surface};
    y_spline_top = sepa_t_y{ct,num_parallel_surface};
    z_spline_top = sepa_t_z{ct,num_parallel_surface};
    
    % backbone from all surfaces
    x_spline = [];
    y_spline = [];
    z_spline = [];
    for i = 1: num_parallel_surface
        x_spline = [x_spline;sepa_t_x{ct,i}];
        y_spline = [y_spline;sepa_t_y{ct,i}];
        z_spline = [z_spline;sepa_t_z{ct,i}];
    end
    
    % get ride of triangles made of points entirely on the top layer
    tri = delaunay(y_spline,z_spline);
    
    % idx of points on the top layer
    idx = [1:length(z_spline)];
    idx_top = idx((end-length(z_spline_top)+1):end);
    
    logical_arr = ones(size(tri,1),1);
    
    for line = 1:size(tri,1)
        
        if sum(ismember(tri(line,:),idx_top))==3
            logical_arr(line) = 0;
        end
        
    end
    
    tri_new = tri(logical(logical_arr),:);
    
    h = trisurf(tri_new, x_spline, y_spline,z_spline,'FaceColor','r','LineStyle','none');
    %===============================================================================================
    xlabel('x')
    ylabel('y')
    zlabel('z')
    %zlim([0 max(x3_all_t{end}(:))])
    title(strcat('$\gamma =k_2$ and separation line at t = ',num2str(t)),'Interpreter','latex')
    colorbar
    box on
    view(30,40)
    
    set(findall(gcf,'-property','FontSize'),'FontSize',24) 
    savename=strcat('lambda_2_analysis_t =',num2str(ct));
    saveas(f1,strcat(folder_4_plot_backbone,'\',savename,'.tif'));
    
    
    
end





%%

save('separation_curved_ridge_v2.mat','-v7.3')
delete(p)