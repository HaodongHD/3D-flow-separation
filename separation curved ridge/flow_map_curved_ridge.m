% using parfor to speed up calculation time

function [t_series, x1_t_grid,x2_t_grid,x3_t_grid] = flow_map_curved_ridge(tspan,F1_0,F2_0,F3_0,cpu_num,id)

t_series = tspan;
num_time_points = length(tspan);

% record the solution in parfor
x1_t_parfor = cell(cpu_num,num_time_points);
x2_t_parfor = cell(cpu_num,num_time_points);

parfor i = 1:cpu_num
    
    
    x1_0_list = F1_0(:);
    x1_this_core = x1_0_list(id(i)+1:id(i+1));
    x2_0_list = F2_0(:);
    x2_this_core = x2_0_list(id(i)+1:id(i+1));
    x3_0_list = F3_0(:);
    x3_this_core = x3_0_list(id(i)+1:id(i+1));
    
    % number of points processed by this core
    num_points_this_core = length(x1_this_core);
    
    
    % in order to use ode45, we need to convert the x1_grid, x2_grid and x3_grid to vector
    x123_t0_this_core = [x1_this_core(:);x2_this_core(:);x3_this_core(:)];
    
    
    % tspan is the time points that we want to evaluate
    % each col is the data for one time points
    % different timepoints are in different rows
    [~,x123_t_this_core] = ode45(@velocity_field_curved_ridge,tspan,x123_t0_this_core);
    
    % each col is the data for one time points
    % different timepoints are in different rows
    for t = 1:num_time_points
        x1_t_parfor{i,t} = x123_t_this_core(t,1:num_points_this_core);
        x2_t_parfor{i,t} = x123_t_this_core(t,num_points_this_core+1:2*num_points_this_core);
        x3_t_parfor{i,t} = x123_t_this_core(t,1+2*num_points_this_core:3*num_points_this_core);
    end
    
    

end


% stiching the solution together
x1_t_grid = cell(num_time_points,1);
x2_t_grid = cell(num_time_points,1);
x3_t_grid = cell(num_time_points,1);

parfor t = 1:num_time_points

    % temp variable to combine results from all cores
    x1_t = [];
    x2_t = [];
    x3_t = [];
    for i = 1:cpu_num
        x1_t = [x1_t,x1_t_parfor{i,t}];
        x2_t = [x2_t,x2_t_parfor{i,t}];
        x3_t = [x3_t,x3_t_parfor{i,t}];      
    end
    
    x1_t_grid{t} = reshape(x1_t(:),size(F1_0));
    x2_t_grid{t} = reshape(x2_t(:),size(F2_0));
    x3_t_grid{t} = reshape(x3_t(:),size(F3_0));
    
end

end








