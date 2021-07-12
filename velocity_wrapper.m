%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate velocity (velocityField function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation_type (velcotyField argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation_type = 1  -> time smoothed instantaneous velocity.
% i.e. to calculate the velocity across 5 time points, we'd do:
% (((t-1)-(t-2)) + ((t)-(t-1)) + ((t+1)-(t)) + ((t+2)-(t+1))) / 4

% calculation_type = 2 -> central differencing of total displacement.
% i.e. to calculate the velocity across 5 time points, we'd do:
% ((t+2)-(t-2))

% NOTE: both velocity calculations are symmetric with time.
% i.e. velocity across 5 time points is (t-2):(t+2)

% specify argument
calculation_type = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displacement_time (velocityField argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The time-interval over which velocity is calculated
% i.e. 5 is (t-2):(t+2)

%%% NOTE: the only options are 3,5,7,9,11,13 %%%

% Because the frame rate is 1frame/5min, this corresponds to 15, 25, 35,
% 55, 65, 75 min.  For reference, one column comes off the MF every 2.5 hours.

displacement_time = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interp_subsample_factor (velocityField argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spatial scale onto which the cell-centroid based velocity field will be
% interpolated

interp_subsample_factor = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocity_scaling_factor (velocityField argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% factor by which velocity vectors will be scaled (for visualization
% purposes)

velocity_scaling_factor = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interp_time_range (velocityField argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interp_time_range = 60;


%%%%%%%%%%%%%%%
% call function
%%%%%%%%%%%%%%%
[Velocity,calculation_time_range] = velocityField(L, cell_area, cell_centroid, tracks ,...
    calculation_type, displacement_time, interp_subsample_factor, ...
    velocity_scaling_factor,interp_time_range);

%%%%%%%%%%%%%%%%%%%%%
% outputs of function
%%%%%%%%%%%%%%%%%%%%%

% x_coordinates/y_coordinates:
% cell array containing a Nx3 matrix for each time point.  Each Nx3 matrix
% is structured as [centroid-component, 

% x_velocities/y_velocities:
% cell array containing a Nx3 matrix for each time point.  Each Nx3 matrix
% is structure as [centroid-component,
























