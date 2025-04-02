% Create rectangular mesh dimensions
n_elec = [8, 1];
xy_size = [42, 60];
xy_size = xy_size + 1;

xvec = linspace(-1,1,xy_size(1));
height = 6; % This can be changed
yvec = linspace(-height/2,height/2,xy_size(2));
fmdl = mk_grid_model([],xvec,yvec);

% Assign electrodes to top surface
options = {'no_meas_current','no_rotate_meas'};
tb_elecs= linspace(1, xy_size(1), 8); 
el_nodes= [];
bdy_nodes= (1:xy_size(1)) + xy_size(1)*(xy_size(2)-1); 
el_nodes= [el_nodes, bdy_nodes(tb_elecs)];
for i=1:n_elec(1)
    n = el_nodes(i);
    fmdl.electrode(i).nodes= n;
    fmdl.electrode(i).z_contact= 0.001; % choose a low value
end
sim_img= mk_image(fmdl,1);

% Define electrode configurations
stim = mk_stim_patterns(8, 1, '{op}', '{ad}', {}, 1);
sim_img.fwd_model.stimulation = stim;

% heights = 68:84:4900;
heights = 4400;
responses = zeros([length(heights), 1]);

for i = 1:length(heights)

    % First measurements: before bending occurs
    % Define homogeneous conductivity field
    for j = 1:84
        sim_img.elem_data(j:84:end) = 1;
    end

    % Simulate effect of film: set conductivity to zero at centre
    stoppoint = 4981-heights(i); % can be changed to simulate height of film
    sim_img.elem_data(4997:-84:stoppoint)= 0;
    sim_img.elem_data(4998:-84:stoppoint)= 0;
    
    % Solve for expected signals
    sim_img.fwd_solve.get_all_meas = 1;
    inh_data=fwd_solve(sim_img);
    prebend = inh_data.meas;
    
    % Second measurements: after bending occurs
    % Define linear conductivity field to represent effect of bending - can be changed
    for j = 1:84
        sim_img.elem_data(j:84:end) = 1 + (j-42)/100;
    end
    
    % Simulate effect of film: set conductivity to zero at centre
    sim_img.elem_data(4997:-84:stoppoint)= 0;
    sim_img.elem_data(4998:-84:stoppoint)= 0;
    
    % Plot FEM conductivity mesh that will be solved
    show_fem(sim_img);
    
    % Solve for expected signals
    sim_img.fwd_solve.get_all_meas = 1;
    inh_data=fwd_solve(sim_img);
    % plot(inh_data.meas);
    % title("Measurements");
    
    % Calculate average magnitude of measured response
    postbend = inh_data.meas;
    responses(i) = mean(abs(postbend-prebend));

end

% plot(responses);
% return


%% Plot Equipotentials and Streamlines
PLANE= [inf,inf,0]; 
sim_img.fwd_model.mdl_slice_mapper.npx = 64;
sim_img.fwd_model.mdl_slice_mapper.npy = 64;
sim_img.fwd_model.mdl_slice_mapper.level = PLANE;

q = show_current(sim_img, inh_data.volt(:,3)); % Change the number here to change which injection pattern being used

sim_img = rmfield(sim_img, 'elem_data');
sim_img.node_data = inh_data.volt(:,3);
sim_img.calc_colours.npoints = 256;
imgs = calc_slices(sim_img,PLANE);

% sx = -0.4;
sy = 2.9; % Change these: starting points of the current streamline
n = 15; % Number of streamlines
xs = linspace(-1, -0.14, n);

% subplot(1,2,1); hold on
for i = 1:n
    hh=streamline(q.xp,q.yp, q.xc, q.yc, xs(i), sy);
    set(hh,'Linewidth',2, 'color', 'b');
    % hh=streamline(q.xp,q.yp,-q.xc,-q.yc, sx+r*cos(thetas(i)), sy-r*sin(thetas(i)));
    % set(hh,'Linewidth',2, 'color', 'b');
end

% pic = shape_library('get','adult_male','pic');
% [x, y] = meshgrid( linspace(pic.X(1), pic.X(2),size(imgs,1)), ...
%                   linspace(pic.Y(2), pic.Y(1),size(imgs,2)));
% hold on;
% contour(x,y,imgs,61, 'color', 1/255*[217 95 2]);