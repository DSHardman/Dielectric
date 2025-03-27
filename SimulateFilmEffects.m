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
% Arsen: you will need to change this line to match the electrode configurations
% you used in experiments (there is good documentation for mk_stim_patterns online):
stim = mk_stim_patterns(8, 1, '{op}', '{ad}', {}, 1);
sim_img.fwd_model.stimulation = stim;

% Define linear conductivity field to represent effect of bending - can be changed
for i = 1:84
    sim_img.elem_data(i:84:end) = 1 + (i-42)/200;
end

% Simulate effect of film: set conductivity to zero at centre
stoppoint = 500; % can be changed to simulate height of film
sim_img.elem_data(4997:-84:stoppoint)= 0;
sim_img.elem_data(4998:-84:stoppoint)= 0;


% Plot FEM conductivity mesh that will be solved
subplot(1,2,1);
show_fem(sim_img);

% Solve for expected signals and plot these
subplot(1,2,2);
sim_img.fwd_solve.get_all_meas = 1;
inh_data=fwd_solve(sim_img);
plot(inh_data.meas);
title("Measurements");


%% Plot Equipotentials and Streamlines
PLANE= [inf,inf,0]; 
sim_img.fwd_model.mdl_slice_mapper.npx = 64;
sim_img.fwd_model.mdl_slice_mapper.npy = 64;
sim_img.fwd_model.mdl_slice_mapper.level = PLANE;

q = show_current(sim_img, inh_data.volt(:,1)); % Change the number here to change which injection pattern being used

sim_img = rmfield(sim_img, 'elem_data');
sim_img.node_data = inh_data.volt(:,1);
sim_img.calc_colours.npoints = 256;
imgs = calc_slices(sim_img,PLANE);

sx = -0.9; sy = 2.9; % Change these: starting points of the current streamline

subplot(1,2,1); hold on
hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',3, 'color', 'b');
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',3, 'color', 'b');
