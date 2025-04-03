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

heights = 0:84:4900;
% heights = 4400;
responses = zeros([length(heights), 1]);
currentmags = zeros([length(heights), 1]);

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
    
    e_curr = calc_elem_current(sim_img, inh_data.volt(:,3));
    currentmags(i) = max(rssq(e_curr.'));

    % Calculate average magnitude of measured response
    postbend = inh_data.meas;
    responses(i) = mean(abs(postbend-prebend));
end

plot((3/4900)*heights, responses, 'linewidth', 2, 'color', 'k');
ylim([0 0.2]);
box off
set(gca, 'linewidth', 2, 'fontsize', 15);
return


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
% sum = 0;
for i = 1:n
    hh=streamline(q.xp,q.yp, q.xc, q.yc, xs(i), sy);
    set(hh,'Linewidth',2, 'color', 'b');
    % if ~isempty(hh)
    %     Vx = interp2(repmat(q.xp, [60 1]), repmat(q.yp, [42 1]).', q.xc, hh.XData, hh.YData);
    %     Vy = interp2(repmat(q.xp, [60 1]), repmat(q.yp, [42 1]).', q.yc, hh.XData, hh.YData);
    %     sum = sum + mean(rssq([Vx; Vy]));
    % end
end
clf
% sum = sum/n
