% test for fast marching and geodesic extraction
%
%   Copyright (c) 2004 Gabriel Peyr?

n = 600;
name = 'road2';
name = 'mountain';
name = 'constant';
[M,W] = load_potential_map(name, n);

rep = 'results/geodesic-2d/';
if not(exist(rep))
    mkdir(rep);
end

warning off;
imwrite(rescale(W), [rep name '-map.png'], 'png');
warning on;

% pick starting point
[start_points] = pick_start_end_point(M);

options.nb_iter_max = Inf;
disp('Performing front propagation.');
[D,S] = perform_fast_marching_2d(W, start_points, options);

npaths = 30;
end_points = floor( rand(npaths,2)*(n-1) )+1;

disp('Extract paths');
paths = {};
for i=1:npaths
    paths{i} = compute_geodesic(D,end_points(i,:)');
end

ms = 30; lw = 3;
% display
A = convert_distance_color(D);
clf; hold on;
imageplot(A); axis image; axis off;
for i=1:npaths
    end_point = end_points(i,:);
    h = plot( paths{i}(:,2), paths{i}(:,1), 'k' );
    set(h, 'LineWidth', lw);    
    h = plot(end_point(2),end_point(1), '.b');
    set(h, 'MarkerSize', ms);    
end
h = plot(start_points(2),start_points(1), '.r');
set(h, 'MarkerSize', ms);
hold off;
colormap jet(256);
axis ij;
saveas(gcf, [rep name '-geodesics.png'], 'png');