% Do a rigid rotation/translation/mirror of a set of n points theta so that
% they are well superimposed on a corresponding set of n points
% theta_actual. The measure of good superimposition is lrms (least root
% mean square deviation)
% Usage:
% [theta_translated, lrms] = superimpose(theta, pattern,prune)
% where:
% theta and prune are the outputs from solveDGP,
% pattern contains a list of coordinates specifying the designed pattern,
% theta_translated is the superimposition of the final embedding that minimizes the RMSD between the designed and reconstructed pattern, and
% lrms is the corresponding RMSD.

function [theta_corrected, lrms] = superimpose(theta, pattern, prune, paint)

if nargin == 3
    paint = 'none';
end

theta_actual = pattern(:,2:3);
pattern_trimmed = pattern(prune,:);
theta_actual_trimmed = theta_actual(prune,:);
theta_actual_pruned = theta_actual(~prune,:);


[U, r, lrms] = Kabsch(theta',theta_actual_trimmed'); % minimize RMSD
%[U2, r2, lrms2] = rigidtform(theta,theta_actual); % minimize RMSD
theta_corrected = U*theta'+r;
theta_corrected = theta_corrected';

% flip left-right and superimpose
theta_fliplr = theta; % copy
theta_fliplr(:,1) = -theta(:,1); % flip
[U_fliplr, r_fliplr, lrms_fliplr] = Kabsch(theta_fliplr',theta_actual_trimmed'); % minimize RMSD of flipped version

if lrms_fliplr < lrms % if the flipped version is a better superimposition then keep that and discard earlier one
    %U = U_fliplr;
    % r = r_fliplr;
    lrms =  lrms_fliplr;
    theta_corrected =  U_fliplr*theta_fliplr' + r_fliplr;
    theta_corrected = theta_corrected';
end

dp = theta_actual_trimmed - theta_corrected;
delta_dp = sqrt(dp(:,1).^2 + dp(:,2).^2);

filter_large_dp = delta_dp > 2.5;
theta_corrected_filtered = theta_corrected(filter_large_dp,:);
dp_filtered = dp(filter_large_dp,:);

% %cmap = hsv(216);
% %cmap_full = linspecer(216);
% saturation_ratio = 0.75;
% cmap_full = distinguishable_colors(216);
% cmap_full = rgb2hsv(cmap_full);
% cmap_full(:,2) = saturation_ratio*cmap_full(:,2);
% cmap_full = hsv2rgb(cmap_full);
% cmap = cmap_full(pattern_trimmed(:,1),:);

%% for color figures
if isequal(paint,'tree')
    cmap = paint_tree(pattern);
    cmap = cmap(prune,:);
elseif isequal(paint,'wheel')
    cmap = paint_wheel(pattern);
    cmap = cmap(prune,:);
end

% define desaturated red
desat_red = rgb2hsv([1 0 0]);
desat_red(2) = 0.75*desat_red(2);
desat_red = hsv2rgb(desat_red);

% figure;
hold on;axis equal;
padding = 5;
xmin = min(theta_actual(:,1))-padding; xmax = max(theta_actual(:,1))+padding;
ymin = min(theta_actual(:,2))-padding; ymax = max(theta_actual(:,2))+padding;
xlim([xmin xmax]); ylim([ymin ymax]);
plot([xmax-2*padding xmax - padding], [ymin+padding ymin+padding],'LineWidth',10,'Color','k'); % scale bar

% old_units = get(gca, 'Units'); % note down units for safekeeping
% set(gca, 'Units', 'points'); % change units to points
% pos_points = get(gca, 'Position'); % get dimensions of the axes
% set(gca, 'Units', old_units); % switch back to old units
% marker_size = pi*(4*pos_points(3)/(xmax-xmin)).^2/4; % marker has diameter 4 nm
marker_size = 100;

scatter(theta_actual_pruned(:,1),theta_actual_pruned(:,2), marker_size,desat_red,'LineWidth',2);
scatter(theta_actual_trimmed(:,1),theta_actual_trimmed(:,2), marker_size,[0.75 0.75 0.75],'filled','LineWidth',2);
if isequal(paint,'none')
    scatter(theta_corrected(:,1),theta_corrected(:,2), marker_size, 'k','LineWidth',2);
else
    scatter(theta_corrected(:,1),theta_corrected(:,2), marker_size, cmap,'filled','MarkerEdgeColor','k');
end
quiver(theta_corrected_filtered(:,1),theta_corrected_filtered(:,2),dp_filtered(:,1),dp_filtered(:,2),0,'color','k','LineWidth',0.5,'ShowArrowHead','off');


%text(1+theta_corrected(:,1), 1+theta_corrected(:,2),num2str([1:size(theta_corrected,1)]'),'FontSize',18);
set(gca,'XColor', 'none','YColor','none')
title(['RMSD = ', num2str(lrms,2), ' nm']);
hold off
end