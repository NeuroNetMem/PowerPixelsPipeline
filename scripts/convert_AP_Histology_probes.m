% This script converts the tracked points from Universal Probe Finder into
% insertion trajectories which can be read in by the IBL alignment GUI.

% Set bregma
vecBregma = [540, 44, 570];

% Load in probe file
[strProbeFile, strProbePath] = uigetfile('.mat', 'Open probe_ccf file');
sLoad = load(fullfile(strProbePath, strProbeFile));
probe_ccf = sLoad.probe_ccf;

% Calculate best fit through points
for curr_probe = 1:length(probe_ccf)

    % make sure the direction goes down in DV - flip if it's going up
    r0 = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    % Fit line to the extend of the placed points in DV
    [~, highest_point] = min(probe_ccf(curr_probe).points(:, 2));
    line_eval = [0, max(probe_ccf(curr_probe).points(:, 2))];
    probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval', histology_probe_direction'), ...
        probe_ccf(curr_probe).points(highest_point, :));

    % Convert fitted trajectory to Bregma coordinates
    matBregmaPoints = ([vecBregma; vecBregma] - probe_fit_line) * 10;  % 10 um voxels
    matBregmaPoints = [matBregmaPoints(:, 3), matBregmaPoints(:, 1), matBregmaPoints(:, 2)];

    % Save to disk as JSON file
    sXYZ = struct('xyz_picks', matBregmaPoints);
    fid = fopen(strcat(strProbePath, ['probe', num2str(curr_probe), '.json']), 'w');
    fprintf(fid, jsonencode(sXYZ)); 
    %writematrix(matBregmaPoints, strcat(strSlicePath , filesep, sTracks.cellNames{i}, '.csv'))
    fprintf('Exported %s as %s\n', ['probe', num2str(curr_probe), '.json'], ...
        strcat(strProbePath, ['probe', num2str(curr_probe), '.json']));
end
