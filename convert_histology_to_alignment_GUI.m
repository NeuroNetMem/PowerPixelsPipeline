% This script converts the tracked points from Universal Probe Finder into
% insertion trajectories which can be read in by the IBL alignment GUI.

% load atlas
sAtlasParams = PF_getAtlasIni();
cellAtlases = {sAtlasParams.name};
[intSelectAtlas, boolContinue] = listdlg('ListSize', [200, 100], 'Name', 'Atlas Selection', ...
    'PromptString', 'Select Atlas:', 'SelectionMode', 'single', 'ListString', cellAtlases);
strAtlasName = sAtlasParams(intSelectAtlas).name; 
strPathVar = sAtlasParams(intSelectAtlas).pathvar;
fLoader = sAtlasParams(intSelectAtlas).loader;
strAtlasPath = PF_getIniVar(strPathVar);
sAtlas = feval(fLoader,strAtlasPath);

% Load in slice file
[strSliceFile, strSlicePath] = uigetfile('.mat', 'Open slice file');
sLoad = load(fullpath(strSlicePath, strSliceFile));
sSlice = sLoad.sSliceData;
sTracks = SF_SliceFile2TracksFile(sLoad.sSliceData, sAtlas);
sTracks.Type = 'native';
sProbeCoords = PH_ExtractProbeCoords(sTracks);

% Correct for the Z axis in the Bregma coordinates
vecBregma = [sAtlas.Bregma(1), sAtlas.Bregma(2), 770];

% Loop over tracks
fprintf('Found %d probe tracks\n', numel(sProbeCoords.cellPoints));
for i = 1:numel(sProbeCoords.cellPoints)
    % Set probe index
    sProbeCoords.intProbeIdx = i;
    
    % Get an inital fit of the vector to calculate the probe length 
    [vecSphereVectorInit, ~, ~] = PH_Points2vec(sProbeCoords, sAtlas);
    dblZdist = max(sProbeCoords.cellPoints{i}(:,3)) - min(sProbeCoords.cellPoints{i}(:,3));
    sProbeCoords.ProbeLength = dblZdist / cosd(vecSphereVectorInit(4));

    % Fit the probe vector and get the cartesian coordinates of the base and tip
    [vecSphereVector, vecLocBrainIntersect, matRefVector] = PH_Points2vec(sProbeCoords, sAtlas);
    probe_vector_cart = PH_SphVec2CartVec(vecSphereVector);

    % Convert the cartesian coordinates into Bregma coordinates
    matBregmaPoints = zeros(size(probe_vector_cart));
    for j = 1:3
        matBregmaPoints(:, j) = (vecBregma(j) - probe_vector_cart(:, j)) * sAtlas.VoxelSize(j);
    end
    matBregmaPoints(:,2) = -matBregmaPoints(:,2);
    matBregmaPoints(:,3) = -matBregmaPoints(:,3);

    % Save to disk as JSON file
    sXYZ = struct('xyz_picks', matBregmaPoints);
    fid = fopen(strcat(strSlicePath , filesep, sTracks.cellNames{i}, '.json'), 'w');
    fprintf(fid, jsonencode(sXYZ)); 
    %writematrix(matBregmaPoints, strcat(strSlicePath , filesep, sTracks.cellNames{i}, '.csv'))
    fprintf('Exported %s as %s\n', sProbeCoords.cellNames{i}, ...
        strcat(strSlicePath, filesep, sTracks.cellNames{i}, '.json'));
end


