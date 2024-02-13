function cerebellarTractography(templateT1, templateFOD, diedrichsenAtlasFolder, ...
                                workDir, outputDir, fslBinPath, antsBinPath, ...
                                freesurferPath, mrtrixBinPath, wholeBrain)

% Setenv for fsl, ants, and mrtrix. Add path to SPM that contains
% cerebellar scripts
setenv('PATH', [getenv('PATH') ':' fslBinPath]);
setenv('PATH', [getenv('PATH') ':' antsBinPath]);
setenv('PATH', [getenv('PATH') ':' mrtrixBinPath]);
setenv('PATH', [getenv('PATH') ':' freesurferPath]);
setenv('FREESURFER_HOME', freesurferPath)
addpath('/Users/ozenctaskin/Documents/MATLAB/toolboxes/freesurferMatlabLibrary')

% Create workDir and outputDir if they don't exist
if ~exist(workDir, 'dir')
    mkdir(workDir);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% T1 - FOD registration calculation
% Extract B0 slice from FOD to use for registrations. Use FSL.
B0Slice = fullfile(workDir, 'b0FOD.nii.gz');
system(['fslroi ' templateFOD ' ' B0Slice ' 0 1']);

% Now rigidly register T1 to FOD with ANTs
system(['antsRegistrationSyN.sh -m ' templateT1 ' -f ' B0Slice ' -n 7 -t r -o' fullfile(workDir, 'T1toFOD')]);
T1toFODMatrix = fullfile(workDir, 'T1toFOD0GenericAffine.mat');

%% Creating cerebellar masks with SUIT
% The SUIT functions expect non-gunzipped files. So extract the T1 image
T1 = fullfile(workDir, 'anatomy.nii');
system(['mrconvert ' templateT1 ' ' T1]);

% We now apply the inverse transform to the B0 slice to move it over to the
% T1 space. We want to use this image with the T1 image for the cerebellar
% operations and we want to do these operations at the T1 resolution
B0SliceinT1Coor = fullfile(workDir, 'b0InT1coor.nii.gz');
system(['antsApplyTransforms -i ' B0Slice ' -r ' templateT1 ' -t [ ' T1toFODMatrix ',1 ] -o ' B0SliceinT1Coor])
B0SliceinT1CoorUnzipped = fullfile(workDir, 'b0InT1coor.nii');
system(['mrconvert ' B0SliceinT1Coor ' ' B0SliceinT1CoorUnzipped]);

% Do cerebellar isolation. You might need to set "bb" variable for isolate 
% function if the frame is not good in c_. 
suit_defaults();
anatomy = {T1, B0SliceinT1CoorUnzipped};
suit_isolate_seg(anatomy);

% Get filename of the first input to anatomy images
[~, fileName, ~] = fileparts(anatomy{1});

% Normalize to SUIT atlas coordinates
job.subjND.gray = {fullfile(workDir, [fileName '_seg1.nii'])};
job.subjND.white = {fullfile(workDir, [fileName '_seg2.nii'])};
job.subjND.isolation = {fullfile(workDir, ['c_' fileName '_pcereb.nii'])};
suit_normalize_dartel(job)

% Get DiedrichAtlas parcellations and copy the segmentations in workDir
atlas = fullfile(workDir, 'atlas.nii');
system(['cp ' fullfile(diedrichsenAtlasFolder, 'atl-Anatom_space-SUIT_dseg.nii') ' ' atlas]);

% Map atlas back to T1 coordinates
job = [];
job.Affine = {fullfile(workDir, ['Affine_' fileName '_seg1.mat'])};
job.flowfield = {fullfile(workDir, ['u_a_' fileName '_seg1.nii'])};
job.resample = {atlas};
job.ref = {T1};
suit_reslice_dartel_inv(job)
resampledAtlas = fullfile(workDir, ['iw_atlas_u_a_' fileName '_seg1.nii']);

% Create a mask folder where we will save the masks we will use for
% tractography.
maskFolder = fullfile(workDir, 'masks');
if ~exist(maskFolder, 'dir')
    mkdir(maskFolder)
end

% Load the T1 image. We will use it as a template
templateT1Loaded = MRIread(templateT1);

% Load the atlas in FOD coordinates. Permute it to match freesurfer
% convention. We do this because the original atlas cannot be read by the
% freesurfer functions but can be with niftiread.
resampledAtlasLoaded = niftiread(resampledAtlas);
resampledAtlasLoaded = permute(resampledAtlasLoaded, [2, 1, 3]);

% Get the label values for left and right hemi gray matter in cerebellum as
% well as dentate and vermis
maskLabels = {'cerebellumLeftGray', 'cerebellumRightGray', ...
              'dentateLeft', 'dentateRight', 'vermis'};
maskValues = {[1,3,5,8,11,14,17,20,23,26], ...
              [2,4,7,10,13,16,19,22,25,28], [29], [30], ...
              [6,9,12,15,18,21,24,27]};

% Assemble the full masks from the labels, save, and register to FOD
% coordinates
for ii = 1:length(maskLabels)
    mask = zeros(size(resampledAtlasLoaded));
    maskIdx = ismember(resampledAtlasLoaded, maskValues{ii});
    mask(maskIdx) = 1;
    templateT1Loaded.vol = mask;
    saveName = fullfile(maskFolder,[maskLabels{ii} '.nii.gz']);
    MRIwrite(templateT1Loaded, saveName);
    system(['antsApplyTransforms -i ' saveName ' -r ' B0Slice ' -t ' T1toFODMatrix ' -n NearestNeighbor -o ' saveName])
end

%% Make thalamus and motor cortex masks 

%% Tractography. Provide options for whole brain or seed based
if istrue(wholeBrain)
    % Make folder for whole brain results 
    saveFolder = fullfile(outputDir, 'wholeBrain');
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    % Make 5tt segmentations
    segmentations5TT = fullfile(workDir, '5tt_image.nii.gz');
    system(['5ttgen fsl ' templateT1, ' ', segmentations5TT ' -force']);
    system(['antsApplyTransforms -i ' segmentations5TT ' -r ' B0Slice ' -t ' T1toFODMatrix ' -n NearestNeighbor -e 3 -o ' segmentations5TT]);
    
    % Perform whole brain tractography
    wholeBrainTracts = fullfile(saveFolder, 'tracks10mil.tck');
    system(['tckgen -act ' segmentations5TT ' -nthreads 7 –backtrack –select 15000000 -seed_dynamic ' templateFOD ' ' templateFOD ' ' wholeBrainTracts]);

    % SIFT it. 
    siftedTracts = fullfile(saveFolder, 'sift_tracks10mil.tck');
    system(['tcksift -act ' segmentations5TT ' -term_number 1500000 ' wholeBrainTracts ' ' templateFOD ' ' siftedTracts])

end


end
