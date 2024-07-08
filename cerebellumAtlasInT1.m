function cerebellumAtlasInT1(inputT1, inputT2, parcellations)

% Warp cerebellum parcellation to anatomical coordinates. This script also
% linearly registers T2 to T1 image using SPM functions. 
%
% Inputs:
%   inputT1: T1 anatomical image path. 
%   inputT2: T2 anatomical path. Optional. Enter 'NA' if you want to do the
%            reconstruction just with the T1 image
%   parcellations: SUIT Cerebellum parcellations. You can download one from
%   here: https://github.com/DiedrichsenLab/cerebellar_atlases/tree/master.
%   e.g. in Diedrichsen, Ji, Buckner. Make sure to input the files with 
%   "SUIT_dseg" extension. 
%
% Outputs are saved in the original input image path. The final product has
% the "iw_atl-" in the name.
%

% Set suit and SPM defaults
suit_defaults();
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Unzip if images are nii.gz so we end up with .nii
[T1path, T1name, T1extension] = fileparts(inputT1);
if strcmp(T1extension, '.gz')
    gunzip(inputT1)
    inputT1 = fullfile(T1path, T1name);
end

if ~strcmp(inputT2, 'NA')
    [T2path, T2name, T2extension] = fileparts(inputT2);
    if strcmp(T2extension, '.gz')
        gunzip(inputT2)
        inputT2 = fullfile(T2path, T2name);
    end
end

% Do a linear registration between T2 and T1 images using SPM. That's if a
% T2 image is inputted.
if ~strcmp(inputT2, 'NA')
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {inputT1};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {inputT2};
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'registered';
    spm_jobman('run', matlabbatch);
    [T2path, T2name, T2extension] = fileparts(inputT2);
    inputT2 = fullfile(T2path, ['registered' T2name T2extension]);
end

% Setup the anatomy cell for SUIT processing 
if ~strcmp(inputT2, 'NA')
    anatomy = {inputT1, inputT2};
else
    anatomy = {inputT1};
end

% Get T1 and T2 and segment cerebellumn
suit_isolate_seg(anatomy);

% Get filename of the first input to anatomy images
[filePath, fileName, ~] = fileparts(anatomy{1});

% Normalize to SUIT atlas coordinates
job.subjND.gray = {fullfile(filePath, [fileName '_seg1.nii'])};
job.subjND.white = {fullfile(filePath, [fileName '_seg2.nii'])};
job.subjND.isolation = {fullfile(filePath, ['c_' fileName '_pcereb.nii'])};
suit_normalize_dartel(job)

% Map atlas back to T1 coordinates
job = [];
job.Affine = {fullfile(filePath, ['Affine_' fileName '_seg1.mat'])};
job.flowfield = {fullfile(filePath, ['u_a_' fileName '_seg1.nii'])};
job.resample = {parcellations};
job.ref = {anatomy{1}};
suit_reslice_dartel_inv(job)
resampledAtlas = fullfile(filePath, ['iw_atlas_u_a_' fileName '_seg1.nii']);

end