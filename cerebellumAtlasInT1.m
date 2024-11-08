function cerebellumAtlasInT1(inputT1, inputT2, MNItemplate, parcellations)

% Warp cerebellum parcellation to anatomical coordinates. This script also
% linearly registers T2 to T1 image using SPM functions. 
%
% Inputs:
%   inputT1: T1 anatomical image path. 
%   inputT2: T2 anatomical path. Optional. Enter 'NA' if you want to do the
%            reconstruction just with the T1 image
%   MNItemplate: Path to an MNI template. One can be found in fsl directory
%   (e.g. /fsl/data/standard/MNI152_T1_1mm.nii.gz. Use the isotropic
%   template that is in the same resolution with you input T1. If your T1
%   is 1mm iso use MNI152_T1_1mm.nii.gz. If it is 2mm iso, use MNI152_T1_2mm.nii.gz
%   if Na is passed to MNI template, the ACPC registration step is skipped.
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

% Create a cerebellum script workdir where T1 image is located
[folderPath, ~, ~] = fileparts(inputT1);
workdir = fullfile(folderPath, 'cerebellumWorkdir');
if ~isfolder(workdir)
    system(['mkdir ' workdir]);
end

% Unzip if images are nii.gz so we end up with .nii
[T1path, T1name, T1extension] = fileparts(inputT1);
if strcmp(T1extension, '.gz')
    gunzip(inputT1, workdir);
    inputT1 = fullfile(workdir, T1name);
end
if ~strcmp(inputT2, 'NA')
    [T2path, T2name, T2extension] = fileparts(inputT2);
    if strcmp(T2extension, '.gz')
        gunzip(inputT2, workdir);
        inputT2 = fullfile(workdir, T2name);
    end
end
if ~strcmp(MNItemplate, 'NA')
    [MNIpath, MNIname, MNIextension] = fileparts(MNItemplate);
    if strcmp(MNIextension, '.gz')
        gunzip(MNItemplate, workdir);
        MNItemplateUzipped = fullfile(workdir, MNIname);
    end
end

% Do a linear registration between T1 and MNI template so that we put the
% T1 image in the ACPC orientation. Do this only if an MNI template is
% passed and not Na
if ~strcmp(MNItemplate, 'Na')
    % Get a copy of the T1 that we will use for reslicing at the end
    inputT1_copy = fullfile(workdir, ['copy_' T1name]);
    system(['cp ' inputT1 ' ' inputT1_copy]);

    % Calculate registration
    inputT1Loaded = spm_vol(inputT1);
    MNItemplateLoaded = spm_vol(MNItemplateUzipped);
    x = spm_coreg(inputT1Loaded, MNItemplateLoaded);
    M = spm_matrix(x);
    spm_get_space(inputT1, M * inputT1Loaded.mat);
    
    % Apply registration to create an ACPC T1
    flags= struct('interp',1,'mask',1,'mean',0,'which',1,'wrap',[0 0 0], 'prefix', 'acpc_');
    files = {MNItemplateUzipped;inputT1};
    spm_reslice(files, flags);
    acpcT1 = fullfile(workdir, ['acpc_' T1name]);
else
    acpcT1 = inputT1;
end

% Do a linear registration between T2 and T1 images using SPM. That's if a
% T2 image is inputted.
if ~strcmp(inputT2, 'Na')
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {acpcT1};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {inputT2};
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'acpc_';
    spm_jobman('run', matlabbatch);
    [T2path, T2name, T2extension] = fileparts(inputT2);
    acpcT2 = fullfile(T2path, ['acpc_' T2name T2extension]);
end

% Setup the anatomy cell for SUIT processing 
if ~strcmp(inputT2, 'NA')
    anatomy = {acpcT1, acpcT2};
else
    anatomy = {acpcT1};
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

% Combine affine matrix we get from T1-MNI registration with the one
% calculated by SUIT so that we can map atlas back to original subject
% coordinates 
load(fullfile(filePath, ['Affine_' fileName '_seg1.mat']));
Affine = Affine * M;
save(fullfile(filePath, ['Affine_' fileName '_seg1.mat']), 'Affine');

% Map atlas back to T1 coordinates
job = [];
job.Affine = {fullfile(filePath, ['Affine_' fileName '_seg1.mat'])};
job.flowfield = {fullfile(filePath, ['u_a_' fileName '_seg1.nii'])};
job.resample = {parcellations};
job.ref = {inputT1_copy};
suit_reslice_dartel_inv(job)

% Move the atlas folder to the main directory
files = dir(fullfile(filePath, 'iw*'));
resampledAtlas = fullfile(filePath, files(1).name);
system(['mv ' resampledAtlas ' ' fullfile(folderPath, files(1).name)]);

end