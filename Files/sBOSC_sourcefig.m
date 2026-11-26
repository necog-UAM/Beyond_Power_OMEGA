function sBOSC_sourcefig(datasource, cfg)

%%
% Plots source data on the inflated cortical surface (both hemispheres, lateral and medial views).
% Fieldtrip required. Can input color lims as cfg.colim.
%   sBOSC_sourcefig(datasource) plots the source
%   sBOSC_sourcefig(datasource, cfg) plots the data applying custom 
%   plot configurations.
%
%   Inputs:
%           datasource - Numeric vector containing exactly 1925 or 3423 voxels.
%
%           cfg        - (Optional) Structure containing visual configurations:
%
%               cfg.interp  : String defining the interpolation method (e.g.,
%                           'nearest', 'linear'). Default is 'nearest'.
%               cfg.colim   : 1x2 numeric array for color limits [min max].
%                           Default is 'auto'.
%               cfg.colmap  : String or N x 3 array defining the colormap.
%                           Default is 'auto'.


% FieldTrip is required
ftPath = fileparts(which('ft_defaults'));
if isempty(ftPath)
    error('sBOSC:MissingDependency', 'FieldTrip not found. Add it to the path and run ft_defaults.');
end

% Defaults
if nargin < 2
    cfg = struct();
end
if isfield(cfg, 'interp')
    interpmethod = cfg.interp;
else
    interpmethod = 'nearest';
end
if isfield(cfg, 'colim')
    funcolorlim = cfg.colim;
else
    funcolorlim = 'auto';
end
if isfield(cfg, 'colmap')
    funcolormap = cfg.colmap;
else
    funcolormap = 'auto';
end

% Correct column orientation
datasource = datasource(:);
Nvox = length(datasource);

if Nvox ~= 1925 && Nvox ~= 3423
    error('Input data must have exactly 1925 or 3423 voxels.');
end

% Load sources
if Nvox == 1925
    t = load('source_template_10mm_1925.mat');
else
    t = load('source_template_10mm_3423.mat');
end
inside = t.source.inverse.inside;

mdl = load(fullfile(ftPath, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm.mat'));
source_mni = mdl.sourcemodel;
source_mni.inside = inside;
source_mni.avg.pow = nan(length(inside), 1);
source_mni.avg.pow(inside) = datasource;
source_mni.time = 1;

templ = load(fullfile(ftPath, 'template', 'headmodel', 'standard_mri.mat'));

% Interpolation
cfg_interp = [];
cfg_interp.parameter  = 'avg.pow';
cfg_interp.downsample = 2;
cfg_interp.interpmethod = interpmethod;
source_interp = ft_sourceinterpolate(cfg_interp, source_mni, templ.mri);

% Plot
figure('WindowState','maximized','Color','w');

cfg_plot = [];
cfg_plot.method        = 'surface';
cfg_plot.funparameter  = 'pow';
cfg_plot.maskparameter = 'pow';
cfg_plot.projmethod    = 'nearest';
cfg_plot.opacity       = 0.8;
cfg_plot.figure        = 'gca';
cfg_plot.camlight      = 'no';
cfg_plot.colorbar      = 'yes';
cfg_plot.funcolorlim = funcolorlim;
cfg_plot.funcolormap = funcolormap;

% Left Hemisphere
cfg_plot.surffile     = fullfile(ftPath, 'template', 'anatomy', 'surface_pial_left.mat');
cfg_plot.surfinflated = 'surface_inflated_left_caret_white.mat';
subplot(2,2,1); ft_sourceplot(cfg_plot, source_interp); view([-90 0]); camlight('left');
subplot(2,2,3); ft_sourceplot(cfg_plot, source_interp); view([90 0]);  camlight('left');

% Right Hemisphere
cfg_plot.surffile     = fullfile(ftPath, 'template', 'anatomy', 'surface_pial_right.mat');
cfg_plot.surfinflated = 'surface_inflated_right_caret_white.mat';
subplot(2,2,2); ft_sourceplot(cfg_plot, source_interp); view([90 0]);  camlight('right');
subplot(2,2,4); ft_sourceplot(cfg_plot, source_interp); view([-90 0]); camlight('right');
