function J = jet_omega_mod(n)
% Jet colormap modified for log-spaced frequencies (natural frequency maps).
% Colors are loaded from jetxmod_freqlog.mat (same folder).
% Input n is ignored; the colormap is returned with its stored length.

load('jetxmod_freqlog.mat', 'cmap');
J = cmap;
