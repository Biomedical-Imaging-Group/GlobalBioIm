%% Script which updates Matlab path to work with the GlobalBioIm Library
[lpath,~,~] = fileparts(which('setGlobalBioImPath'));
addpath(genpath(lpath));



