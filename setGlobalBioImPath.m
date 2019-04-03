%% Script which update matlab path to work with the  GlobalBioIm Library
[lpath,~,~] = fileparts(which('setGlobalBioImPath'));
addpath(genpath(lpath));



