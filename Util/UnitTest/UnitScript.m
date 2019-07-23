%--------------------------------------------------------------
%            Unitary Tests for the GlobalBioIm Library
%
% USAGE :
%    - set generateDataUnitTests=1; to generate the data to compare (e.g.
%    from an older version of the Library
%    - set generateDataUnitTests=0; to run the scripts and compare their
%    results to those saved (e.g. from the current version of the Library)
%    - use the structure 'test' to choose which scripts are used
%
%  Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%--------------------------------------------------------------
% Initializations
clc; clear; close all;
global generateDataUnitTests  stateTest message
generateDataUnitTests=0; 
mess={'FAIL !','OK !'};

% Select scripts to run
test.Deconv_LS_NoReg=1;
test.Deconv_LS_NonNeg_NoReg=1;
test.Deconv_LS_TV=1;
test.Deconv_LS_TV_NonNeg=1;
test.Deconv_LS_HessSchatt=1;
test.Deconv_LS_HessSchatt_NonNeg=1;
test.Deconv_KL_NonNeg_NoReg=1;
test.Deconv_KL_TV_NonNeg=1;
test.Deconv_KL_HessSchatt_NonNeg=1;
test.TestProxL2SumConv=1;
test.TestProxL2DownSampledConv=1;
test.TestsSummationLinOps=1;
test.TestsCompositionLinOps=1;

% Put 1 to run test when test.??=0
revert=0;

[lpath,~,~] = fileparts(which('setGlobalBioImPath'));
pth = cd;
cd(lpath);

warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir('Util/UnitTest/Data');
warning('on','MATLAB:MKDIR:DirectoryExists');

tt=tic;
fnames=fieldnames(test);
for idx=1:length(fnames)
    message=[];
    if (test.(fnames{idx}) && ~revert) || (~test.(fnames{idx}) && revert) 
        evalc(['run ',fnames{idx}]);close all;
        if generateDataUnitTests
            fprintf(' %-30s  --> results saved \n',fnames{idx});
        else
            if stateTest
                fprintf(' %-30s  --> %s \n',fnames{idx},mess{stateTest+1});
            else
                fprintf(' %-30s  --> %s \n \t %s\n',fnames{idx},mess{stateTest+1},message);
            end
        end
    end
end
toc(tt);
cd(pth);
