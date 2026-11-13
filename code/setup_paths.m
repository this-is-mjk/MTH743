% setup_paths: add lib/ and experiments/ to the MATLAB path
function setup_paths()
    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'lib'));
    addpath(fullfile(root,'tests'));
    addpath(fullfile(root,'experiments'));
    fprintf('paths ready: %s\n', root);
end
