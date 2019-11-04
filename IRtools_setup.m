function IRtools_setup
%IRtools_setup Set up search paths to IR Tools with low-rank solvers
%
%  Run this function to setup IR Tools with low-rank solvers. 
%
%  In most cases, this script should only need to be run once, and then the 
%  path will be permanently saved. Note the following:
%    - It forces all user paths in the current session to set permanently.
%    - On systems with shared MATLAB licenses, a permanent save may need to
%      be done manually. In this case, see SAVEPATH for more information.
%  
%  If you want to remove the paths associated with IR Tools, use RMPATH.
%
% See also: addpath, rmpath, savepath

% Silvia Gazzola, University of Bath

addpath(genpath(fileparts(mfilename('fullpath'))));

status = savepath;
if status == 1
    warning('IR Tools was added to the MATLAB search path for the current session only. Adding it permanently failed, probably due to a write permission issue. It is possible to manually add IR Tools permanently to your search path, but may require consulting your system administrator. Alternatively, you can re-run this installation function in each new MATLAB session where you want to use IR Tools.')
end
