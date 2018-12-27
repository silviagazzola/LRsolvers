%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IR Tools
% December 2018
% 
% This file is part of the IR Tools package and is distributed under the
% 3-Clause BSD License. A separate license file should be provided as
% part of the package. 
% 
% Copyright 2018 Silvia Gazzola, University of Bath, Per Christian Hansen,
% Technical University of Denmark and James G. Nagy, Emory University.
%
% The present version is extended by J. Chung, Virginia Tech, and 
% S. Gazzola, University of Bath
% 
% Contact: jmchung@vt.edu, s.gazzola@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Installating
%   IRtools_setup       - Set up search paths to IRtools
%
% Example scripts:
%   Example 1
%   Example 2
%
% Iterative methods:
%   IRart              - Algebraic Reconstruction Technique
%   IRcgls             - Conjugate Gradient algorithm for Least Squares problems
%   IRconstr_ls        - Least squares solver with box and energy constraints
%   IRell1             - Least squares solver with 1-norm penalization term
%   IRenrich           - Enriched CGLS method
%   IRfista            - FISTA algorithm for constrained least squares problems 
%                        and for 1-norm penalization
%   IRget              - Get options for IR Tools functions
%   IRhtv              - Least squares solver with heuristic total variation penalization
%   IRhybrid_fgmres    - Hybrid version of FGMRES for enforcing 1-norm penalization
%   IRhybrid_flsqr     - Hybrid version of FLSQR for enforcing 1-norm penalization
%   IRhybrid_gmres     - Hybrid version of GMRES algorithm for square systems
%   IRhybrid_lsqr      - Hybrid version of LSQR algorithm
%   IRirn              - Least squares solver with 1-norm penalization term
%   IRmrnsd            - Modified Residual Norm Steepest Descent method
%   IRnnfcgls          - Modified flexible CGLS for nonnegatively constrained LS problems
%   IRrestart          - Restarted Krylov subspace methods
%   IRrrgmres          - Range Restricted GMRES for square systems
%   IRset              - Set options for IR Tools functions
%   IRsirt             - Simuletaneous Iterative Reconstruction Technique
%
% Operators (functions) for some test problems:
%   OPdiffusion        - The forward computation and its adjoint for PRdiffusion
%   OPinvinterp2       - The forward computation and its adjoint for PRinvinterp2
%   OPnmr              - The forward computation and its adjoint for PRnmr
%   OPblur_transform   - The forward computation and its adjoint for PRblur, incorporating sparsity transform
