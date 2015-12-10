% aquaDOM Toolbox for MATLAB (v1.0.1)
% December 2015
%
% Copyright Urban Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Kavalergården 6
% 2920 Charlottenlund, Denmark
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License <http://www.gnu.org/licenses/>
%
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
% Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
%-----------------------------------------------
% GENERAL DESCRIPTION
% aquaDOM is a MATLAB toolbox for the calculation of apparent
% fluorescence quantum yields of dissolved organic matter. It contains
% functions that can utilize EEMs that were initially processed by the
% drEEM toolbox (Murphy et al., 2013). It contains functions for the
% assembly of data, the quantum yield calculation, data exploration and
% data postprocessing to enable users to feed the data back into the drEEM
% toolbox.
%
%-----------------------------------------------
% FUNCTIONS
%-----------------------------------------------
% MAIN
% assembleMetadata         Assemble metadata describing your experiment
% aqy                      Calculation of AQYs, cross-calibration as well as sample calculations
% aqyview                  view AQY of samples including errorbars if possible
% epsilon                  Calculate molar fluorescence and absorbance 
% mergeDS                  merge and export aquaDOM datasets in drEEM-compatible structures  
% splitDS                  Split EEM datasets into AQY standards and samples
%-----------------------------------------------
% AUXILLIARY FUNCTIONS
% qy_plot_errBar           plot results within aqy.m
% qy                       calculation of individual qy's
% AbsSNR                   Calculate absorbance noise
% formatCheck              Initial check of dataset for aquaDOM compatibility
%
%-----------------------------------------------
% REFERENCES
%-----------------------------------------------
% 
% drEEM toolbox
% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 
%
%% END of contents.m