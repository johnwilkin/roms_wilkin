function varlistout = roms_varlist(varargin)
% get cell array of subsets of roms variable names
%
% varlist = roms_varlist(categories...) % will be concatenated
%
% categories are: 
%
%   physics, physics2d, physics3d, mixing3d, 
%   s-param, s-coord, grid, 
%   fennel (incl oxygen), fennelN (N only), fennelC (C only)
%   usecos, usecosC (new C variables only), dom (synonym)
%   biodiags (all), fenneldiags (just fennel), usecosdiags (just usecos)
%   bulkflux (all inputs to bulk fluxes)
%   fluxes (stress, net heat flux and net shortwave)
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_varlist.m 526 2019-05-15 18:59:30Z wilkin $

if nargin == 0
  help(mfilename)
end

model = 'roms';

for k=1:nargin
  
  category = varargin{k};
  
  switch model
    
    case 'roms'
      switch category
        case 'time'
          varlist = 'ocean_time';
        case 'physics'
          varlist = {'temp','salt','u','v','zeta','ubar','vbar'};
        case 'physics2d'
          varlist = {'zeta','ubar','vbar'};
        case 'physics3d'
          varlist = {'temp','salt','u','v'};
        case 'mixing3d'
          varlist = {'AKv','AKt','AKs'};
        case 's-param'
          varlist = {'theta_s','theta_b','Tcline','hc','Vtransform','Vstretching'};
        case 's-coord'
          varlist = {'s_rho','s_w','Cs_r','Cs_w'};
        case 'grid'
          varlist = {'h','f','pm','pn','angle','lon_rho','lat_rho',...
            'lon_u','lat_u','lon_v','lat_v','lon_psi','lat_psi',...
            'mask_rho','mask_u','mask_v','mask_psi',...
            'dmde','dndx'};
        case 'fennel'
          varlist = {'NO3','NH4','chlorophyll','phytoplankton','zooplankton',...
            'LdetritusN','SdetritusN','TIC','alkalinity','LdetritusC',...
            'SdetritusC','oxygen'};
        case 'fennelC'
          varlist = {'NO3','NH4','chlorophyll','phytoplankton','zooplankton',...
            'LdetritusN','SdetritusN','TIC','alkalinity','LdetritusC',...
            'SdetritusC'};
        case 'fennelN'
          varlist = {'NO3','NH4','chlorophyll','phytoplankton','zooplankton',...
            'LdetritusN','SdetritusN'};
        case {'dom','usecos','usecosC'}
          varlist = {'semilabileDON','refractoryDON',...
            'semilabileDOC','refractoryDOC'};
        case {'co2sys'}
          varlist = {'OmegaCa','OmegaAr','pH'};
        case 'usecosall'
          varlist = {'NO3','NH4','chlorophyll','phytoplankton',...
            'zooplankton',...
            'LdetritusN','SdetritusN','TIC','alkalinity','LdetritusC',...
            'SdetritusC','oxygen','semilabileDON','refractoryDON',...
            'semilabileDOC','refractoryDOC','OmegaCa','OmegaAr','pH'};
        case 'fenneldiags'
          varlist = {'denitrification','CO2_airsea','pCO2',...
            'O2_airsea','P_Production','NO3_uptake'};
        case 'usecosdiags'
          varlist = {'nitrogen_buried','carbon_bottom','carbon_buried',...
            'C_excess_uptake'};
        case 'biodiags'
          varlist = {'denitrification','CO2_airsea','pCO2',...
            'P_Production','NO3_uptake','nitrogen_buried',...
            'carbon_bottom','carbon_buried','C_excess_uptake'};
        case 'bulkflux'
          varlist = {'Uwind','Vwind','Pair','Tair','Qair','swrad',...
            'lwrad','lwrad_down','rain'};
        case 'fluxes'
          varlist = {'sustr','svstr','shflux','swrad','swflux'};
        otherwise
          error("Input category "+category+" is not an allowed option")
      end
      
    case 'ncom'
      switch category
        case 'physics'
          varlist = {'water_temp','salinity','water_u','water_v','surf_el',[],[]};
        case 'physics2d'
          varlist = {'surf_el',[],[]};
        case 'physics3d'
          varlist = {'water_temp','salinity','water_u','water_v'};
      end
  end
  
  if k==1
    varlistout = varlist;
  else
    varlistout = [varlistout varlist];
  end
  
end

