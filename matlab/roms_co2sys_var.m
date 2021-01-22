function co2_data = roms_co2sys_var(varname,temp,salt,alkalinity,TIC,press)
% co2_data = roms_co2sys_var(varname,temp,salt,alkalinity,TIC,press)
%
% varname = 'OmegaCa', 'OmeagAr', 'pHtotal' or just 'pH'
% temp,salt,alkalinity,TIC = ROMS BGC model output
% press = pressure in db (negative z_r from grid structure)
%
% John Wilkin - March 2019
%
% Thie function requires CO2SYS.m from
% https://www.mathworks.com/matlabcentral/fileexchange/78378-co2sysv3-for-matlab

% We have to make some assumptions about some parameters that are not
% simulated in ROMS biogeochemical model:
sil = 50;    % Concentration of silicate in the sample (in umol/kg)
po4 =  2;    % Concentration of phosphate in the sample (in umol/kg)
pHscale = 1; % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c = 4;   % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c = 1;   % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% additional inputs for CO2SYS version 3.1
nh4 = 0;
h2s = 0;
kfconstant = 2; % see CO2SYSv3 KF of Perez & Fraga, 1987  (PREFERRED)
boron = 2; % TB of Lee 2010

% Something is messed up with handling NaNs and 999s in version 3
% I'll figure this out some other day
CO2SYSversion = 1;

switch CO2SYSversion
  case 1
    A = CO2SYS(alkalinity(:),TIC(:),1,2,...
      salt(:),temp(:),temp(:),press(:),press(:),...
      sil,po4,...
      pHscale,k1k2c,kso4c);
    iCO2 = [30 31 37]; % version 1.1 outputs order
  case 3
    A = CO2SYSv3(alkalinity(:),TIC(:),1,2,...
      salt(:),temp(:),temp(:),press(:),press(:),...
      sil,po4,...
      nh4,h2s,...
      pHscale,k1k2c,kso4c,...
      kfconstant,boron);
    iCO2 = [34 35 41]; % version 3.1 outputs order
end

%%
OmegaCa = reshape(A(:,iCO2(1)),size(salt));
OmegaAr = reshape(A(:,iCO2(2)),size(salt));
pHtotal = reshape(A(:,iCO2(3)),size(salt));

switch lower(varname)
  case 'omegaca'
    co2_data = OmegaCa;
  case 'omegaar'
    co2_data = OmegaAr;
  case {'phtotal','ph'}
    co2_data = pHtotal;
  otherwise
    warning([varname ' not a CO2 sys variable option'])
end


