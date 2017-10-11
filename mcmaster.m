% /*--------------------------------------------------------------- 
%  * mcmaster
%  *    given a photon energy and the fit coefficients, calculate 
%  *    a mcmaster x-section
%  *
%  * boyan boyanov 2/95
%  *---------------------------------------------------------------*/ 
% Adapted for Matlab - Chris Hall 2011
%
function mcm=mcmaster(ephot, fit)
   
%   /* ephot = 1 need special handling */
%   /* no it doesn't really! CUS 11/02/2005 */
% /*  log_e = (ephot == 1.0) ? ephot : log(ephot); */
log_e = log(ephot);
xsec=0;
for i=1:4
    xsec=xsec+(fit(i)*log_e^(i-1));
end
mcm=exp(xsec);