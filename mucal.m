%%  * mucal
% Given an element name and a photon energy, calculate a whole slew of
% x-ray data.
%  * This is a program to calculate x-sections using McMaster
%  * data in May 1969 edition.
%  *
%  * NOTE: this program has data for all the elements from 1 to 94
%  * with the following exceptions:
%  *     84  po \
%  *     85  at |
%  *     87  fr |
%  *     88  ra  > Mc Master did not have any data for these
%  *     89  ac |
%  *     91  pa |
%  *     93  np /
%  *
%  * Usage:
%  *
%  * [xsec, errmsg, err]=mucal(name, ZZ, ephot, pflag)
%  *
%  * boyan boyanov 2/95
%  *---------------------------------------------------------------*/
% Adapted for Matlab - Chris Hall 2011
%
function [xsec, errmsg, err]=mucal(name, ZZ, ephot, pflag)
mucalData; % Input the data matricies.
global ZMAX
global element k_edge l1_edge l2_edge l3_edge m_edge
global k_alpha1 k_beta1 l_alpha1 l_beta1 at_weight density
global l1_jump l2_jump l3_jump conv_fac k_yield l_yield
global xsect_coh xsect_ncoh k_fit l_fit m_fit n_fit
Z=0;
xsec=zeros(3,1);
% the return codes for mucal */
No_error = 0;         % no error */
No_input=-1;          % no name, no Z, no service */
No_zmatch=-2;         % Z does not match name */
No_data=-3;           % data not avaialble for requested material */
Bad_z=-4;             % bad Z given as input */
Bad_name=-5;          % invalid element name */
Bad_energy=-6;        % negative or zero photon energy */
Within_edge=-7;       % photon energy within 1 eV of an edge */
M_edge_warn=-8;       % M-edge data for a Z<30 element requested */

stderr=2;
errmsg = [];       % no errors yet */
err = No_error;
nameSz=length(name);
%   % either name or Z must be given */
if (nameSz==0 && ZZ==0)
    errmsg='mucal: no shirt/name, no shoes/Z, no service';
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg); end
    err= No_input;         % this is a terminal error */
    Z=0; xsec=0;
    return
end

%   % ZZ must not be negative */
if (ZZ < 0)
    errmsg='mucal: Z must be non-negative';
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg); end
    err=Bad_z;         % this is a terminal error */
    Z=0; xsec=0;
    return
end

%   % determine material Z, if necessary */
if (nameSz)
    Z = name_z(name);
    if (ZZ>0 && Z ~= ZZ) % Z and name, if both given, must agree */
        errmsg='mucal: Z and element name are not consistent';
        if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg), end; %#ok<*PRTCAL>
        err=No_zmatch;  % this is a terminal error */
        Z=0; xsec=0;
        return
    end
else
    Z = ZZ;
end

% make sure material is available */
if Z==84 || Z==85 || Z==87 || Z==88 || Z==89 || Z==91 || Z==93
    errmsg=...
        'mucal: no data is avaialble for Po, At, Fr, Ra, Ac, Pa, Np';
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg),end;
    err=No_data;         % this is a terminal error */
    Z=0; xsec=0;
    return
end
% Z must be less than ZMAX */
if (Z > ZMAX)
    errmsg=sprintf('mucal: no data for Z>%d', ZMAX);
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg), end;
    err=No_data;  % this is a terminal error */
    Z=0; xsec=0;
    return
end

%   % name must be a valid element symbol */
if (~Z)
    errmsg=sprintf('mucal: invalid element name %s', name);
    if (pflag), fprintf(errmsg, '\n%s\a\n\n', errmsg), end;
    err=Bad_name; % this is a terminal error */
    Z=0; xsec=0;
    return
end

% OK, input is fine */
if (nameSz)
    name = element(Z); %#ok<NASGU>
end
% cannot calculate at negative energies */
if (ephot < 0.0)
    errmsg= 'mucal: photon energy must be non-negative';
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg), end;
    err=Bad_energy; % this is a terminal error */
    Z=0; xsec=0;
    return
end

% stuff the energy-independent parts of all arrays */

xsec(1) = Z;
xsec(2) = at_weight(Z);

% is ephot=0 return physical constants and x-ray energies only */
if (ephot == 0.0)
    err=Bad_energy;
    return
end


%   % determine shell being ionized */
if (ephot >= k_edge(Z))                % K shell */
    shell = 1;
elseif (ephot >= l3_edge(Z))          % L shell */
    shell = 2;
elseif (ephot >= m_edge(Z))           % M1 subshell */
    shell = 3;
else                                % everything else */
    shell = 4;
end

% calculate photo-absorption barns/atom x-section */
switch (shell)
    case 1        % K shell */
        barn_photo = mcmaster(ephot, k_fit(Z,:));
    case 2        % L shell */
        barn_photo = mcmaster(ephot, l_fit(Z,:));
        if (ephot >= l1_edge(Z))   % above L1-no corrections */
        elseif (ephot >= l2_edge(Z)) % between L1 and L2 */
            barn_photo = barn_photo/l1_jump;
        elseif (ephot >= l3_edge(Z))   % between L2 and L3 */
            barn_photo = barn_photo/(l1_jump * l2_jump);
        end
    case 3        % M1 subshell */
        barn_photo = mcmaster(ephot, m_fit(Z,:));
    case 4       % all other shells  */
        barn_photo = mcmaster(ephot, n_fit(Z,:));
    otherwise       % this should never happen */
        errmsg= 'mucal: congratualtions, you have just found a bug';
        if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg);
            err=-9;
            energy=0; xsec=0; fluo=0;
            return
        end
end


% calculate coherent, incoherent x-sections, and total */
barn_coh = mcmaster(ephot, xsect_coh(Z,:));
barn_ncoh = mcmaster(ephot, xsect_ncoh(Z,:));
barn_tot = barn_photo + barn_coh + barn_ncoh;

% stuff the x-section array with the barn/atom data */
xsec(3) = barn_tot * density(Z) / conv_fac(Z);  % absorption  coef */

end