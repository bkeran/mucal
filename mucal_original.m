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
%  * [energy, xsec, fluo, mu, errmsg, mu] = mucal(element, Z, e_phot,
%  unit, print_flag)
%  *
%  * boyan boyanov 2/95
%  *---------------------------------------------------------------*/
% Adpated for Matlab - Chris Hall 2011
%
function [energy, xsec, fluo, errmsg, err]=mucal(name, ZZ, ephot, unit, pflag)
mucalData; % Input the data matricies.
global ZMAX
global element k_edge l1_edge l2_edge l3_edge m_edge
global k_alpha1 k_beta1 l_alpha1 l_beta1 at_weight density
global l1_jump l2_jump l3_jump conv_fac k_yield l_yield
global xsect_coh xsect_ncoh k_fit l_fit m_fit n_fit
energy=zeros(9,1);
xsec=zeros(11,1);
fluo=zeros(4,1);
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
    energy=0; xsec=0; fluo=0;
    return
end

%   % ZZ must not be negative */
if (ZZ < 0)
    errmsg='mucal: Z must be non-negative';
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg); end
    err=Bad_z;         % this is a terminal error */
    energy=0; xsec=0; fluo=0;
    return
end

%   % determine material Z, if necessary */
if (nameSz)
    Z = name_z(name);
    if (ZZ>0 && Z ~= ZZ) % Z and name, if both given, must agree */
        errmsg='mucal: Z and element name are not consistent';
        if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg), end; %#ok<*PRTCAL>
        err=No_zmatch;  % this is a terminal error */
        energy=0; xsec=0; fluo=0;
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
    energy=0; xsec=0; fluo=0;
    return
end
% Z must be less than ZMAX */
if (Z > ZMAX)
    errmsg=sprintf('mucal: no data for Z>%d', ZMAX);
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg), end;
    err=No_data;  % this is a terminal error */
    energy=0; xsec=0; fluo=0;
    return
end

%   % name must be a valid element symbol */
if (~Z)
    errmsg=sprintf('mucal: invalid element name %s', name);
    if (pflag), fprintf(errmsg, '\n%s\a\n\n', errmsg), end;
    err=Bad_name; % this is a terminal error */
    energy=0; xsec=0; fluo=0;
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
    energy=0; xsec=0; fluo=0;
    return
end

% stuff the energy-independent parts of all arrays */
energy(1) = k_edge(Z);
energy(2) = l1_edge(Z);
energy(3) = l2_edge(Z);
energy(4) = l3_edge(Z);
energy(5) = m_edge(Z);
energy(6) = k_alpha1(Z);
energy(7) = k_beta1(Z);
energy(8) = l_alpha1(Z);
energy(9) = l_beta1(Z);

xsec(7) = at_weight(Z);
xsec(8) = density(Z);
if (Z+1 > 27)
    xsec(9) = l1_jump;
    xsec(10) = l2_jump;
else
    xsec(9) = 0.0;
    xsec(10) = 0.0;
end
xsec(11) = l3_jump(Z);

fluo(1) = k_yield(Z);
fluo(2) = l_yield(Z,1);
fluo(3) = l_yield(Z,2);
fluo(4) = l_yield(Z,3);

% is ephot=0 return physical constants and x-ray energies only */
if (ephot == 0.0)
    err=Bad_energy;
    return
end

%   % check for middle of edge input */
if ( (abs(k_edge(Z) - ephot)  <= 0.001) ||...    % data within K edge */
        (abs(l1_edge(Z) - ephot) <= 0.001) ||... % data within L1 edge */
        (abs(l2_edge(Z) - ephot) <= 0.001) ||... % data within L2 edge */
        (abs(l3_edge(Z) - ephot) <= 0.001) ||...  % data within L3 edge */
        (abs(m_edge(Z) - ephot) <= 0.001) )      % data within M edge */
    errmsg=sprintf('%s\%s',...
        'mucal:  photon energy  is within 1 eV of edge',...
        ' fit results may be inaccurate');
    if (pflag), fprintf(stderr, '\n%s\a\n\n', errmsg);
        err=Within_edge;     % non-terminal error */
    end
end
%
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
            err=Satan_rules;
            energy=0; xsec=0; fluo=0;
            return
        end
end
% M edges for Z<30 are unreliable */
if ((shell > 2) && (Z+1 < 30))
    errmsg=sprintf('%s\n%s',...
        'mucal: McMaster et al. use L-edge fits for the M edges for Z<30',...
        'WARNING: results may be inaccurate');
    if (pflag), fprintf(errmsg, '\n%s\a\n\n', errmsg); end;
    err=M_edge_warn;
end

% calculate coherent, incoherent x-sections, and total */
barn_coh = mcmaster(ephot, xsect_coh(Z,:));
barn_ncoh = mcmaster(ephot, xsect_ncoh(Z,:));
barn_tot = barn_photo + barn_coh + barn_ncoh;

% stuff the x-section array with the barn/atom data */
xsec(1) = barn_photo;     % Photo electric
xsec(2) = barn_coh;       % Coherent
xsec(3) = barn_ncoh;      % Incoherent
xsec(4) = barn_tot;       % Total
xsec(5) = conv_fac(Z);    % conversion factor */
xsec(6) = barn_tot * density(Z) / conv_fac(Z);  % absorption  coef */

% convert to cm^2/g if necessary */
if (upper(unit) == 'C')
    for i=1:4
        xsec(i)= xsec(i)/xsec(5);
    end
    % we are done */
    err=0;
end
