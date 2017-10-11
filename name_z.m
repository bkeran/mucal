function z=name_z(name) 
% Given an element name, return its atomic number
% boyan boyanov 2/95
% Adapted for Matlab - Chris Hall 2011
mucalData;
global NELEM element
% copy only first two non-blank characters */
ename=sscanf(name, '%2c');
             
%  convert the name to appropriate format */
  ename(1) = upper(ename(1));
  if exist('ename(2)')
      ename(2) = lower(ename(2));
  end
  
% search element list */  
  for i=1:NELEM
    if (strcmp(ename, element(i)))
        z=i;
        return
    end
  end
% can't find name in list */
  z=0;