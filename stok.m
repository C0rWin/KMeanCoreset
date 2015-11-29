function [istart, iend] = stok(str, delim)
%******************************************************************************
% function stok
% Purpose: Find the addresses of tokens in the strings
% Input:  
%    str    - character string to be searched for tokens
%    delim  - character string holding the token delimiters
%    maxtok - maximum allowed number of tokens
% Output: 
%    strtok - return the number of tokens in the string
%    istart - integer array holding the token starting positions in str
%    iend   - integer array holding the token ending positions in str
% Usage example:
% str='ab;cd de=0'; delim='; ='; [istart, iend] = stok(str, '; =')
% istart =
%     1     4     7    10
% iend =
%     2     5     8    10
% str(istart(1):iend(1))
% ans =
% ab
% Revision history:
%    15/06/2005 created (Michael Yoshpe)
% Remarks:
%******************************************************************************
ns = length(str);
nd = length(delim);
ntok  = 0;
intok  = false;
intokp = false;

% memory preallocation
istart = zeros(1,ns);
iend   = istart;

% check all characters in str
for i=1:ns,
   % check delimiters until the string character is equal to 
   % delimiter, or until all delimiters are checked
   for j=1:nd,
      if(str(i)==delim(j)) 
         intok = false;
         break;
      else
         % we are inside the token if the current string's character is 
         % not equal to any of the delimiters
         intok = true;
      end
   end
   
   %  here the token start and end is identified
   if((intokp==false)&&(intok==true))
      ntok = ntok + 1;
      istart(ntok) = i;
   elseif((intokp==true)&&(intok==false))
      iend(ntok) = i - 1;
   end

   if((i==ns) && (intok==true))
      iend(ntok) = i;
   end
   intokp = intok;
end

istart(ntok+1:end) = [];
iend(ntok+1:end) = [];


