function [allfcns,msg] = prefcnchk(funstr,caller,lenVarIn,gradflag)
%PREFCNCHK Pre- and post-process function expression for FUNCHK.
%   [ALLFCNS,MSG] = PREFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   Since for ATTGOAL, MINIMAX, and CONSTR, the expression may have 
%   two parts (the objective and constraint function), a cell array 
%   is returned in ALLFCNS.  GRADFLAG specifies what strings to look for
%   when expressions are passed in.  
%
%   ALLFCNS is a cell array: the first cell contains the objective function, 
%   the next contains the constraint function, and the last contains a flag 
%   that says if the objective and constraints are together in one function 
%   (strtype==1) or in two functions (strtype==2).  The fourth cell contains
%   the string CALLER.
%
%   Note: ATTGOAL does not have constraints and SEMINF does not support 
%   the expression syntax.  

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 1997/11/29 01:23:22 $
%   Mary Ann Branch 9-1-96

% Initialize
msg='';
allfcns = {};
funfcn = [];
confcn = [];
strtype = 1;

if gradflag
  fdelim='df=';
  gdelim=';dg=';
else
  fdelim='f=';
  gdelim=';g=';
end

switch caller
  case {'attgoal', 'goal'}
    % Try to inline: no need to search for "f=" type expression
    [allfcns{1}, msg] = fcnchk(funstr,lenVarIn);
    if ~isempty(msg)  % return error msg
      return
    end
    allfcns{3} = strtype;
    allfcns{4} = caller;
    allfcns{2} = confcn;
  case {'constr', 'minimax'}
    % Check for "f=" type expression
    ii = find(~isspace(funstr));
    temp = findstr(funstr(ii),fdelim);
    if isempty(temp) 
      % Not type "f=", so must be function name.  
      % But fcnchk should only return a string, not an inline object.
      [funfcn, msg] = fcnchk(funstr,lenVarIn);
      if ~isempty(msg)  % return error msg
        return
      end
      confcn = [];
      if isa(funfcn,'inline')
        msg = ('Objective/constraint expression not properly defined.');
        return
      end
    else % "f=" type expression
      % Find and remove fdelim and gdelim
      strtype = 2;
      startf = ii(temp(1)+length(fdelim));
      temp = findstr(funstr(ii),gdelim);
      if isempty(temp)
        msg = ('Objective/constraint expression not properly defined.');
        return
      end
      endf = ii(temp(1)-1);
      startg = ii(temp(1)+length(gdelim));
      ffunstr = funstr(startf:endf);
      gfunstr = funstr(startg:end);
      % get rid of last semi-colon
      for j = length(gfunstr):-1:1
        if strncmp(gfunstr(j), ';',1)
          break
        elseif ~isspace(gfunstr(j))
          msg = ('Constraint expression must end with a semi-colon');
          return
        end
      end
      gfunstr = gfunstr(1:max(1,j-1));
      % fcnchk may return a string or an inline object
      [funfcn, msg]  = fcnchk(ffunstr,lenVarIn);
      if ~isempty(msg)  % return error msg
        return
      end
      [confcn, msg]  = fcnchk(gfunstr,lenVarIn);
      if ~isempty(msg)  % return error msg
        return
      end
    end % isempty(temp)     
    allfcns{1} = funfcn;
    allfcns{3} = strtype;
    allfcns{4} = caller;
    if strtype == 1
      allfcns{2}= [];
    elseif strtype == 2
      allfcns{2} = confcn;
    end
  case {'seminf'}
    ii = find(~isspace(funstr));
    temp = findstr(funstr(ii),'f=');
    if isempty(temp) 
      % Not type "f=", so must be function name.  
      % But fcnchk should only return a string, not an inline object
      % (seminf does not take expressions, only m-file names).
      [funfcn, msg] = fcnchk(funstr,lenVarIn);
      if ~isempty(msg)  % return error msg
        return
      end
      confcn = [];
      if isa(funfcn,'inline')
        msg = ('Expression syntax not supported for seminf.');
        return
      end
    else
      msg = ('Expression syntax not supported for seminf.');
      return
    end
    allfcns{1} = funfcn;
    allfcns{4} = caller;
    allfcns{3} = strtype;
    allfcns{2} = confcn;
  otherwise
    disp('unknown caller')
end




