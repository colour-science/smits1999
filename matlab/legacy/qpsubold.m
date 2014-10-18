function [X,lambda,how]=qpsubold(H,f,A,B,vlb,vub,X,neqcstr,verbosity,caller,ncstr,nvars)
%QPSUBOLD Quadratic programming subproblem. Handles qp and constrained
%   linear least-squares as well as subproblems generated from NLCONST.
%
%   X=QPSUBOLD(H,f,A,b) solves the quadratic programming problem:
%
%            min 0.5*x'Hx + f'x   subject to:  Ax <= b 
%             x    
%

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1998/03/03 00:24:20 $
%   Andy Grace 7-9-90. Mary Ann Branch 9-30-96.

% Define constant strings
NewtonStep = 'Newton';
SteepDescent = 'steepest descent';
Conls = 'conls';
Lp = 'lp';
Qp = 'qp';
Qpsub = 'qpsub';

msg = nargchk(12,12,nargin);
if isempty(verbosity), verbosity = 0; end
if isempty(neqcstr), neqcstr = 0; end

LLS = 0;
if strcmp(caller, Conls)
   LLS = 1;
   [rowH,colH]=size(H);
   nvars = colH;
end
if strcmp(caller, Qpsub)
   normalize = -1;
else
   normalize = 1;
end

simplex_iter = 0;
if  norm(H,'inf')==0 | isempty(H), is_qp=0; else, is_qp=1; end
how = 'ok'; 

if LLS==1
   is_qp=0;
end

normf = 1;
if normalize > 0
   % Check for lp
   if ~is_qp & ~LLS
      normf = norm(f);
      if normf > 0
         f = f./normf;
      end
   end
end

% Handle bounds as linear constraints
lenvlb=length(vlb);
if lenvlb > 0     
   A=[A;-eye(lenvlb,nvars)];
   B=[B;-vlb(:)];
end
lenvub=length(vub);
if lenvub>0
   A=[A;eye(lenvub,nvars)];
   B=[B;vub(:)];
end 
ncstr=ncstr+lenvlb+lenvub;

errcstr = 100*sqrt(eps)*norm(A); 
% Used for determining threshold for whether a direction will violate
% a constraint.
normA = ones(ncstr,1);
if normalize > 0 
   for i=1:ncstr
      n = norm(A(i,:));
      if (n ~= 0)
         A(i,:) = A(i,:)/n;
         B(i) = B(i)/n;
         normA(i,1) = n;
      end
   end
else 
   normA = ones(ncstr,1);
end
errnorm = 0.01*sqrt(eps); 

tolDep = 100*nvars*eps;      
lambda=zeros(ncstr,1);
aix=lambda;
ACTCNT=0;
ACTSET=[];
ACTIND=0;
CIND=1;
eqix = 1:neqcstr; 

%------------EQUALITY CONSTRAINTS---------------------------
Q = zeros(nvars,nvars);
R = []; 
indepInd = 1:ncstr; 

if neqcstr>0
   % call equality constraint solver
   [Q,R,A,B,CIND,X,Z,actlambda,how,...
         ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove]= ...
      eqnsolv(A,B,eqix,neqcstr,ncstr,nvars,LLS,H,X,f,normf,normA,verbosity, ...
      aix,how);   
   
   if ~isempty(remove)
      indepInd(remove)=[];
      normA = normA(indepInd);
   end
   
   if ACTCNT >= nvars - 1  
      simplex_iter = 1; 
   end
   [m,n]=size(ACTSET);
   
   if strcmp(how,'infeasible')
      % Equalities are inconsistent, so X and lambda have no valid values
      % Return original X and zeros for lambda.
      return
   end
   
   err = 0;
   if neqcstr > nvars
      err = max(abs(A(eqix,:)*X-B(eqix)));
      if (err > 1e-8)  % Equalities not met
         how='infeasible';
         if verbosity > -1
            disp('Warning: The equality constraints are overly stringent;')
            disp('         there is no feasible solution.')
         end
         % Equalities are inconsistent, X and lambda have no valid values
         % Return original X and zeros for lambda.
         return
      else % Check inequalities
         if (max(A*X-B) > 1e-8)
            how = 'infeasible';
            if verbosity > -1
               disp('Warning: The constraints or bounds are overly stringent;')
               disp('         there is no feasible solution.')
               disp('         Equality constraints have been met.')
            end
         end
      end
      if is_qp
         actlambda = -R\(Q'*(H*X+f));
      elseif LLS
         actlambda = -R\(Q'*(H'*(H*X-f)));
      else
         actlambda = -R\(Q'*f);
      end
      lambda(indepInd(eqix)) = normf * (actlambda ./normA(eqix));
      return
   end
   if isempty(Z)
      if is_qp
         actlambda = -R\(Q'*(H*X+f));
      elseif LLS
         actlambda = -R\(Q'*(H'*(H*X-f)));
      else
         actlambda = -R\(Q'*f);
      end
      lambda(indepInd(eqix)) = normf * (actlambda./normA(eqix));
      if (max(A*X-B) > 1e-8)
         how = 'infeasible';
         if verbosity > -1
            disp('Warning: The constraints or bounds are overly stringent;')
            disp('         there is no feasible solution.')
            disp('         Equality constraints have been met.')
         end
      end
      return
   end
   
   
   % Check whether in Phase 1 of feasibility point finding. 
   if (verbosity == -2)
      cstr = A*X-B; 
      mc=max(cstr(neqcstr+1:ncstr));
      if (mc > 0)
         X(nvars) = mc + 1;
      end
   end
else
   Z=1;
end

% Find Initial Feasible Solution
cstr = A*X-B;
mc=max(cstr(neqcstr+1:ncstr));
if mc>eps
   A2=[[A;zeros(1,nvars)],[zeros(neqcstr,1);-ones(ncstr+1-neqcstr,1)]];
   [XS,lambdas] = qpsubold([],[zeros(nvars,1);1],A2,[B;1e-5], ...
      [],[],[X;mc+1],neqcstr,-2,Qpsub,size(A2,1),nvars+1);
   
   X=XS(1:nvars);
   cstr=A*X-B;
   if XS(nvars+1)>eps 
      if XS(nvars+1)>1e-8 
         how='infeasible';
         if verbosity > -1
            disp('Warning: The constraints are overly stringent;')
            disp('         there is no feasible solution.')
         end
      else
         how = 'overly constrained';
      end
      lambda(indepInd) = normf * (lambdas((1:ncstr)')./normA);
      return
   end
end

if (is_qp)
   gf=H*X+f;
   %  SD=-Z*((Z'*H*Z)\(Z'*gf));
   [SD, dirType] = compdir(Z,H,gf,nvars,f);
   
   % Check for -ve definite problems:
   %  if SD'*gf>0, is_qp = 0; SD=-SD; end
elseif (LLS)
   HXf=H*X-f;
   gf=H'*(HXf);
   HZ= H*Z;
   [mm,nn]=size(HZ);
   if mm >= nn
      %   SD =-Z*((HZ'*HZ)\(Z'*gf));
      [QHZ, RHZ] =  qr(HZ);
      Pd = QHZ'*HXf;
      % Now need to check which is dependent
      if min(size(RHZ))==1 % Make sure RHZ isn't a vector
         depInd = find( abs(RHZ(1,1)) < tolDep);
      else
         depInd = find( abs(diag(RHZ)) < tolDep );
      end  
   end
   if mm >= nn & isempty(depInd) % Newton step
      SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
      dirType = NewtonStep;
   else % steepest descent direction
      SD = -Z*(Z'*gf);
      dirType = SteepDescent;
   end
else % lp
   gf = f;
   SD=-Z*Z'*gf;
   dirType = SteepDescent; 
   if norm(SD) < 1e-10 & neqcstr
      % This happens when equality constraint is perpendicular
      % to objective function f.x.
      actlambda = -R\(Q'*(gf));
      lambda(indepInd(eqix)) = normf * (actlambda ./ normA(eqix));
      return;
   end
end

oldind = 0; 

% The maximum number of iterations for a simplex type method is:
% maxiters = prod(1:ncstr)/(prod(1:nvars)*prod(1:max(1,ncstr-nvars)));

%--------------Main Routine-------------------
while 1
   % Find distance we can move in search direction SD before a 
   % constraint is violated.
   % Gradient with respect to search direction.
   GSD=A*SD;
   
   % Note: we consider only constraints whose gradients are greater
   % than some threshold. If we considered all gradients greater than 
   % zero then it might be possible to add a constraint which would lead to
   % a singular (rank deficient) working set. The gradient (GSD) of such
   % a constraint in the direction of search would be very close to zero.
   indf = find((GSD > errnorm * norm(SD))  &  ~aix);
   
   if isempty(indf) % No constraints to hit
      STEPMIN=1e16;
      dist=[]; ind2=[]; ind=[];
   else % Find distance to the nearest constraint
      dist = abs(cstr(indf)./GSD(indf));
      [STEPMIN,ind2] =  min(dist);
      ind2 = find(dist == STEPMIN);
      % Bland's rule for anti-cycling: if there is more than one 
      % blocking constraint then add the one with the smallest index.
      ind=indf(min(ind2));
      % Non-cycling rule:
      % ind = indf(ind2(1));
   end
   %-----Update X-------------
   
   % Assume we do not delete a constraint
   delete_constr = 0;   
   
   if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
      if strcmp(dirType, NewtonStep)
         % Newton step and hit a constraint: LLS or is_qp
         if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
            STEPMIN = 1;
            delete_constr = 1;
         end
         X = X+STEPMIN*SD;
      else
         % Not a Newton step and hit a constraint: is_qp or LLS or maybe lp
         X = X+STEPMIN*SD;          
      end              
   else %  isempty(indf) | ~isfinite(STEPMIN)
      % did not hit a constraint
      if strcmp(dirType, NewtonStep)
         % Newton step and no constraint hit: LLS or maybe is_qp
         STEPMIN = 1;   % Exact distance to the solution. Now delete constr.
         X = X + SD;
         delete_constr = 1;
      else % Not a Newton step: is_qp or lp or LLS
         if is_qp
            % Is it semi-def, neg-def or indef?
            eigoptions.disp = 0;
            ZHZ = Z'*H*Z;
            [ev,smallRealEig,flag] = eigs(ZHZ,1,'sr',eigoptions);
            if flag  % Call to eigs failed
               [VV,DD] = eig(ZHZ);
               [smallRealEig, eigind] = min(diag(DD));
               ev = VV(eigind(1),:);
            end
         else % define smallRealEig for LLS
            smallRealEig=0;
         end
         
         if (~is_qp & ~LLS) | (smallRealEig < -eps) % LP or neg def: not LLS
            % neg def -- unbounded
            if norm(SD) > errnorm
               if normalize < 0
                  STEPMIN=abs((X(nvars)+1e-5)/(SD(nvars)+eps));
               else 
                  STEPMIN = 1e16;
               end
               X=X+STEPMIN*SD;
               how='unbounded'; 
            else % norm(SD) <= errnorm
               how = 'ill posed';
            end
            if verbosity > -1
               if norm(SD) > errnorm
                  disp('Warning: The solution is unbounded and at infinity;')
                  disp('         the constraints are not restrictive enough.') 
               else
                  disp('Warning: The search direction is close to zero; ')
                  disp('      the problem is ill-posed.')
                  disp('      The gradient of the objective function may be zero')
                  disp('         or the problem may be badly conditioned.')
               end
            end % if verbosity > -1
            return
         else % singular: solve compatible system for a solution: is_qp or LLS
            if is_qp
               projH = Z'*H*Z; 
               Zgf = Z'*gf;
               projSD = pinv(projH)*(-Zgf);
            else % LLS
               projH = HZ'*HZ; 
               Zgf = Z'*gf;
               projSD = pinv(projH)*(-Zgf);
            end
            
            % Check if compatible
            if norm(projH*projSD+Zgf) > 10*eps*(norm(projH) + norm(Zgf))
               % system is incompatible --> it's a "chute": use SD from compdir
               % unbounded in SD direction
               if norm(SD) > errnorm
                  if normalize < 0
                     STEPMIN=abs((X(nvars)+1e-5)/(SD(nvars)+eps));
                  else 
                     STEPMIN = 1e16;
                  end
                  X=X+STEPMIN*SD;
                  how='unbounded'; 
               else % norm(SD) <= errnorm
                  how = 'ill posed';
               end
               if verbosity > -1
                  if norm(SD) > errnorm
                     disp('Warning: The solution is unbounded and at infinity;')
                     disp('         the constraints are not restrictive enough.') 
                  else
                     disp('Warning: The search direction is close to zero; ')
                     disp('      the problem is ill-posed.')
                     disp('      The gradient of the objective function may be zero')
                     disp('         or the problem may be badly conditioned.')
                  end
               end % if verbosity > -1
               return
            else % Convex -- move to the minimum (compatible system)
               SD = Z*projSD;
               dirType = 'singular';
               % First check if constraint is violated.
               GSD=A*SD;
               indf = find((GSD > errnorm * norm(SD))  &  ~aix);
               if isempty(indf) % No constraints to hit
                  STEPMIN=1;
                  delete_constr = 1;
                  dist=[]; ind2=[]; ind=[];
               else % Find distance to the nearest constraint
                  dist = abs(cstr(indf)./GSD(indf));
                  [STEPMIN,ind2] =  min(dist);
                  ind2 = find(dist == STEPMIN);
                  % Bland's rule for anti-cycling: if there is more than one 
                  % blocking constraint then add the one with the smallest index.
                  ind=indf(min(ind2));
               end
               if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                  STEPMIN = 1;
                  delete_constr = 1;
               end
               X = X + STEPMIN*SD; 
            end
         end % if ~is_qp | smallRealEig < -eps
      end % if strcmp(dirType, NewtonStep)
   end % if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
   
   if delete_constr
      % Note: only reach here if a minimum in the current subspace found
      if ACTCNT>0
         if ACTCNT>=nvars-1, 
            % Avoid case when CIND is greater than ACTCNT
            if CIND <= ACTCNT
               ACTSET(CIND,:)=[];
               ACTIND(CIND)=[]; 
            end
         end
         if is_qp
            rlambda = -R\(Q'*(H*X+f));
         elseif LLS
            rlambda = -R\(Q'*(H'*(H*X-f)));
            % else: lp does not reach this point
         end
         actlambda = rlambda;
         actlambda(eqix) = abs(rlambda(eqix));
         indlam = find(actlambda < 0);
         if (~length(indlam)) 
            lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
            return
         end
         % Remove constraint
         lind = find(ACTIND == min(ACTIND(indlam)));
         lind=lind(1);
         ACTSET(lind,:) = [];
         aix(ACTIND(lind)) = 0;
         [Q,R]=qrdelete(Q,R,lind);
         ACTIND(lind) = [];
         ACTCNT = ACTCNT - 2;
         simplex_iter = 0;
         ind = 0;
      else % ACTCNT == 0
         return
      end
      delete_constr = 0;
   end
   
   % Calculate gradient w.r.t objective at this point
   if is_qp
      gf=H*X+f;
   elseif LLS % LLS
      gf=H'*(H*X-f);
      % else gf=f still true.
   end
   
   
   % Update X and calculate constraints
   cstr = A*X-B;
   cstr(eqix) = abs(cstr(eqix));
   % Check no constraint is violated
   if normalize < 0 
      if X(nvars,1) < eps
         return;
      end
   end
   
   if max(cstr) > 1e5 * errnorm
      if max(cstr) > norm(X) * errnorm 
         if verbosity > -1
            disp('Warning: The problem is badly conditioned;')
            disp('         the solution is not reliable') 
            verbosity = -1;
         end
         how='unreliable'; 
         if 0
            X=X-STEPMIN*SD;
            return
         end
      end
   end
   
   if ind % Hit a constraint
      aix(ind)=1;
      ACTSET(CIND,:)=A(ind,:);
      ACTIND(CIND)=ind;
      [m,n]=size(ACTSET);
      [Q,R] = qrinsert(Q,R,CIND,A(ind,:)');
   end
   if oldind 
      aix(oldind) = 0; 
   end
   if ~simplex_iter
      % Z = null(ACTSET);
      [m,n]=size(ACTSET);
      Z = Q(:,m+1:n);
      ACTCNT=ACTCNT+1;
      if ACTCNT == nvars - 1, simplex_iter = 1; end
      CIND=ACTCNT+1;
      oldind = 0; 
   else
      rlambda = -R\(Q'*gf);
      
      if isinf(rlambda(1)) & rlambda(1) < 0 
         fprintf('         Working set is singular; results may still be reliable.\n');
         [m,n] = size(ACTSET);
         rlambda = -(ACTSET + sqrt(eps)*randn(m,n))'\gf;
      end
      actlambda = rlambda;
      actlambda(eqix)=abs(actlambda(eqix));
      indlam = find(actlambda<0);
      if length(indlam)
         if STEPMIN > errnorm
            % If there is no chance of cycling then pick the constraint 
            % which causes the biggest reduction in the cost function. 
            % i.e the constraint with the most negative Lagrangian 
            % multiplier. Since the constraints are normalized this may 
            % result in less iterations.
            [minl,CIND] = min(actlambda);
         else
            % Bland's rule for anti-cycling: if there is more than one 
            % negative Lagrangian multiplier then delete the constraint
            % with the smallest index in the active set.
            CIND = find(ACTIND == min(ACTIND(indlam)));
         end
         
         [Q,R]=qrdelete(Q,R,CIND);
         Z = Q(:,nvars);
         oldind = ACTIND(CIND);
      else
         lambda(indepInd(ACTIND))= normf * (rlambda./normA(ACTIND));
         return
      end
   end %if ACTCNT<nvars
   if (is_qp)
      Zgf = Z'*gf; 
      if ~isempty(Zgf) & (norm(Zgf) < 1e-15) 
         SD = zeros(nvars,1); 
      else
         [SD, dirType] = compdir(Z,H,gf,nvars,f);
      end
   elseif (LLS)
      Zgf = Z'*gf;
      HZ = H*Z;
      if (norm(Zgf) < 1e-15)
         SD = zeros(nvars,1);
      else
         HXf=H*X-f;
         gf=H'*(HXf);
         [mm,nn]=size(HZ);
         if mm >= nn
            [QHZ, RHZ] =  qr(HZ);
            Pd = QHZ'*HXf;
            % SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
            % Now need to check which is dependent
            if min(size(RHZ))==1 % Make sure RHZ isn't a vector
               depInd = find( abs(RHZ(1,1)) < tolDep);
            else
               depInd = find( abs(diag(RHZ)) < tolDep );
            end  
         end
         if mm >= nn & isempty(depInd) % Newton step
            SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
            dirType = NewtonStep;
         else % steepest descent direction
            SD = -Z*(Z'*gf);
            dirType = SteepDescent;
         end
      end
   else % LP
      if ~simplex_iter
         SD = -Z*(Z'*gf);
         gradsd = norm(SD);
      else
         gradsd = Z'*gf;
         if  gradsd > 0
            SD = -Z;
         else
            SD = Z;
         end
      end
      if abs(gradsd) < 1e-10  % Search direction null
         % Check whether any constraints can be deleted from active set.
         % rlambda = -ACTSET'\gf;
         if ~oldind
            rlambda = -R\(Q'*gf);
         end
         actlambda = rlambda;
         actlambda(1:neqcstr) = abs(actlambda(1:neqcstr));
         indlam = find(actlambda < errnorm);
         lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
         if ~length(indlam)
            return
         end
         cindmax = length(indlam);
         cindcnt = 0;
         newactcnt = 0;
         while (abs(gradsd) < 1e-10) & (cindcnt < cindmax)
            cindcnt = cindcnt + 1;
            if oldind
               % Put back constraint which we deleted
               [Q,R] = qrinsert(Q,R,CIND,A(oldind,:)');
            else
               simplex_iter = 0;
               if ~newactcnt
                  newactcnt = ACTCNT - 1;
               end
            end
            CIND = indlam(cindcnt);
            oldind = ACTIND(CIND);
            
            [Q,R]=qrdelete(Q,R,CIND);
            [m,n]=size(ACTSET);
            Z = Q(:,m:n);
            
            if m ~= nvars
               SD = -Z*Z'*gf;
               gradsd = norm(SD);
            else
               gradsd = Z'*gf;
               if  gradsd > 0
                  SD = -Z;
               else
                  SD = Z;
               end
            end
         end
         if abs(gradsd) < 1e-10  % Search direction still null
            return;
         end
         lambda = zeros(ncstr,1);
         if newactcnt 
            ACTCNT = newactcnt;
         end
      end
   end
   
   if simplex_iter & oldind
      % Avoid case when CIND is greater than ACTCNT
      if CIND <= ACTCNT
         ACTIND(CIND)=[];
         ACTSET(CIND,:)=[];
         CIND = nvars;
      end
   end 
end % while 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Q,R,A,B,CIND,X,Z,actlambda,how,...
      ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove]= ...
   eqnsolv(A,B,eqix,neqcstr,ncstr,nvars,LLS,H,X,f,normf,normA,verbosity, ...
   aix,how)
% EQNSOLV Helper function for QPSUB.
%    Finds a feasible point with respect to the equality constraints.
%    If the equalities are dependent but not consistent, warning
%    messages are given. If the equalities are dependent but consistent, 
%    the redundant constraints are removed and the corresponding variables 
%    adjusted.

% set tolerances
tolDep = 100*nvars*eps;      
tolCons = 1e-10;

actlambda = [];
aix(eqix)=ones(neqcstr,1);
ACTSET=A(eqix,:);
ACTIND=eqix;
ACTCNT=neqcstr;
CIND=neqcstr+1;
Z=[]; Anew=[]; Bnew=[]; remove =[];

% See if the equalities form a consistent system:
%   QR factorization of A
[Qa,Ra,Ea]=qr(A(eqix,:));
% Now need to check which is dependent
if min(size(Ra))==1 % Make sure Ra isn't a vector
   depInd = find( abs(Ra(1,1)) < tolDep);
else
   depInd = find( abs(diag(Ra)) < tolDep );
end
if neqcstr > nvars
   depInd = [depInd; ((nvars+1):neqcstr)'];
end      

if ~isempty(depInd)
   if verbosity > -1
      disp('The equality constraints are dependent.')
   end
   how='dependent';
   bdepInd =  abs(Qa(:,depInd)'*B(eqix)) >= tolDep ;
   
   if any( bdepInd ) % Not consistent
      how='infeasible';   
      if verbosity > -1
         disp('The system of equality constraints is not consistent.');
         if ncstr > neqcstr
            disp('The inequality constraints may or may not be satisfied.');
         end
         disp('  There is no feasible solution.')
      end
   else % the equality constraints are consistent
      numDepend = nnz(depInd);
      % delete the redundant constraints:
      % By QR factoring the transpose, we see which columns of A'
      %   (rows of A) move to the end
      [Qat,Rat,Eat]=qr(ACTSET');        
      [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
      remove = i(depInd);
      if verbosity > -1
         disp('The system of equality constraints is consistent. Removing');
         disp('the following dependent constraints before continuing:');
         disp(remove)
      end
      A(eqix(remove),:)=[];
      B(eqix(remove))=[];
      neqcstr = neqcstr - nnz(remove);
      ncstr = ncstr - nnz(remove);
      eqix = 1:neqcstr;
      aix=[ones(neqcstr,1); zeros(ncstr-neqcstr,1)];
      ACTIND = eqix;
      ACTSET=A(eqix,:);
      
      CIND = neqcstr+1;
      ACTCNT = neqcstr;
   end % consistency check
end % dependency check


if ~strcmp(how,'infeasible')
   % Find a feasible point
   if max(abs(A(eqix,:)*X-B(eqix))) > tolCons
      X = A(eqix,:)\B(eqix);  
   end
end

[Q,R]=qr(ACTSET');
Z = Q(:,neqcstr+1:nvars);

% End of eqnsolv.m




