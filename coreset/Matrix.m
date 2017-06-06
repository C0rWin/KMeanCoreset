classdef Matrix < HandleObject
  % A matrix of size n-by-d. Used to store a set P of n points in R^d, or a
  % subspace (derived class).
  properties
    % Used in makeRandomNoise that constructs a random set of points near
    % a random subspace. The subspace is saved in this property for
    % analysis purposes. Type=Cell<Subspace>
    %randSubspace was changed to randSubspaces due to clusters
    %each cluster has a randSubspace of it's own, all saved in this
    %cell array
    randSubspaces
    
    regularSVD=false;
  end % properties
  
  properties(Dependent)
    % Contains a reference to privateMatrix, which is the matrix in
    % a matlab (non-OOP) format.
    % The reason for this wrapper is that in order of a property to be used by derived class,
    % such as LowDimCoreset, it must be depedant.
    %
    % Note: To make Matrix(x,y) call Matrix.matrix(x,y)
    matrix;
    
    m % The same as matrix ( a shorthand)
    nRows % the number of rows in privateMatrix
    nCols % the number of cols in privateMatrix
    n % shorthand for nRows
    d % shorthand for nCols
  end % properties (Dependent)
  
  properties (Access = protected)
    % the j-subspace that minimizes the sum of squared distances to the
    % set of points represented by Matrix. Updated only when
    % getOptSqSubspace is called.
    % j is a propery of optSqSubspace
    optSqSubspace;
    
    %Singular Values of the SVD decomposition of the matrix
    singularVals
    
    % The sum of squard distances from P to optSqSubspace
    optSqDistances;
    
    % The same, for some of squared projections
    optSqProjections;
    
    % see property Matrix
    privateMatrix;
  end % properties (Access = protected)
  
  methods (Access=public)
    function result=isEmpty(obj)
      result= obj.n==0;
    end
    function  [Q]= orth(obj, full)
      % find an orthonormal base subspace=Q for the columns of matrix, similar to orth.m in matlab.
      % The matrix R is the coordinate of of the points on Q, i.e, P=QR.
      % As written in the code of ORTH.M, qr is less accurate. However, it is
      % faster.
      % note: Maybe use qrinsert and stop when the rank is full for faster code.
      % Input:
      %      maxR  (optional)- an upper bound for the rank of A. Show error
      %      message if the rank of the orthogonalization is higher.
      %
      % Output:
      %      Q - n-by-r matrix, where r is the rank of A
      %      R - r-by d matrix such that A=QR
      A=obj.m;
      % use singular value decomposition, SVD, instead of orthogonal
      % factorization, QR.  This doubles the computation time, but
      % provides more reliable and consistent rank determination.
      if issparse(A)
        [U,S] = svds(A,min(size(A)));
        U = sparse(U);
        S = sparse(S);
      else
        [U,S] = svd(A,0);
      end
      [m,n] = size(A);
      if m > 1, s = diag(S);
      elseif m == 1, s = S(1);
      else s = 0;
      end
      tol = max(m,n) * max(s) * eps(class(A));
      if tol<inf
        r = sum(s >= tol);
      else
        r=length(s);
      end
      Q = U(:,1:r);
      Q = Matrix(Q);
      
      
      % Here is the old code.  Use it for faster computation, or for
      % generating the same results as earlier versions of MATLAB.
      %             if ~exist('maxR','var')
      %                 maxR=inf;
      %             end %if
      %             [Qfull,Rfull,~]=qr(A,0); %[Q R] = qr(A,0) and [Q R E] = qr(A,0) produces different results
      %             [r,minError]=obj.rank(Qfull,Rfull);
      %             if (r>maxR)
      %                 r=maxR;
      %                 disp('rank more than desired upper bound..')
      %                 disp(minError);
      %             end %if
      %             Q=Matrix(Qfull(:,1:r)); % remove zero vectors
      %             R=Matrix(Rfull(1:r,:)'); % compute coordination on new base
    end % function orth
    
    % get entries from the matrix as a Matrix object
    function result = Mat(obj,rindexes,cindexes)
      result = Matrix(obj.m(rindexes,cindexes));
    end
    
    %get function for the property matrix
    function mat = getRawMatrix(obj)
      mat = obj.m;
    end
    % Compute the j-subspace that minimizes sum of squared distances
    % to P (rows)
    function result = computeOptSqSubspace(obj, j, regSVD)
      % Compute the j-subspace that minimizes the sum of squared distances to P,
      % over all j-subspaces in R^d.
      %
      % Input:
      %        j - an integer such that 0<j<d
      % Output: optSubspace (Type=Subspace) a subspace that spans the optimal j-subspace of P,
      if ~exist('regSVD', 'var')
        regSVD=obj.regularSVD;
      end
      switch regSVD
        case 2
          %[~, ~, result]=svd(obj.m );
          [tempA, tempB, result]=svd(obj.m );
          result=Subspace(Matrix(result));
        case 1
          %[~, ~, result]=svds(obj.m ,j);
          [tempA, tempB, result]=svds(obj.m ,j);
          result=Matrix(result);
        case 0
          try
            if(obj.d<=obj.n)
              opt = obj.DsmallerThanN(obj.m,j);
            else
              Vj = obj.DsmallerThanN(obj.m',j);
              DVjTrans=obj.m'*Vj;
              DVjTrans = DVjTrans(:,sum(DVjTrans,1)~=0); %eliminate zero columns
              for i=1:size(DVjTrans,2)
                opt(:,i) = DVjTrans(:,i)/norm(DVjTrans(:,i)); %some column may be all zeros
              end
            end
            result = Matrix(sparse(opt));
          catch exception
            throw(exception);
          end % catch
          result=Subspace(result,false);
      end % end
      obj.optSqSubspace = result;
      
      %             [~,~,V,flag] = svds(P,j);
      %             result = Matrix(V);
      %             if(flag==1)
      %                 warning('svds results are unstable');
      %                 try
      %                     [~,~,V2] = svd(full(P),0);
      %                     V2 = V2(:,1:j);
      %                     result = Matrix(V2);
      %                 catch exception
      %                     warning('matrix to large for full svd');
      %                 end
      %             end
      
      
      %This implementation maybe faster:
      %             options.issym=true;
      %             options.isreal=true;
      %
      %             % Here we use the fact that we don't need the singular vectors, only their
      %             % span (by first eigen vectors of P'*P)
      %             if (obj.d<=obj.n)
      %                 [nonOrthogonalOptSubspace, ~]=eigs(P'*P,j,'lm',options);
      %                 % Note: Due to some numerical problems, it seems that optCost1
      %                 % is quite different than the sum of last d-j singular values of P'*P
      %                 % Maybe svd will improve equality?
      %             else
      %                 % Compute a span of the left sinuglar vectors of P, instead
      %                 % of right.
      %                 [Uj,~]=eigs(P*P',j,'lm',options);
      %                 nonOrthogonalOptSubspace=P'*Uj;
      %             end % if
      
      % Here we compute an orthogonal base for the optimal subspace.
      % note: 'orth' can be used instead of QR. The running time will be longer, but
      % results are more accurate. see "orth.m"
      %             [result,~]=Matrix(nonOrthogonalOptSubspace).orth(j);
      %             result=Subspace(result,true); %no need for orth, result is allready orthogonal
      %             obj.optSqSubspace=result;
    end % computeOptSubspace
    
    % Make the matrix an empty matrix ([])
    function clear(obj)
      obj.m=[];
    end % function clear
    
    function r=horzcat(obj1,obj2)
      if isempty(obj1)
        r=Matrix(obj2.m);
      else
        r=Matrix([obj1.m obj2.m]);
      end % if
    end
    function r=vertcat(obj1,obj2)
      r=Matrix([obj1.m; obj2.m]);
    end
    
    function r = plus(obj1,obj2)
      % MTIMES   Implement obj1 * obj2 for Matrix.
      r = Matrix(obj1.m+obj2.m);
    end % plus
    
    function r = minus(obj1,obj2)
      % MTIMES   Implement obj1 * obj2 for Matrix.
      r = Matrix(obj1.m-obj2.m);
    end % minus
    
    % compute epsilon-grid around the point of the Matrix.
    % works only for 1-line matrix.
    % Put the matrix points of G in the grid and choose representatives
    function [Q idxs]=epsilonGrid (obj, epsilon, G)
      if obj.n~=1
        error('only one point matrix allowed');
      end
      H=Matrix(bsxfun(@minus, G.m, obj.m));
      ma=max(max(abs(H.m)));
      eps2 = epsilon*ma;
      %[~,idxs] = H.gridPoints(eps2);
      [temp,idxs] = H.gridPoints(eps2);
      Q = G.getRows(idxs);
    end
    % choose eps-net such that every point has a representative of
    % distance at most epsilon
    function [Q idxs rep]=gridPoints(obj, epsilon)
      if epsilon==0 || obj.n==0
        Q = obj.getRows(1:obj.n);
        idxs = 1:obj.n;
      else
        % todo - remove factor 2 that deals with first cell
        epsilon2=epsilon/(2*sqrt(obj.d));
        grid=fix(obj.m/epsilon2);
        %[~, idxs, rep]=unique(grid,'rows');
        [temp, idxs, rep]=unique(grid,'rows');
        Q = obj.getRows(idxs);
        M=Matrix(obj.m-Q.m(rep,:));
        A=sqrt(max(M.squaredRows));
        if A>epsilon
          error('bad grid');
        end
      end
    end
    %
    %           function [Q idxs centers]=gridPoints(obj, s, isCount)
    %                   D=1;%max(obj.m);
    %
    %
    %                   % convert s to sample
    %                   if ~isCount
    %                        s=s/sqrt(obj.d); % error per direction
    %                        k=ceil(D/(2*s)); % number of cells in each direction
    %                   else
    %                        k=s;
    %                   end
    %                   if k==0 || obj.n==0
    %                     Q = obj.getRows(1:obj.n);
    %                     idxs = 1:obj.n;
    %                     centers=[];
    %                   else
    %                    % find j s.t. there is a point that goes to the jth cell
    %                    % [D*(j-1)/(2i)  D*j/(2i)). cells are between 0 to i-1
    %                    cells=floor(obj.m/(k*D));
    %                    cells(i)=i-1;
    %                    centers = D*(1/(2*k)+cells/(2*k));
    %                    [~, idxs]=unique(cells,'rows');
    %                    Q = obj.getRows(idxs);
    %                  end
    %           end
    
    function m=maxDistance(obj,centers)
      P=PointFunctionSet(obj);
      m=P.maxDistance(centers);
    end
    
    function r = mtimes(obj1,obj2)
      % MTIMES   Implement obj1 * obj2 for Matrix.
      r = Matrix(obj1.m*obj2.m);
    end % mtimes
    
    function r=ctranspose(obj)
      r=Matrix(obj.m');
    end
    
    % wrapper for reshape.m
    % A new object is returned. The existing matrix does not change!
    function result=Reshape(obj,varargin)
      result=Matrix(reshape(obj.m,varargin{:}));
    end
    
    function result=Numel(obj)
      result=numel(obj.m);
    end
    
    function result=sub(obj,varargin)
      result=obj.m(varargin{:});
    end
    
    %         function ans = subsasgn(obj, S, B)
    %             if S(1).type ~= '.'
    %                ans=Matrix(subsasgn(obj.m,S,B.m));
    %             else
    %                 ans = builtin('subsasgn',obj,S,B);
    %             end
    %         end % ans
    
    % Support commands such as Matrix (2:1,3:4).
    % see "help subsrefs".
    %function display(obj)
    %    obj.m
    %end
  end % public methods
  
  methods (Access=protected)
    
    %function [r minError]= rank(~, Q, R)
    function [r minError]= rank(temp, Q, R)
      % (The effective) rank of a matrix Q that returned from QR decomposition, based on
      % R. That is, almost zeroes entries are considered zeroes.
      % stolen from "orth.m"
      tol = 100*eps*norm(Q,'fro');
      absDiagR=abs(diag(R));
      
      % Note: tol might be too small and yield high-dim sets, check on debug
      r = sum(absDiagR> tol);
      r = r(1); % fix for case where R is vector.
      
      % Check that the minimum eigenvalue that is greater than tol
      % is much greater than tol, otherwise tol need to be redefined.
      smallTolIndexes=find(absDiagR> tol);
      if isempty(smallTolIndexes)
        minError=1; % all the values are greater than tol.
        % May occur if tol too big
      else
        minError=min(absDiagR(smallTolIndexes'))/tol;
      end % if isempty
      if (minError<1000)
        %                warning(' It seems that the eigenvalue that is higher than eps, is very similar to value that is smaller than eps')
        %                 minError
      end %if (minError<1000)
    end % function rank
    
    % Check if we already computed the optimal j-subspace for P
    function result=isUpdatedOptSqSubspace(obj, j)
      result=false;
      if ~isempty(obj.optSqSubspace)
        if obj.optSqSubspace.j ==j
          result=true;
        end % if obj.opstSqSubspace
      end %if ~isempty
    end % function isUpdatedOptSqSubspace
    
    % Check if we already computed the sum of squared distances to the
    % opt j-subspace of P
    function result=isUpdatedOptSqDistances(obj, j)
      if ~isempty(obj.optSqDistances) && isUpdatedOptSqSubspace(obj, j)
        result=true;
      else result=false;
      end %if
    end% function isUpdated...
    
    % Same as previous function, for sum of squared projections
    function result=isUpdatedOptSqProjections(obj, j)
      if isUpdatedOptSqSubspace(obj, j) && ~isempty(obj.optSqProjections)
        result=true;
      else result=false;
      end % if
    end % function
    
    %computes optimalSqSubspace for matrices with d<=n
    function result = DsmallerThanN(obj,M,j)
      %            %dense calculation
      %            R = full(qr(M,0)); %maybe unstable
      %            R = triu(full(qr(M,0)));
      %            [~,R] = qr(M,0); %use the line above for faster but less numerical stable
      %            R = full(R);
      if(issparse(M))
        R = full(qr(M,0)); %when M is sparse then qr return triangular matrix
      else
        R = triu(qr(M,0)); %when M is dense qr returns GEQRF output.
      end
      [V,DSquared] = eig(R'*R); %probably more stable than eigs, for speed use eigs
      %sorting and taking j highest magnitude
      Ddiag = diag(DSquared)';
      %[~,indices] = sort(Ddiag,'descend');
      [temp,indices] = sort(Ddiag,'descend');
      
      result = V(:,indices);
      result = result(:,1:j);
      %             sparse calculation
      %             R=qr(M,0);
      %             if(size(R,1) == 1)
      %                R = full(R);
      %                [V,DSquared] = eig(R*R); %probably more stable than eigs, for speed use eigs
      %                %sorting and taking j highest magnitude
      %                Ddiag = DSquared;
      %                result = V;
      %             else
      %                options.issym=true;
      %                options.isreal=true;
      %                [result,DSquared,flag] = eigs(R'*R,j,'lm',options);
      %                if flag==1
      %                    warning('eigs is unstable')
      %                end
      %             end
      result = sparse(result);
    end
    
    
    % Clear all properties.
    % Used for properties that are computed once and then stored in
    % memory. So they need to be cleared when the matrix is changed.
    function clearMatrixInfo(obj)
      obj.optSqSubspace=[];
      obj.optSqProjections=[];
      obj.optSqDistances=[];
    end
    
    % Pad with zeroes the input matrix otherMat
    % or obj's matrix have so that they'l the same number of columns.
    % Used by mergeSparse/Dense
    function matrixSizeMatch(obj,otherMat)
      if(size(obj.m,2) > size(otherMat.m,2))
        otherMat.m(1,size(obj.m,2)) = 0;
      elseif (size(obj.m,2) < size(otherMat.m,2))
        obj.m(1,size(otherMat.m,2)) = 0;
      end %function sizematch
    end % function matrixSizeMatch(obj,otherMat)
    
    
    function check(Q,E)
      j=E.nCols;
      m=Q.nCols;
      Q1=Matrix(Q.m(:,1:j));
      if isDifferent(Q'*Q,Matrix(eye(m,m))) || isDifferent(Q1,E)
        t=Q'*Q;
        max(max(t.m-eye(m,m)))
        max(max(Q1.m-E.m))
        error ('bug');
      end
    end
    
    function result=isDifferent(M1,M2)
      if max(max(M1.m-M2.m))>1e-13%1000*eps
        result=true;
      else result=false;
      end
    end
  end % protected methods
  
  methods
    function result=isreal(obj)
      result=1;
    end
    
    % returns row number ind from the matrix
    function result=getRow(obj, ind)
      result=Matrix(obj.m(ind,:));
    end
    
    function result = getRows(obj,indexes)
      result = Matrix(obj.m(indexes,:));
    end
    
    function show3d(obj)
      params.lineWidth =  3   ;    %Line width                      1
      params.viewAngle =[40 17.5];       %Viewing angle                   [40 17.5]
      params.depthDegree =3;     %3D depth degree                 3
      params.nFrames =1;         %# frames to step when evolving data
      params.bgColor  ='w';        %Background color                'w'
      params.animate  =0;        %Animation on / off              1
      params.nRotations = [0 0 1];      %# of rotations about axes       [1 1 1]
      params.dotsOrLine = '-';      %plot dots or lines              '-'
      params.threeD =1;          %3D on / off                     1
      close all;
      figure;
      plot3AniAnaglyph(obj.m, params);
      pause(0.5);
      figure;
      params.animate  =1;
      params.nFrames =500;
      for i=1:10
        plot3AniAnaglyph(obj.m, params);
        pause(1);
      end
    end
    function [result,indexes] = sampleRows(obj,sampleSize,repetition)
      samplingAlg = RandomSamplingAlg();
      [tempM,indexes] = samplingAlg.compute(obj.m,sampleSize,repetition);
      result =  Matrix(tempM);
    end
    
    function result=get.m(obj)
      result=obj.matrix;
    end % function set.matrix
    function set.m(obj, in_matrix)
      obj.matrix=in_matrix;
    end % function set.matrix
    
    function result=get.matrix(obj)
      result=obj.getMatrix();
    end % function get.matrix
    
    function result=getMatrix(obj)
      result=obj.privateMatrix;
    end % function set.matrix
    
    % Set new matrix, and remove all the
    % pre-computed values for the old matrix
    function set.matrix(obj, in_matrix)
      obj.clearMatrixInfo();
      obj.privateMatrix=in_matrix;
    end % function set.matrix
    function result=get.nRows(obj)
      result=size(obj.matrix,1);
    end
    function result=get.nCols(obj)
      result=size(obj.matrix,2);
    end
    function result=get.n(obj)
      result=size(obj.matrix,1);
    end
    function result=get.d(obj)
      result=size(obj.matrix,2);
    end
    
    %adds a randsubspace to randSubspaces
    %randSubspace was changed to randSubspaces due to clusters
    function [] = addRandSubspace(obj,randSubspace)
      obj.randSubspaces = [obj.randSubspaces {randSubspace}];
    end
    
    % varargin = n,d,j,noise
    % Make a random matrix of size n*d, where all the points are lying
    % on a j-dimensional subspace, with additional Gaussian noise added
    % of magnitude 'noise'
    function makeRandomNoise(obj,varargin)
      [myn, myd, myj, mynoise] = ...
        Utils.complete({100, 2, 1, 0.5},varargin);
      proj=randn(myn,myj);
      %         obj.randSubspace=Subspace(myd,myj);
      obj.addRandSubspace(Subspace(myd,myj));
      obj.matrix=...
        proj*obj.randSubspaces{1}.matrix'+mynoise*randn(myn,myd); %using this function only when there's 1 cluster, u.e 1 randsubspace
    end % function makeRandomNoise
    
    % get subset of columns. Indexes is an array of column numbers.
    % Output is a Matrix
    function result=getCols(obj, indexes)
      result=Matrix(obj.m(:,indexes));
    end
    function makeUniformRandomNoise(obj,varargin)
      [myn, myd] = ...
        Utils.complete({100, 2},varargin);
      obj.matrix=...
        randn(myn,myd);
    end % function makeRandomNoise
    
    function result=isSparse(obj)
      result=or(issparse(obj.m),isempty(obj.m)); %fixed for the case that m= []
    end %function isSparse
    
    % Allocate space and create a new sparse matrix, according to the matrices
    % in varargin. Used before merging varargin with obj's matrix in the constructor.
    function allocBeforeMerge(obj,varargin)
      nnzAll=0; % number of nnz in all the matrices
      dMax=0; % maximum dimension (columns) of an input matrix
      nAll=0; % total number of rows (in the output matrix)
      for i=1:(nargin-1)
        nnzAll=nnzAll+ nnz(varargin{i}.matrix);
        dMax=max(dMax,varargin{i}.d);
        nAll=nAll+varargin{i}.n;
      end % for
      obj.m=spalloc(nAll,dMax,nnzAll);
    end % function allocBeforeMerge(varargin)
    
    % Constructors:
    % Matrix('sparse',n,d,j,nzpr)- create a a sparse random matrix.
    %  (all parameters after 'sparse' are optional). See MakeSparseRandomNoise.
    % Matrix('merge',M1,M2..) - merge the rows of the matrices M1, M2
    %   to obj. Pad with zeros if sizes do not fit.
    % Matrix(n,d,..)=Matrix('dense',n,d,..) - Construct a random dense matrix
    % Matrix (m) construct a Matrix from a matlab matrix m.
    function obj=Matrix(varargin)
      if nargin>0
        if strcmp(varargin{1},'dense')
          obj.makeRandomNoise(varargin(2:end));
          %updating distances if singular values input
        elseif strcmp(varargin{1},'uniform')
          obj.makeUniformRandomNoise(varargin(2:end));
        elseif strcmp(varargin{1},'sparse')
          obj.makeSparseRandomNoise(varargin(2:end));
        elseif strcmp(varargin{1},'merge')
          if varargin{2}.isSparse
            % allocate memory for input matrices
            %                         obj.allocBeforeMerge(varargin{2:end});
          end
          obj.merge(obj,varargin{2:end});
        elseif nargin==1
          obj.matrix= varargin{1};
        else
          obj.makeRandomNoise(varargin);
        end % if strcmp
      end %if nargin
    end % constructor
    
    
    % Compute the optimal j-subspace (LMS) or retrive it
    % if was already computed
    function result = getOptSqSubspace(obj,j)
      if j > min(obj.nRows,obj.nCols)
        j = min(obj.nRows,obj.nCols);
      end
      if ~isUpdatedOptSqSubspace(obj, j)
        computeOptSqSubspace(obj, j);
      end % if
      result=obj.optSqSubspace;
    end % function getOptSqSubspace
    
    function [result,selectedIndices] = getApproxSubspace(varargin)
      %improving the optSqSubspace
      %input: varargin{1} = P (input points), varargin{2} = epsilon,
      %varargin{3} = optSqSubspace
      [result,selectedIndices] = ApproxSubspaceAlg.run1(varargin{2:end});
    end
    
    % Compute squared distances to Optimal Subspace,
    % or retreive if were already computed
    function result = getOptSqDistances(obj, j)
      if isUpdatedOptSqDistances(obj, j)
        result=obj.optSqDistances;
      else
        if(isempty(obj.singularVals))
          subspace= obj.getOptSqSubspace(j);
          obj.optSqDistances= subspace.sqDistances(obj);
          result=obj.optSqDistances;
        else  %The singular values are input
          obj.optSqDistances = Subspace(Matrix([]),j);
          eigVals = obj.singularVals.^2;
          obj.optSqDistances = sum(eigVals(j+1:end));
          result = obj.optSqDistances;
        end
      end
    end % function result=getOpt..
    
    % Compute squared projections on OptSubspace or retrieved.
    function result = getOptSqProjections(obj, j)
      if isUpdatedOptSqProjections(obj, j)
        result=obj.optSqProjections;
      else
        subspace= getOptSqSubspace(obj,j);
        obj.optSqProjections= subspace.sqProjections(obj);
        result=obj.optSqProjections;
      end % if
    end % function
    
    % compute sum of squared of the entries in the matrix.
    % Known as frobinues norm
    function [sumOfSquared] = sumOfSquaredEntries(obj)
      % Returns the sum of squared distances of the entries of the matrix P
      sumOfSquared=(norm(obj.matrix,'fro'))^2;
    end % function sumOfSquaredEntries
    
    % Returns a column vector where the i'th entry is the sum of squared
    % entries in the i'th row
    function [result] = squaredRows(obj)
      result= sum(obj.m.^2,2);
    end
    
    % retun a number which is the sum of all entries
    function [sumX] = sumAll(obj)
      sumX=sum(sum(obj.matrix));
    end
    
    function [translatedP] = translate(obj)
      % Translate the rows of P to their center of mass.
      % The new sum of rows is thus zero.
      % Output:
      %   translatedP - n-by-d matrix
      P=obj.matrix;
      mean=sum(P)/obj.n;
      translatedP=bsxfun(@minus,P,mean);
    end
    
    
    % Merge the rows of the input sparse matrices to this (obj's) matrix.
    % varargin is a list of sparse matrices.
    % If you want pre-allocation of space, construct a new merged
    % matrix; see constructor
    function obj=merge(obj, varargin)
      for i=1:nargin-1
        %                 rows1=[];
        %                 rows2=[];
        %                 cols1=[];
        %                 cols2=[];
        %                 vals1=[];
        %                 vals2=[];
        
        
        % All input matrices should be sparse
        %                 if(or(obj.nCols ~=0,obj.nRows~=0))
        %                     if(varargin{i}.isSparse~=obj.isSparse)
        %                         error('input sparsity and obj should be the same')
        %                     end
        %                 end
        
        if(isempty(varargin{i}.m))
          continue;
        end
        
        %bugs with sparsity in matlab R2008 fixed
        %                 if(~issparse(obj.m))
        %                     obj.m=sparse(obj.m);
        %                 end
        %                 if(~issparse(varargin{i}.m))
        %                     varargin{i}.m = sparse(varargin{i}.m);
        %                 end
        
        % Number of columns should be the same. Otherwise - pad
        % with zeroes
        if(~isempty(obj.m))
          obj.matrixSizeMatch(varargin{i});
        end
        obj.m = [obj.m; varargin{i}.m];
        
        %sparse merge maybe faster:
        %                 [rows1,cols1,vals1] = find(obj.m);
        %                 [rows2,cols2,vals2] = find(varargin{i}.m);
        %
        %                 if obj.nRows == 1
        %                    rows1 = rows1';
        %                    cols1 = cols1';
        %                    vals1 = vals1';
        %                 end
        %                 if varargin{i}.nRows == 1
        %                     rows2 = rows2';
        %                    cols2 = cols2';
        %                    vals2 = vals2';
        %                 end
        %                 rows2 = obj.nRows + rows2;
        %                 try
        %                     tempM = sparse([rows1 ; rows2],[cols1 ; cols2],[vals1 ; vals2],obj.nRows+varargin{i}.nRows,max(obj.nCols,varargin{i}.nCols));
        %                     obj.m = tempM;
        %                 catch exception
        %                     do = 1;
        %                 end
        
        
        
      end % for
    end % merge
    
    function makeSparseRandomNoise(obj,varargin)
      % varargin = n,d,j,nzpr (all optional)
      % defaults = 50, 10^6, 1, 5, respectively.
      % Return a random sparse matrix of n rows, d columns, and
      % approximately nzpr non-zeroes per row. The rank of the matrix is
      % j and all the singular values are equal to 1.
      [n, d, j, nzpr] = ...
        Utils.complete({50, 10^6, 1, 5},varargin);
      obj.matrix=...
        sprand(n,d, nzpr/d, ones(1,j));
    end % function makeSparseRandomNoise
    
    % Type = Subspace
    function [sqDists projP]= sqDistances(obj, subspace)
      [sqDists projP]= subspace.sqDistances(obj);
    end
    %
    %         function subspaces = AllJsubspaces(obj,j)
    %             subspaceIndexes = nchoosek(1:obj.nRows,j);
    %             subspaces = [];
    %             for i=1:size(subspaceIndexes,1)
    %                 tempSubset = (obj.m(subspaceIndexes(i,:),:))';
    %                 tempSubspace = Subspace(Matrix(tempSubset));
    %                 subspaces = [subspaces tempSubspace];
    %             end
    %         end
    
    function subspace = getSpannedSubspace(obj)
      %this function calculates the standard basis of obj.m (i.e v_i will be all
      %zeros except the i'th element iff there's atlist 1 non zero p in P in
      %the i'th dimension
      %obj.m is n*dim
      colSum = sum(obj.m,1);
      indexes = find(colSum ~= 0);
      subspace = sparse(1:length(indexes),indexes,ones(1,length(indexes)),length(indexes),length(colSum));
      
      subspace = subspace';
      subspace = Subspace(Matrix(subspace),false); %wrapping the subspace matrix in subspace class
    end
    
    function []=setSingVals(obj,singVals)
      obj.singularVals = singVals;
    end
    
    function [Q,U,R] = newRotate(Vs,Es)
      %V= Matrix(full(V.m));
      %E=Matrix(full(E.m));
      
      % find base for V and E
      B=Matrix(orth([full(Vs.m) full(Es.m)]));
      V= Matrix(full(B.m'*Vs.m));
      E= Matrix(full(B.m'*Es.m));
      
      Q1=V-E*(E'*V); %E*E' will cause a temp d*d matrix
      Q2=orth(Q1);
      Q=[E Q2];
      
      M1=E-V*(V'*E); %V*V' will cause a temp d*d matrix
      %M2=Subspace(M1);
      M2=orth(M1);
      M=[V M2];
      %               check(M,V);
      U=B*Q'*M*B';
      %              Q=B*Q;
      %             M=B*M;
      %Us=U;
      if size(Q.m,2)~=size(M.m,2)
        error('err');
      end
      
      if nargout >= 3
        if issparse(V.m)
          R = (Matrix(speye(V.nRows,V.nRows)) - Q*Q' ) + Q*U'*Q';
        else
          R = (Matrix(eye(V.nRows,V.nRows)) - Q*Q' ) + Q*U'*Q';
        end
      end
      
    end
    
    %Q&U are the rotation matrices for the method obj.rotate
    %R = I-QQ'+QU'Q'  - the rotation matrix given.
    function [Q,U,R] = getRotateMatrix(V,E)
      %V= Matrix(full(V.m));
      %E=Matrix(full(E.m));
      
      % find base for V and E
      Q1=V-E*(E'*V); %E*E' will cause a temp d*d matrix
      Q1.m=full(Q1.m);
      Q2=orth(Q1);
      Q=[E Q2];
      
      M1=E-V*(V'*E); %V*V' will cause a temp d*d matrix
      %M2=Subspace(M1);
      M1.m=full(M1.m);
      M2=orth(M1);
      M=[V M2];
      %               check(M,V);
      U=Q'*M;
      if size(Q.m,2)~=size(M.m,2)
        error('err');
      end
      
      if nargout >= 3
        if issparse(V.m)
          R = (Matrix(speye(V.nRows,V.nRows)) - Q*Q' ) + Q*U'*Q';
        else
          R = (Matrix(eye(V.nRows,V.nRows)) - Q*Q' ) + Q*U'*Q';
        end
      end
      
    end
    
    function rotated = rotate(obj,Q,U)
      P = obj';
      projCoordinates=Q'*P;
      residual=P-Q*projCoordinates;
      clear P;
      fullProj=Q*U'*projCoordinates;
      clear projCoordinates;
      clear Q;
      clear U;
      rotated=residual+fullProj;
      clear residual;
      clear fullProj;
      rotated = rotated';
    end
    
    function distance = distancesTopoint(obj,c)
      %Compute the distances from the points in obj to a specific
      %point.
      %Input: obj, point c
      %Output: the squared distances from each points in obj to this point c
      DifferenceMatrix = bsxfun(@minus, obj.m, c);
      distance = sum( DifferenceMatrix.^2, 2);
    end
    
    % find indexes of the given input points
    function idxs=indexesFromData(obj, subset)
      a=1; % pointer on the full set
      b=1; % pointer on the small set
      j=1;
      sizeA=obj.n;
      sizeB=subset.n;
      A=obj.m;
      B=subset.m;
      idxs=zeros(sizeB,1);
      while a<=sizeA && b<=sizeB
        if (any(isnan(B(b,:))))
          a = find(any(isnan(A(a:end,:)),2),1,'first') + a - 1;
        else
          while ~isequal(A(a,:),B(b,:)) && a<=sizeA
            a=a+1;
          end
        end
        idxs(j)=a;
        j=j+1;
        b=b+1;
        a=a+1;
      end
      if (a>sizeA && b<=sizeB)  || idxs(sizeB)==0
        b
        error('index not found');
      end
    end % indexesFromData
  end % methods
end % class Matrix
