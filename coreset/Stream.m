classdef Stream < handle
% A general class for maintaining a coreset on-line.
% The input coreset algorithm needs to suppo.mergedCoreset and .computeCoreset
% methods, as in AbstractCoresetAlg. 
% The stream collects the original points, untill we have a set lastLeaf of size
% leafSize. Then, we compute coreset in the first level from this set. If there is already a
% coreset in the same level, we continue to merge coresets recursively.
%
% An item in the streaming tree ( that is implemented using stack) consists
% of two properties: coreset, which is the coreset itself, and level, which
% is the level of the coreset in the tree. There are no two coresets in the
% same level. The coreset of minimum level is on the top of the stack.
    properties
        % don't delete data while constructing tree
        saveTree; % boolean
        
        % saved list of coresets constructured over time.
        % Only if saveTree==true
        coresetsList;
        
        % maximum size of stack
        maxN=0;

        %size
        n=0;
        
        % number of points that have been fed to the stream
        % (useful for figuring out the next "time" value
        numPointsStreamed=0;
        
        % Coreset that support .mergedCoreset and .computeCoreset 
        % (usually of type coresetAlg)
        coresetAlg
        
        %BiCriteriaAlg supports .compute
        %of type AbstractBicriteriaAlg
        BicriteriaAlg
        
        % Size of original points to collect before constructing coreset.
        leafSize;
        
        % stack to implement the coresets tree
        stack=Stack();
        
        % The set of last original points collected, since the last coreset
        % construction. Size at most leafSize.
        %Type = AbstractFunctionSet
        lastLeaf;
        
        %will the coreset be using Bicriteria or not
        usingBiCrit = false;
        
        %will the stream use hard disk for saving the stack or not
        paging = false;
        %the paging directory
        pagePath = pwd;
    end % properties

    methods (Access=protected)
        
        % Returns false if integer 'level' is the same as the level of the
        % last coreset that was inserted to the stack.
        % Used in the recursion to decide whether to merge the current coreset 
        % with a previous one, and go up another level.
        function result=isCorrectLevel(obj, level)
            mystack=obj.stack;
            if mystack.isEmpty()
                result=true;
            elseif mystack.top().level>level 
                result=true;
            elseif  mystack.top().level==level 
                result=false;
            else % mystack.top().level > level 
                error('should not occur');
            end % if
        end % function isCorrectLevel
        
        % Get a coreset and insert it to the streaming tree.
        % Merge with existing coresets untill there are no two coresets in
        % the same level of the tree. 
        % 'coreset' is the output of coresetAlg, that i 
        % accepted by coresetAlg.mergedCoreset
        function addCoreset(obj, coreset)
            obj.addCoreset2List(coreset);
            mystack=obj.stack;
            level=1;
            while not(obj.isCorrectLevel(level))
                stackItem = obj.stack.pop();
                 if obj.paging
                    %saving the stack item and freeing memory
                    pagedStack = obj.stack;
                    save([obj.pagePath '\pagedStack.mat'],'pagedStack');
                    obj.stack = [];
                    clear pagedStack;
                end
                coreset=obj.coresetAlg.mergedCoreset(stackItem.coreset,coreset);
%                 load('coresetlog')
%                 c{end+1} = coreset;
%                 save('coresetlog','c');
                
                obj.addCoreset2List(coreset);
                if obj.paging
                    %loading the stack from HardDisk
                     loadedStruct = load([obj.pagePath '\pagedStack.mat']); %matlab's load returns a structure of the data
                     obj.stack = loadedStruct.pagedStack;
                end
                level=level+1;
            end % while
            newStackItem.level=level;    
            newStackItem.coreset=coreset;
            obj.stack.push(newStackItem);
        end % function addCoreset

        
        % Construct a corest and add it to the tree
        function addLeaf(obj,P)
            if obj.paging
                %saving the stack item and freeing memory
                pagedStack = obj.stack;
                save([obj.pagePath '\pagedStack.mat'],'pagedStack');
                obj.stack = [];
                clear pagedStack;
            end
            % the line that computes coreset
            coreset = obj.coresetAlg.computeCoreset(P);
            
%             if (exist('coresetlog.mat','file')>0)
%                 load('coresetlog')
%                 c{end+1} = coreset;
%             else
%                 c = {coreset};
%             end
%             save('coresetlog','c');
            
            if obj.paging
                %loading the stack from HardDisk
                 loadedStruct = load([obj.pagePath '\pagedStack.mat']); %matlab's load returns a structure of the data
                 obj.stack = loadedStruct.pagedStack;
            end
            obj.addCoreset(coreset);
        end
    end % protected methods
    
    methods     
        function addCoreset2List(obj, coreset)
                if obj.saveTree
                    l=length(obj.coresetsList)+1;
                    obj.coresetsList{l}=coreset;
                end
        end
        
        % return streamingTree, if obj.saveTree is true
        function root=getTree(obj)
            if ~obj.saveTree
                error('You should turn on obj.saveTree');
            end
            root=Node.list2Tree(obj.coresetsList);
        end
        
        function updateN(obj)
            % scan coresets and sum num of points
             allC = obj.getUnifiedCoreset;
             obj.n=allC.n;   
             obj.maxN = max(obj.n, obj.maxN);
             %disp('stream size:');
             %obj.maxN
        end
        
        % to avoid recursion
        function addPointSet(obj, P)
            obj.addPointSetBatches(P);
            %obj.updateN(); 
            obj.numPointsStreamed = obj.numPointsStreamed + P.M.n;
        end
        % Add a set of points to the stream. Add AbstractFunctionSet.
        % If the set is larger than leafSize, we should slive it to several
        % sets and construct a coreset on each set.
        % Type of P AbstractFunctionSet.
        function addPointSetBatches(obj, P)
            for i = 1:obj.leafSize:P.M.n
                idx=i:min(i+obj.leafSize,P.M.n);
                if (isempty(idx))
                    break;
                end
                obj.addLeaf(P.subset(idx));
            end
        end % function addPointSet
        
        function addPointSetBatches_old(obj, P)
            if P.M.n == obj.leafSize
                 obj.addLeaf(P)
            elseif P.M.n > obj.leafSize
                obj.addLeaf(P.subset(1:obj.leafSize));
                obj.addPointSetBatches(P.subset(obj.leafSize+1:P.M.n)); 
                % 'missing' is how much points are needed in the last leaf in
                % order to construct a coreset from it
            else
                 if isempty(obj.lastLeaf)
                     % create a leaf of class P
                     obj.addLeaf(P);
                 else
                     missing=obj.leafSize-obj.lastLeaf.size; %P.n<obj.leafSize
                     if P.M.n < missing
                         obj.lastLeaf.merge(P);
                         % Construct a leaf from the first points in P.
                         % Then continue recursively with the rest of points.
                     else
                         obj.lastLeaf.merge(P.subset(1:missing));
                         obj.addLeaf(obj.lastLeaf);
                         obj.lastLeaf = []; %clear lastLeaf after we've added it to stack
                         if P.M.n > missing 
                              obj.addPointSetBatches(P.subset(missing+1:P.M.n))
                         end % if P.n>missing
                     end % if P.n < missing
                 end
            end % if P.n==obj.leafSize
        end % function addPointSet

        % constructor.
        function obj=Stream(varargin)
            if nargin>0
                [obj.coresetAlg, obj.leafSize]=varargin{1:2};
            end
            if nargin > 3
                obj.usingBiCrit = varargin{3};
                obj.BicriteriaAlg = varargin{4};
            end
            if nargin > 4
                obj.paging = varargin{5};
                obj.pagePath = varargin{6};
            end
            obj.stack=Stack();
            obj.lastLeaf=[];
        end % constructor

        % Unite all the coresets in the tree to a single coreset .
        function allC = getUnifiedCoreset(obj)
            % turn the last set of given points into a coreset
                 s = Stack(obj.stack);
                 isFirst=true;
                while not(s.isEmpty())
                    C=s.pop().coreset;
                    if isFirst
                        isFirst=false;
                        allC=C;
                    else
                        allC.merge(C);
                    end
                end % while
                if isempty(obj.lastLeaf)
                    nLeaf=0;
                else
                    nLeaf = obj.lastLeaf.M.n;
                end
                if nLeaf>0
                    allC.merge(obj.lastLeaf);
                end
                
        end % function allC

        % gets an object of type DataStream and compute a coreset till the
        % end of the stream
        function coreset=computeCoresetFromDataStream(obj, dataStream)
            obj.clear();
            dataStream.init();
            [P, eof]=dataStream.getNext();
            PP=PointFunctionSet(P);
            obj.addPointSet(PP);
            iter=0;
            while not(eof)
                clear P;
                [P, eof]=dataStream.getNext();
                PP=PointFunctionSet(P);
                obj.addPointSet(PP);
                iter=iter+1;
                disp(iter);
            end % while
            coreset=obj.getUnifiedCoreset();
        end % computeCoresetFromDataStream
        % Compute a coreset for a set P of points using the streaming
        % tree.
        function coreset=computeCoreset(obj, P)
            obj.clear();
            obj.addPointSet(P);
            coreset=obj.getUnifiedCoreset();
        end % function computeCoreset
        
        % empty stream
        function clear(obj)
            obj.lastLeaf=[];
            obj.stack.clear();
            obj.coresetsList=[];
        end % function clear
    end % methods

    methods (Static)
        function testKMedianStream(n, d, k, iterations, leafSize)
            
            calg=KMedianCoresetAlg();            
            calg.setParameters(k, 1, 0.1, 250, 1);
            S=Stream(calg, leafSize);

            allP=PointFunctionSet();
            
            for i=1:iterations
                P = PointFunctionSet(Matrix(n, d, d, 0.5) );
                allP.merge(P);
                S.addPointSet(P);
            end
            
            allC=S.getUnifiedCoreset();
            [ct, t, ce, e] = compareCoresetKMeans(allP, allC, calg.k, inf, inf);
            % [errorCostPUsingOptC, errorOptCostC]=calg.isGoodCoreset(allP,allC)
        end
        
        function oldtest(n,d,j,iterations, leafSize)
            a = LowDimOriginCoresetAlg();
            a.j=j;
            a.originSampleSize=10;
            a.lowDimSampleSize=20;
            a.isProjection=false;

            s=Stream();
            s.coresetAlg=a;
            s.leafSize=leafSize;

            allP=Matrix();
            for i=1:iterations
                P =Matrix(n,d,j,0.5);
                allP.merge(P);
                s.addPointSet(P);
            end % for
            allC=s.getUnifiedCoreset();
            [errorCostPUsingOptC, errorOptCostC]=a.isGoodCoreset(allP,allC)
        end

        function test(n,d,j,leafSize, coresetSize)
            P =Matrix(n,d,j,0.5);
            a = LowDimOriginCoresetAlg();
            a.j=j;
            a.originSampleSize=coresetSize;
            a.lowDimSampleSize=coresetSize;
            a.isProjection=false;

            s=Stream(a,leafSize);
            coreset=s.computeCoreset(P);
            [errorCostPUsingOptC, errorOptCostC]=a.isGoodCoreset(P,coreset)
        end % test
                
    end % methods (Static)
end % class Stream
