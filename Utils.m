classdef Utils
% statuc global functions that are used by all classes.
    properties (Constant)
        inputColor ='b';
    end
    
    methods (Access=public, Static);
        function tol=findTol(k,foo, initTol, minTol, maxTol)
            if ~exist('minTol','var')
                minTol=0;
            end
            if ~exist('maxTol','var')
                maxTol=inf;
            end
            if ~exist('initTol', 'var')
                tol=[];
            else
                tol=initTol;
            end
            mink = 0; maxk = inf;
            
            while (maxk - mink > 0)
                [l, ~] = foo(tol);
                %l=length(lat1)-1;
                
                if l<k
                    mink = l;
                    maxTol=tol;
                    tol=minTol+(tol-minTol)/2;
                elseif  l>k
                    minTol=tol;
                    maxk = l;
                    if maxTol<inf
                        tol=tol+(maxTol-tol)/2;
                    else
                        % "binary" search on unbounded size
                        tol=minTol*2;
                    end
                else
                    return;
                end % if
            if l<k+1 && maxTol-minTol<eps% dead-lock
                return;
            end
                
            end
        end
        
        function [var value idx]=readParam(var, name, default)
            f=@(x) strcmp(x,name);
            idx=find(cellfun(f,var));
            if isempty(idx)
                value=default;
            else
                value=var{idx+1};
                var([idx idx+1])=[];
            end
        end
        function c=color( i)
               colors=setdiff('rgbcmyk', Utils.inputColor);
               n=length(colors);
               c=colors(mod(i,n)+1);
        end
        
        function varargin=show(varargin)
            [varargin title]= Utils.readParam(varargin, 'title', []);
            [varargin func]= Utils.readParam(varargin, 'func', '');
            [varargin animate]= Utils.readParam(varargin, 'animate', '');
            [varargin pau]= Utils.readParam(varargin, 'pause', false);
            [varargin d3d]= Utils.readParam(varargin, 'd3d', false);
            [varargin mclear]= Utils.readParam(varargin, 'clear', false);
            [varargin P]= Utils.readParam(varargin, 'P',[]);
            if isa(P,'PointFunctionSet')
                P=P.M;
            end
            if mclear==true
                        close all;
            elseif mclear==false
               if ~isempty(get(0,'CurrentFigure'))
                hold on;
               end
            else
                figure;
            end
            
            switch func
            case ''
            case 'plot3'
               [varargin indexes]= Utils.readParam(varargin, 'indexes', 1:P.n);
               if d3d
                  P=P.getRows(indexes);
                  L=Plot3d();
                  L.showDouble(P);
               else
                P=P.m(indexes,:,:); % convert to unweighted points
                if animate
                    Utils.animateSets(P);
                else
                    plot3(P(:,3), P(:,2), P(:,1), varargin{:});
                    grid on;
                end
               end
            otherwise 
                error('undefined');
            end % switch
            if ~isempty(title)
                set(gcf, 'Name', title);
            end
            if pau
                disp('press any key to continue...');
                pause;
            end
        end % show

        function breakpoints (save)
            mlock;
            persistent s;
            if strcmp(save,'save');       
                s=dbstatus('-completenames'); 
            else
                dbstop(s);  
                munlock;
                %# do some cleanup 
                clear s;
            end
        end

        % clear all while keeping breakpoints
        function clearAll()
            s=dbstatus('-completenames'); 
            save('myBreakpoints.mat','s'); 

            %# if you're clearing, you may as well just clear everything 
            %# note that if there is stuff stored inside figures (e.g. as callbacks), not all of  
            %# it may be removed, so you may have to 'close all' as well 
            clear all; %hidden
            clear classes;
            load('myBreakpoints.mat') ;
            dbstop(s);  

            %# do some cleanup 
            clear s;
            delete('myBreakpoints.mat');
        end
        function S=fieldName2Struct(fieldName)
                 fieldName=['.' fieldName];
                 [s e]=stok(fieldName, '.(){}');
                 nTokens=length(s);
                 fieldNames=cell(1,nTokens);
                 for i=1:nTokens
                     fieldNames{i}=fieldName(s(i):e(i));
                 end
                 c={};
                 %s(nTokens+1)=' ';
                 for i=1:nTokens
                      if strcmp(fieldName(s(i)-1),'.')
                        c{2*i-1}='.';
                        c{2*i}=fieldName(s(i):e(i));
                      elseif strcmp(fieldName(s(i)-1),'(')
                        c{2*i-1}='()';
                        c{2*i}={str2num(fieldName(s(i):e(i)))};
                      elseif strcmp(fieldName(s(i)-1),'{')
                        c{2*i-1}='{}';
                        c{2*i}={str2num(fieldName(s(i):e(i)))};
                      end % if
                 end
                 try
                 S = substruct(c{:});
                 catch
                     error('er');
                 end
        end % fieldName2struct        
        function fieldValue=getFieldValue(obj, fieldName)
                 try
                    S=Utils.fieldName2Struct(fieldName);
                     fieldValue=subsref(obj,S);
                 catch exception
                     fieldName
                     throw(exception)
                 end % try
        end % getFieldValue
        
        % Does this method exists for this object.
        % method is a string.
        function result=methodExists(object, method)
            mc = metaclass(object);
            if length(findobj(mc.MethodList,'Name',method))>=1
                result=true;
            else result=false;
            end
        end 
        % return string representation of an object,
        % or just the varialbe if it's not an object
        function s=asString(v)
                if isobject(v)
                    %if Utils.methodExists(v, 'asString')
                    %    s=v.asString;
                    %else
                        s=class(v);
                    %end
                else 
                    if iscell(v)
                        s=v{1};
                    else
                        s=v;
                    end
                end % outer if
        end
        function obj=setFieldValue(obj, fieldName, fieldValue)
                 S=Utils.fieldName2Struct(fieldName);
                 try
                    obj=subsasgn(obj, S, fieldValue);
                 catch exception
                     fieldName
                     fieldValue
                     throw(exception)
                 end % try
        end
        
        function C= setProd(cellArray)
            C = Matrix(setprod(cellArray));
        end
        function error=ratio(opt,approx, absolute)
            if nargin<3
                absolute=true;
            end
            if absolute
                error=abs(approx)/opt;
            else
                error=approx/opt;
            end
            error= Utils.numericCheck (opt,approx,error);
        end
        function error= numericCheck (opt,approx,error)
            if (abs(error)<100*eps)
               error=0;
           end
            if(and(opt < 100*eps,approx < 100*eps))
                error = 0;
            end
        end
        function error = relativeError (opt, approx, absolute)
            if nargin<3
                absolute=true;
            end
            % return a number between 0 to 1, represents the relative error of
            % approximating approx by opt
            if absolute
                error=abs(approx-opt)/opt;
            else
                error=(approx-opt)/opt;
            end
            error= Utils.numericCheck (opt,approx,error);
        end % function relativeError
        
        function varargout =complete(optargs, varin)
            % Some function that requires 2 inputs and has some optional inputs.
            % only want 3 nonempty optional inputs at most
            varin=varin{1};
            numvarargs = find(~cellfun('isempty',varin));
            if length(numvarargs) > length(optargs)
                error(['requires at most'  length(optargs) 'optional inputs']);
            end
            % now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(numvarargs) = varin(numvarargs);
            % Place optional args in memorable variable names
            varargout = optargs;
        end
        function ownerObj=completeProps(ownerObj, optargs, varin)
               s=[];
               for i=1:length(varin)
                   s=char(strcat(s, ' ownerObj.', varin(i)));
               end
               b=['[' s ']=Utils.complete(optargs,{{', s,'}});'];
               eval(b);
        end
        
        % animate set of points
        function animateSets(P)
            close all;
            % replace with an image of your choice
            h=figure;
            hold on;
            %axis ([min(x) max(x) min(y) max(y)])
            img = imread('C:\papers\old\ClusteringMotion\sf.bmp');
            imagesc(img);
            ax=axis;
             x = P(:,2);
             y = P(:,1);
             t=[min(x) max(x) min(y) max(y)];
             % normalize to 0..1
             x=(x-t(1))/(t(2)-t(1));
             y=(y-t(3))/(t(4)-t(3));

             % normalize to ax
             x=ax(1)+x*(ax(2)-ax(1));
             y=ax(3)+y*(ax(4)-ax(3));
             
            hold on;
            
%            colormap(gray);
            % set the range of the axes
            % The image will be stretched to this.
             pid=P(:,4);
             paus=0.001;
            % h=figure; % force view
             hold on;
             nPoints=2;
                 p=cell(1,nPoints);
                 grid off;
                 set(gca, 'YTick', []);
                 set(gca, 'XTick', []);
            
                 for ii=1:numel(x),
                       nPoint = pid(ii);
                       if isempty(p{nPoint})
                            p{nPoint}=plot(x(1),y(1),'ro') ;
                       else
                           set(p{nPoint},'xdata',x(ii),'ydata',y(ii)) ;
                       end 
                       %pause(paus) ;
                       drawnow ; % visibility
                 end
             end
          
        % animate set of points
        function animate(P)
            close all;
             x = 1:size(P,1);
             y = P(:,1);
             z = P(:,2);
             v{1}= [-21.5,2];
             v{2}=[90 0];
             paus=0.001;
              pos{1}=[ 9     5   409   327];
              pos{2}=[435   -14   429   346];
             h{1}=figure; % force view
             hold on;
             h{2}=figure;
             
             for i=1:2
                 figure(h{i});
                 set(h{i}, 'Position', pos{i});
                 p{i}=plot3(x(1),y(1),z(1),'b-',x(1),y(1),z(1),'ro') ;
                 axis ([min(x) max(x) min(y) max(y) min(z) max(z) ])
                 grid off;
                 set(gca, 'YTick', []);
                 set(gca, 'XTick', []);
                 set(gca, 'ZTick', []);
                 view(v{i});
             end
%             pause;
   
             %axis([-3 3 -4 14]) ; % specify a strange axis, that is not changed
            
             for ii=2:numel(x),
                 for i=1:2
                     figure(h{i});
                       % update plots
                       set(p{i}(1),'xdata',x(1:ii),'ydata',y(1:ii), 'zdata',z(1:ii)) ;
                       set(p{i}(2),'xdata',x(ii),'ydata',y(ii), 'zdata',z(ii)) ;
                       %pause(paus) ;
                       drawnow ; % visibility
                 end
             end
          end
      end % methods
    
end

