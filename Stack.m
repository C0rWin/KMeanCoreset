classdef Stack < handle
    % implements stack data structure. The itmes are stored in a cell
    % arrray, and can be of any type.
    properties
        stack={};
    end
    
    methods
        function obj=Stack(inputStack)
            if nargin>0
                obj.stack = inputStack.stack;
            end
        end
        function clear(obj)
            obj.stack={};
        end
       function result=isEmpty(obj)
           if isempty(obj.stack)
                result=true;
           else result=false;
           end
        end % function isEmpty
        function item=top(obj)
          if obj.isEmpty() 
              error('top from empty stack');
          end
          item= obj.stack{end};
        end
        function [item]=pop(obj)
          if obj.isEmpty() 
              error('pop from empty stack');
          end
          item= obj.stack{end};
          obj.stack(end)=[];
        end % function pop       
        function push(obj, item)
            obj.stack{end+1}=item;
        end
    end % methods
    
end

