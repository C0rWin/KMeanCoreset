classdef HandleObject < handle
    properties
        
    end
    
    methods
        function saveObj(obj, filename) %#ok<INUSD
            loadedObj=obj;
            save(filename, 'loadedObj');
        end
    end
    methods (Static)
        function loadedObj=loadObj (filename) %#ok<MANU,INUSD>
            load(filename, 'loadedObj'); %#ok<NASGU>
        end
    end
    
end

