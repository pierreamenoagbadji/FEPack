classdef FEPackObject < matlab.mixin.Copyable
    %>@brief FEPackObject Attributes and methods shared by all FEPack objects.

    properties

        % name for error/log/debug-messages
        name = 'FEPack_Object';

    end

    properties (Constant)

      % path to folder in cpp-friendly format (space handling, ...)
      pathCpp = 'home/pierre/Documents/Recherche/Code/FEPack';
      
      % path to folder in bash-friendly format (space handling, ...)
      pathBash = 'home/pierre/Documents/Recherche/Code/FEPack';

    end

    methods

        function FEObj = FEPackObject

          FEObj.randomName

        end

        % Sets the name to a random string of length 10
        function randomName(FEobj)

          symbols = ['a':'z' 'A':'Z' '0':'9'];
          nums = randi(numel(symbols),[1 10]);
          FEobj.name = symbols (nums);

        end
        
    end

end
