classdef FEPackObject < matlab.mixin.Copyable
    % FEPackObject Attributes and methods shared by all FEPack objects.

    properties

        % name for error/log/debug-messages
        name = 'FEPack Object';

    end

    methods

        % Sets the name to a random string of length 10
        function randomName(obj)

            symbols = ['a':'z' 'A':'Z' '0':'9'];
            nums = randi(numel(symbols),[1 10]);
            obj.name = symbols (nums);

        end

    end

end
