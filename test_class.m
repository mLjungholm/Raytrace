classdef test_class < handle
    properties
        a = [6,1];
        b = 1;
        c;
        d; e;
    end
    
    methods
        function this = test_class(varargin)
            this.a = varargin{1};
            this.c = cell2mat(reshape(varargin,[],1));
            this.d = nargin;
            this.e = size(varargin{1},1);
        end
        function print_a(this)
            disp(this.a)
        end
        function print_b(this)
            disp(this.b)
        end
        function mod_b(this)
            this.b = this.b + this.a(2);
        end
    end
end