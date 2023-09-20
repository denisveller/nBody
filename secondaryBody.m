classdef secondaryBody < massiveBody
    
    methods  
        function f = forceExertedBy(obj, body2)
            m = body2(4);
            r = body2(1:3) - obj.position;
            f = ((obj.G*obj.mass*m).*r)/(norm(r)^3);
        end
        function obj = netAcceleration(obj,bodyArr,varargin)
            if size(varargin) == 0
                varargin = 0;
            end
            fNet = varargin{1};
            for i = 1:size(bodyArr,1)
                fNet = fNet + obj.forceExertedBy(bodyArr(i,:));
            end
            obj.acceleration = fNet/obj.mass;
        end
    end
end

