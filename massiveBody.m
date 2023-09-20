classdef massiveBody
    %MASSIVEBODY Summary of this class goes here
    %Detailed explanation goes here
    properties
        mass;
        name;
        position;
        velocity;
        acceleration;
    end
    properties (Constant)
        G = 6.6743 * 10 ^ -11;
    end
    
    methods
        function obj = massiveBody(name, m,r,v)
            obj.name = name;
            obj.mass = m;
            obj.position = r*1000;
            obj.velocity = v*1000;
            obj.acceleration = [0,0,0];
        end
        function f = forceExertedBy(obj, body2)
            r = body2.position - obj.position;
            f = ((obj.G*obj.mass*body2.mass).*r)/(norm(r)^3);
        end

        function obj = netAcceleration(obj,bodyArr,varargin)
            if size(varargin) == 0
                varargin = 0;
            end
            fNet = varargin;
            for i = 1:length(bodyArr)
                fNet = fNet + obj.forceExertedBy(bodyArr(i));
            end
            obj.acceleration = fNet/obj.mass;
        end

        function obj = integrate(obj,tStep)

            obj.velocity = obj.velocity + obj.acceleration*tStep;
            obj.position = obj.position + obj.velocity*tStep;
        end
    end
end

