function [posArr,velArr] = nBodyMatrixGen(tMax,tStep,objArr)   
    t = 0;
    n = 1;
    posArr1 = [objArr.position];
    posArr = zeros((tMax/tStep + 1)*9,1);
    posArr(1:9) = posArr1;

    velArr1 = [objArr.velocity];
    velArr = zeros((tMax/tStep + 1)*9,1);
    velArr(1:9) = velArr1;

    %bodiesVec = objArr;
    while t < tMax
        for i = 1:max(size(objArr))
            objArr(i) = objArr(i).netAcceleration([objArr(1:i-1) objArr(i+1:end)]);
        end

        for i = 1:max(size(objArr))
            objArr(i) = objArr(i).integrate(tStep);
        end
        posArrNew = [objArr.position];
        posArr(n*9+1:(n+1)*9) = posArrNew;

        velArrNew = [objArr.velocity];
        velArr(n*9+1:(n+1)*9) = velArrNew;
        %bodiesVec(end+1,:) = objArr;
        t = t + tStep;
        n = n + 1;
    end
end