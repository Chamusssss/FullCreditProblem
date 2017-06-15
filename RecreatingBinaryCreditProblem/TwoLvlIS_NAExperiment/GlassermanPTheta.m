function [pTheta,thetaVec] = GlassermanPTheta(pncz,LGC,EAD,tail)

    [~,~,NMC] = size(pncz);
    pTheta = pncz;
    thetaVec = zeros(1,NMC);
    
    weights = LGC.*EAD;
    
    qnc = @(theta,pnc) pnc.*exp(weights.*theta)./sum(pnc.*exp(weights.*theta),2);
    mu = @(theta,pnc) computeMu(weights,qnc(theta,pnc));
    sigma = @(theta,pnc) computeSigma(LGC,EAD,qnc(theta,pnc));
    psi = @(theta,pnc) sum(log(sum(pnc.*exp(weights.*theta),2)),1);

    for i=1:NMC
         pnc = pncz(:,:,i);
         energy = @(theta) ((1-normcdf((tail - mu(theta,pnc)) / sigma(theta,pnc)))*...
             exp(2*(sigma(theta,pnc)^2)*theta^2 - 2*mu(theta,pnc)*theta) - ...
             ((1-normcdf((tail - mu(theta,pnc)) / sigma(theta,pnc)))^2)*...
              exp((sigma(theta,pnc)^2)*theta^2 - 2*mu(theta,pnc)*theta))*...
              exp(2*psi(theta,pnc))...
         ;
         options = optimset('LargeScale','off', 'display', 'off');
         intialGuess = 0;
         if(i > 1); intialGuess = thetaVec(i-1); end
         [theta, fval, exitflag, output] = fminunc(energy, intialGuess, options);
         pTheta(:,:,i) = qnc(theta,pnc);
         thetaVec(i) = theta;
     end
    
    function [mu] = computeMu(weights,p)
        mu = sum(sum(weights.*p));
    end

    function [sigma] = computeSigma(LGC,EAD,p)
        [N,C] = size(LGC);
        A = zeros(N,(C-1)*C/2);
        index = 1;
        for a=1:C
            for b=1:(a-1)
                A(:,index) = ((LGC(:,a) - LGC(:,b)).^2).*p(:,a).*p(:,b);
                index = index + 1;
            end
        end
        B = sum(A,2);
        sigma = sqrt(((EAD.^2)')*B);
    end
end