function [theta, rho]=sendToPolar(x,y)
    rho=[];
    theta=[];
    theta=atan2(x,y);
    theta=mod(theta,2*pi);
    % theta(1)=0;
    % theta(length(theta)+1) = theta(1)+2*pi;
    rho=sqrt(x.^2+y.^2);
    % rho(length(rho)+1)=rho(1);
    
    %Sort them
    thetaRho = [theta' rho'];
    sorted_thetaRho = sortrows(thetaRho,1);
    theta = sorted_thetaRho(:,1);
    rho = sorted_thetaRho(:,2);
    theta = theta';
    rho = rho';
    
end
