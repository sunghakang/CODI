function [diffusedIm,objval,Itercount,CPUcount] = AlgIEdge(SeedIm,EdgeIm,X_D,mu,theta,eta, perc)
% This algorithm is the Edge-weighted harmonic variational optimization model
% proposed in the paper
% AlgIEdge(SeedIm,EdgeIm,X_D,mu,theta, eta, perc)
%     -inputs:
%      SeedIm       : Seed Image
%      EdgeIm       : Mask Image of the image that $g$ applied onto
%      mu           : multiplier on the constraint V−U=0.
%      lambda       : λ, the Lagrange multiplier associated with the linear constraint V−U=0.
%      eta          : η, fidelity parameter

%     -output :
%      diffusedIm  : Diffused image
%      objval       : Energy Function
%      Itercount    : Iteration Number
%      CPUcount     : CPU time

objval = [];
Itercount =[];
CPUcount = [];
X_D = X_D.*EdgeIm;
u0 = SeedIm;
u=u0;
v =u0;
lambda = u0*0;
G0 = max(max(EdgeIm));
C = EdgeIm-G0;
[n,m]=size(u0);

A = zeros(n, m);
for i = 1:n
    for j = 1:m
        A(i,j)= theta+mu- 4*G0 *(cos(2*pi*(i-1)/n)+ cos(2*pi*(j-1)/m)-2);
    end
end

tic

obj= sum(sum( EdgeIm.*(dxf(u).^2+dyf(u).^2)+0.5*eta*X_D.*(u-SeedIm)));
objval = [objval,obj];
Itercount =[Itercount,0];
CPUcount = [CPUcount,toc];
t = 1;
err = 100;
while err>perc

    ux = dxf(u);
    uy = dyf(u);
    g = theta*u + 2* dxb(C.* ux)+2*dyb(C.* uy)+ mu *v +lambda;
    u = real(ifftn(fftn (g)./ A));
    v =(eta*X_D.*u0+mu*u-lambda)./(eta*X_D+mu);
    lambda = lambda + mu*(v-u);

    obj= sum(sum( EdgeIm.*(dxf(u).^2+dyf(u).^2)+0.5*eta*X_D.*(u-SeedIm)));
    objval = [objval,obj];
    Itercount =[Itercount,t];
    CPUcount = [CPUcount,toc];
    err = (objval(end-1)-objval(end))/objval(end-1);
    t = t+1;
end

diffusedIm=u;

end

function v=dyb(u)
v = u(:,:) - u(:,[end 1:end-1]);
end


function v=dxb(u)
v = u(:,:) - u([end 1:end-1],:);
end


function v=dxf(u)
v = u([2:end 1],:) - u(:,:);
end


function v=dyf(u)
v = u(:,[2:end 1]) - u(:,:);
end
