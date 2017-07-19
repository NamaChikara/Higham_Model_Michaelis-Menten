function [out,num]=newton_higham(k,u,tol,iter,K1,K2,K3)

% k := time mesh size
% m := # of spatial grid points
% u := concentration after advective step
% tol := accuracy requirement for exiting Newton loop
% iter := max number of Newton iterations
% K1, K2, K3 := reaction rate constants

% Reaction equations

f1=@(x) -K1.*x(1)*(x(2)')+K2.*x(3);
f2=@(x) -K1.*x(1)*(x(2)')+(K2+K3).*x(3);
f3=@(x) K1.*x(1)*(x(2)')-(K2+K3).*x(3);
f4=@(x) K3.*x(3);

% Generate initial guess

uG=u+0.01*u;

% Initalize residual

res=zeros(size(u));

% Initialize outputs

out=NaN;
num=NaN;

% Solve for updated u value

for i=1:iter
    
    res(1)=uG(1)-u(1)-k.*f1(uG);
    res(2)=uG(2)-u(2)-k.*f2(uG);
    res(3)=uG(3)-u(3)-k.*f3(uG);
    res(4)=uG(4)-u(4)-k.*f4(uG);
    
    % Tolerance applied to Absolute Error
    %%{
    if norm(res,inf)<tol
        out=uG;
        num=i;
        break
    end
    %%}
    %{
    % Tolerance applied to Relative Error
    
    if i==1 
        res0=norm(res,inf);     % Set scaling factor for relative tolerance
    end
    
    if i>1 && norm(res,inf)/res0<tol
        out=u;
    
        num=i;
        if num<4
            k=k*2;
        else if num>8
            k=k/2;
            end
        end
    
        if time+k>tend
            time=tend;
        else time=time+k;
        end
        return
    end
    %}

    Jac11=ones(size(uG(1,:)))+(k*K1).*uG(2,:);       % dr1/duG1
    Jac21=(k*K1).*uG(2,:);                          % dr2/duG1
    Jac31=-(k*K1).*uG(2,:);                         % dr3/duG1
    Jac41=zeros(size(uG(1,:)));                      % dr4/duG1
    
    Jac12=(k*K1).*uG(1,:);                          % dr1/duG2
    Jac22=ones(size(uG(2,:)))+(k*K1).*uG(1,:);      % dr2/duG2
    Jac32=-(k*K1).*uG(1,:);                         % dr3/duG2
    Jac42=zeros(size(uG(2,:)));                      % dr4/duG2
    
    Jac13=-(k*K2).*ones(size(uG(3,:)));             % dr1/duG3
    Jac23=-(k*(K2+K3)).*ones(size(uG(3,:)));        % dr2/duG3
    Jac33=(1+k*(K2+K3)).*ones(size(uG(3,:)));       % dr3/duG3
    Jac43=-(k*K3).*ones(size(uG(3,:)));              % dr4/duG3
    
    Jac14=zeros(size(uG(4,:)));                      % dr1/duG4
    Jac24=zeros(size(uG(4,:)));                      % dr2/duG4
    Jac34=zeros(size(uG(4,:)));                      % dr3/duG4
    Jac44=ones(size(uG(4,:)));                       % dr4/duG4
    
    % Update guess
    
        Jac=[Jac11,Jac12,Jac13,Jac14;...
            Jac21,Jac22,Jac23,Jac24;...
            Jac31,Jac32,Jac33,Jac34;...
            Jac41,Jac42,Jac43,Jac44];
        step=-Jac\res;
        uG=uG+step
    
end

end

%}