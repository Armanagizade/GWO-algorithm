lowerb = 20;
uperb = 100;
npop  = 10;
nvar = 3;
maxiter = 50;
pop = zeros(2,nvar,npop);
result = zeros(1,maxiter+1);
a = 2 - linspace(0,2,maxiter);
x = zeros(3,nvar);
xnew = zeros(1,nvar);
pos = 1;
fit = 2;
alpha = 1.5 * ones(1,3);
beta = 1* ones(1,3);
delta = 0.5* ones(1,3);
%random pop-------------------------------------------------------------------------------------------------------
pop(pos,:,:) = unifrnd(lowerb,uperb,1,nvar,npop);
for i=1:npop
    pop(fit,1,i) = pid(pop(pos,:,i)); 
end
[pppp ,in] = sort(pop(fit,1,:));
pop = pop(:,:,in);
result(1) = pop(fit,1,1);
disp(pop);
%main loop -------------------------------------------------------------------------------------------------------
for i=1:maxiter    
    alpha= pop(:,:,1);
    beta = pop(:,:,2);
    delta = pop(:,:,3);
    a1 = i^3;
    a2 = i^2;
    a3 = i;
    for k=1:npop        
        A = (2 * a(i) * rand(3,3)) - a(i);
        C = 2*rand(3,3);
        x(1,:) = alpha(pos,:) - (A(1,:) .* abs(C(1,:).*alpha(pos,:) - pop(pos,:,k)));
        x(2,:) = beta(pos,:) - (A(2,:) .* abs(C(2,:).*beta(pos,:) - pop(pos,:,k)));
        x(3,:) = delta(pos,:) - (A(3,:) .* abs(C(3,:).*delta(pos,:) - pop(pos,:,k)));
        xnew = (a1 * (x(1,:)) +  (a2 * x(2,:)) +  (a3 * x(3,:)))/ (a1+a2+a3) ;
        if (all(xnew > lowerb)) & (all(xnew < uperb)) 
            fnew = pid(xnew); 
            if pop(fit,1,k) > fnew 
                pop(fit,1,k) = fnew;
                pop(pos,:,k) = xnew;
            
            end
        end
    end
    [pppp ,in] = sort(pop(fit,1,:));
    pop = pop(:,:,in);
    result(i+1) = pop(fit,1,1);
    disp(i)
    disp(result(i));
end
disp(pop(:,:,1))
plot([1:maxiter+1],result,LineWidth=1.3)
xlabel('iteration')
ylabel('IAE')
%function --------------------------------------------------------------------------------------------------------
function[IAE] = pid(x)
    syms t  positive;
    syms s  ;
    kp = x(1);
    ki = x(2);
    kd = x(3);
    l = sym(1);
    u = sym(1);
    ys =@(s) 1/s -(((0.913242 / (1.39*s^2 + 1.215*s + 0.913242))*(kd*s+kp+ki/s)/((0.913242 / (1.39*s^2 + 1.215*s + 0.913242))*(kd*s+kp+ki/s)+1))*1/s);
    et =ilaplace(ys, s, t);
    f1 = vpa(et);
    f1 = abs(f1);
    IAE = vpaintegral(f1, t, 0, inf);
end