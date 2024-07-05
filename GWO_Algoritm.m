%ingredients
lowerbound = -1;
uperbound = 10;
n_population  = 50;
n_variable = 3;
maxiter = 100;
pop = zeros(2,n_variable,n_population);
result = zeros(1,maxiter);
a = 2 - linspace(0,2,maxiter);
x1 = zeros(1,n_variable); 
x2 = zeros(1,n_variable);
x3 = zeros(1,n_variable);
xnew = zeros(1,n_variable);
position = 1;
fitness = 2;
%random pop-------------------------------------------------------------------------------------------------------
for i=1:n_population
    pop(position,:,i) = rand(1,n_variable)*(uperbound-lowerbound)+lowerbound; 
    pop(fitness,1,i) = sphere(pop(position,:,i));
end
[pppp ,in] = sort(pop(fitness,1,:));
pop = pop(:,:,in);
%main loop -------------------------------------------------------------------------------------------------------

for i=1:maxiter    
    alpha= pop(:,:,1);
    beta = pop(:,:,2);
    delta = pop(:,:,3);
    for k=1:n_population        
        for j=1:n_variable
            A1 = (2 * a(i) * rand()) - a(i);
            A2 = (2 * a(i) * rand()) - a(i);
            A3 = (2 * a(i) * rand())- a(i);
            C1 = 2*rand();
            C2 = 2*rand();
            C3 = 2*rand();
            x1(j) = alpha(position,j) - (A1 * abs(C1*alpha(position,j) - pop(position,j,k)));
            x2(j) = beta(position,j) - (A2 * abs(C2*beta(position,j) - pop(position,j,k)));
            x3(j) = delta(position,j) - (A3 * abs(C3*delta(position,j) - pop(position,j,k)));
            xnew(j) = (x1(j) + x2(j) + x3(j))/3 ;
        end
        fnew = sphere(xnew);
        if pop(fitness,1,k) > fnew & (xnew > lowerbound) & all(xnew < uperbound)
            pop(fitness,1,k) = fnew;
            pop(position,:,k) = xnew;
        end
    end
    [pppp ,in] = sort(pop(fitness,1,:));
    pop = pop(:,:,in);
    result(i) = pop(fitness,1,1);
    disp(result(i));
end
plot([1:maxiter],result)
%function --------------------------------------------------------------------------------------------------------
function[y] = sphere(x)
    y = sum(x.^2);
end