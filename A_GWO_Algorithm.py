#p1 : ingredients
import numpy as np
import random
import matplotlib.pyplot as plt
lowerbound , uperbound = 0 , 10
n_population , n_variables , maxiter = 50 , 3 , 50
#population
pop = np.zeros((n_population,2,n_variables))
result = np.zeros(maxiter)
a = 2-np.linspace(0,2,maxiter)
x1,x2,x3 = np.zeros((n_variables)) , np.zeros((n_variables)) , np.zeros((n_variables))
xnew = np.zeros((n_variables))
#indexes
positon , fitness = 0 , 1
#p2: function---------------------------------------------------------------------------------------------------------
def sphere(x):
    return np.sum(x**2)
#p3 : random population-----------------------------------------------------------------------------------------------
for i in range(n_population):
    pop[i][positon] = np.random.uniform(lowerbound,uperbound,n_variables) 
    pop[i][fitness][0] = sphere(pop[i][positon])
pop = pop[np.argsort(pop[: , fitness , 0])]
#main loop------------------------------------------------------------------------------------------------------------
for i in range(maxiter):    
    alpha= pop[0]
    beta = pop[1]
    delta = pop[2]
    for k in range(n_population):        
        A1, A2, A3 = (2 * a[i] * random.random()) - a[i] , (2 * a[i] * random.random()) - a[i], (2 * a[i] * random.random())- a[i]
        C1, C2, C3 = 2 * random.random(), 2*random.random(), 2*random.random()
        x1 = alpha[positon] - (A1 * abs(C1*alpha[positon] - pop[k][positon]))
        x2 = beta[positon] - (A2 * abs(C2*beta[positon] - pop[k][positon]))
        x3 = delta[positon] - (A3 * abs(C3*delta[positon] - pop[k][positon]))
        xnew = (x1 + x2 + x3)/3
        if random.random() > (float(i/maxiter)) :
            index = random.randint(0,n_variables-1)
            A1, A2, A3 = (2 * a[i] * random.random()) - a[i] , (2 * a[i] * random.random()) - a[i], (2 * a[i] * random.random())- a[i]
            C1, C2, C3 = 2 * random.random(), 2*random.random(), 2*random.random()
            x1[index] = alpha[positon][index] - (A1 * abs(C1*alpha[positon][index] - pop[k][positon][index]))
            x2[index] = beta[positon][index] - (A2 * abs(C2*beta[positon][index] - pop[k][positon][index]))
            x3[index] = delta[positon][index] - (A3 * abs(C3*delta[positon][index] - pop[k][positon][index]))
            xnew[index] = (x1[index] + x2[index] + x3[index])/3               
        fnew = sphere(xnew)    
        if pop[k][fitness][0] > fnew and all(xnew > lowerbound) and all(xnew < uperbound):
            pop[k][fitness][0] = fnew
            pop[k][positon] = xnew            
    pop = pop[np.argsort(pop[: , fitness , 0])] 
    result[i] = pop[0][fitness][0]
    if i % 9 == 0 :print(f"iter = {i} --- best = {result[i]} --- alpha = {alpha[positon]}")    
fig = plt.figure(figsize=(10,5))
axs = fig.add_axes([0,0,0.9,0.9])
axs.plot(np.arange(maxiter),result )