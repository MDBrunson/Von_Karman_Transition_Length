import numpy as np
import time

diam_s = 4
diam_l = 6
accuracy = 1e-6
def transition_length(diam_s, diam_l):
    start = time.process_time()
    x = 0
    h = accuracy
    nextx = 0
    dx = 1

    def func(input, D, T):
        y = ((D/(2*np.sqrt(3.1415926358979)))*np.sqrt(np.arccos(1-((2*input)/(3*D))) - ((np.sin(2*(np.arccos(1-((2*input)/(3*D))))))/2))-(T/2))
        return y
    
    def deriv(x_pos, D_l, D_s):
        return (func(x_pos+h, D_l, D_s)-func(x_pos, D_l, D_s))/h
    
    count = 0
    while dx >= accuracy:
        nextx = (-func(x, diam_l, diam_s)/deriv(x,diam_l, diam_s))+nextx
        dx = nextx - x
        x = nextx 
        count +=1
        
    trans_length = 3*diam_l - x
    end = time.process_time()
    return trans_length, x, count, (end - start)


print('Transition Length, dist from tip/x-int of small radius intercept, iterations, time:')
print(transition_length(diam_s, diam_l))