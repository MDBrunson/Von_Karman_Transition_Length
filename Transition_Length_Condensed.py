import numpy as np

#values predefined in FS_Pro or whatever
diam_s = 4
diam_l = 6


"""
--------------------------------------------- Find Transition Length --------------------------------------------------------------------
"""
# finds transition length from diameter of small with diameter of large. Assumes von karman nose cone length is 6* larger diameter
# desmos for future reference https://www.desmos.com/calculator/nhrjfdwia1


def transition_length(diam_s, diam_l):
    x = 0 
    accuracy = 1e-6 #corresponds to max dx for root algorithm and infinitessimal value for numerical derivative
    z = 0 # next step for x value on the Von Karman function,
    dx = 1 # change in x, updated later. 1 right now to allow while loop to begin. 

# Von Karman equation from x input, diameter, translation downward by smaller RADIUS
    def func(input, D, T):
        # Try to follow this:
        y = ((D/(2*np.sqrt(3.1415926358979)))*np.sqrt(np.arccos(1-((2*input)/(3*D))) - ((np.sin(2*(np.arccos(1-((2*input)/(3*D))))))/2))-(T/2))
        return y
    
# find numerical derivative with 
    def deriv(x_pos, D_l, D_s):
        return (func(x_pos+accuracy, D_l, D_s)-func(x_pos, D_l, D_s))/accuracy
    
#Newton-Rapshon root finding algorithm, finds 0 since the original function is translated down 
    while dx >= accuracy:
        z = (-func(x, diam_l, diam_s)/deriv(x,diam_l, diam_s))+z
        dx = z - x
        x = z 
        
    trans_length = 3*diam_l - x
    return trans_length

