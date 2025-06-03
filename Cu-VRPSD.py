### This code provides all numerical calculations and generates Figures 2-6 in the paper.
### The code is organized into four parts

### PART I ############################################################################################################################
import matplotlib.pyplot as plt
import numpy as np
from gurobipy import *
plt.rcParams['font.family'] = 'DejaVu Sans'
############################################################################################################################
def R11(gamma):   ### APPROX.1(\lambda, 0.5, p)
    alpha = 1.5
    
    labda = min(1, 4*gamma/alpha)

    theta = 0.5

    p = (1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))/(1/(2*labda)+1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))
    
    ### R(1)
    L = gamma*(alpha+p*2/labda+(1-p)*1/(theta*labda)*(2*labda-theta*labda)/(labda-theta*labda))+(alpha*(p*labda+(1-p)*theta*labda)/2+p/2+(1-p)/2*(2*labda-theta*labda)/(labda-theta*labda))
    L = L/(gamma+0.5)
    
    ### R(\infty)
    R = gamma*alpha+alpha*(p*labda+(1-p)*theta*labda)/2
    R = R/gamma

    return max(L,R)
############################################################################################################################
def R12(gamma):   ### APPROX.1(\lambda, 0.6677, p)
    alpha = 1.5
    
    labda = min(1, 4*gamma/alpha)

    theta = 0.6677

    p = (1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))/(1/(2*labda)+1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))
    
    ### R(1)
    L = gamma*(alpha+p*2/labda+(1-p)*1/(theta*labda)*(2*labda-theta*labda)/(labda-theta*labda))+(alpha*(p*labda+(1-p)*theta*labda)/2+p/2+(1-p)/2*(2*labda-theta*labda)/(labda-theta*labda))
    L = L/(gamma+0.5)
    
    ### R(\infty)
    R = gamma*alpha+alpha*(p*labda+(1-p)*theta*labda)/2
    R = R/gamma

    return max(L,R)
############################################################################################################################
def R13(gamma):   ### APPROX.1(\lambda, \theta, p)
    alpha = 1.5
    
    labda = min(1, 4*gamma/alpha)
    
    temvalue = 111
    temtheta = 111
    
    for i in range(9999):
        
        ### \theta
        theta = (i+1)/10000
        
        p = (1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))/(1/(2*labda)+1/(2*(labda-theta*labda))+gamma/(theta*labda*(labda-theta*labda)))
       
        ### R(1)
        L = gamma*(alpha+p*2/labda+(1-p)*1/(theta*labda)*(2*labda-theta*labda)/(labda-theta*labda))+(alpha*(p*labda+(1-p)*theta*labda)/2+p/2+(1-p)/2*(2*labda-theta*labda)/(labda-theta*labda))
        L = L/(gamma+0.5)
       
        ### R(\infty)
        R = gamma*alpha+alpha*(p*labda+(1-p)*theta*labda)/2
        R = R/gamma
        
        if max(L,R) < temvalue:
            temvalue = max(L,R)
            temtheta = theta
    
    return temvalue
############################################################################################################################
def R21(gamma):   ### APPROX.2
    alpha = 1.5
    
    ### R(1)
    L = gamma*(alpha+1.5)+2/3*alpha+1
    L = L/(gamma+0.5)
    L = L + (6*gamma-1)/(24*gamma+4)
    
    ### R(\infty)
    R = alpha + 2/3*alpha/gamma + (6*gamma-1)/(24*gamma+4)
    
    return max(L,R)
############################################################################################################################
def R22(gamma):   ### APPROX.2(\sigma=1) the maximum value obtained from the linear programs
    return max(lp1(gamma, 1), lp2(gamma, 1))
############################################################################################################################
### CASE 1
def lp1(gamma, sigma):   ### the value of the first linear program
    ############################################################################################################################
    m = Model('calculate ratio')
    y = m.addVar(vtype=GRB.CONTINUOUS, name="y"); m.setParam('OutputFlag', 0)
    ############################################################################################################################
    LB = gamma*max(sigma, 1) + 0.5
    
    alpha = 1.5; 
    N=300
    dx = []
    for i in range(N):
        dx.append(i+1)
    ############################################################################################################################    
    r0 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    r1 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    r2 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    A = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    B = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    m.addConstr( y <= 0.5*A+ 0.5*B )
    ############################################################################################################################
    m.addConstr( A <= gamma*(alpha*sigma + 1.5*r1[1/3*N] 
                             + 3*(r1[2/3*N]-r1[1/3*N])-0.5*(r0[2/3*N]-r0[1/3*N]) 
                             + 1.5*(r1[N]-r1[2/3*N])+0.5*(r0[N]-r0[2/3*N])) / LB 
                
                + (2/3*alpha*sigma + r1[1/3*N]-0.75*r2[1/3*N] 
                   + 0.75*(r2[2/3*N]-r2[1/3*N])+0.5*(r1[2/3*N]-r1[1/3*N]) 
                   + 1.5*(r2[N]-r2[2/3*N])-(r1[N]-r1[2/3*N])+2/3*(r0[N]-r0[2/3*N])) / LB )
    ############################################################################################################################
    m.addConstr( B <= gamma*(alpha*sigma + 1.5*r1[1/3*N]) / LB 
                
                + (2/3*alpha*sigma + r1[1/3*N]-0.75*r2[1/3*N]) / LB 
                
                + (gamma*(r0[N]-r0[1/3*N]) + 0.5*(r1[N]-r1[1/3*N])) / LB )
    ############################################################################################################################
    m.addConstr( (gamma*(r0[N]-r0[1/3*N]) + 0.5*(r1[N]-r1[1/3*N])) <= LB )
    ############################################################################################################################
    for i in dx:
        if i > 1:
            m.addConstr( r2[i]-r2[i-1] >= (i-1)/N * (r1[i]-r1[i-1]) ) 
            m.addConstr( r2[i]-r2[i-1] <= i/N * (r1[i]-r1[i-1]) )
            m.addConstr( r1[i]-r1[i-1] >= (i-1)/N * (r0[i]-r0[i-1]) )
            m.addConstr( r1[i]-r1[i-1] <= i/N * (r0[i]-r0[i-1]) )
        else:
            m.addConstr( r2[i] >= (i-1)/N * r1[i] ) 
            m.addConstr( r2[i] <= i/N * r1[i] )
            m.addConstr( r1[i] >= (i-1)/N * r0[i] )
            m.addConstr( r1[i] <= i/N * r0[i] )
    ############################################################################################################################        
    m.addConstr( r1[N] >= 1 ); m.addConstr( r1[N] <= 1 );
    ############################################################################################################################
    
    m.setObjective(y, GRB.MAXIMIZE); m.optimize()
    if m.status == GRB.OPTIMAL: return m.objVal
    else: return -1
############################################################################################################################
### CASE 2
def lp2(gamma, sigma):   ### the value of the second linear program
    ############################################################################################################################
    m = Model('calculate ratio')
    y = m.addVar(vtype=GRB.CONTINUOUS, name="y"); m.setParam('OutputFlag', 0)
    ############################################################################################################################
    LB = gamma*max(sigma, 1) + 0.5
    
    alpha = 1.5; 
    N=300
    dx = []
    for i in range(N):
        dx.append(i+1)
    ############################################################################################################################    
    r0 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    r1 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    r2 = m.addVars(dx, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    A = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    B = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    m.addConstr( y <= 0.5*A+ 0.5*B )
    ############################################################################################################################
    m.addConstr( A <= gamma*(alpha*sigma + 1.5*r1[1/3*N] 
                             + 3*(r1[2/3*N]-r1[1/3*N])-0.5*(r0[2/3*N]-r0[1/3*N]) 
                             + 1.5*(r1[N]-r1[2/3*N])+0.5*(r0[N]-r0[2/3*N])) / LB 
                
                + (2/3*alpha*sigma + r1[1/3*N]-0.75*r2[1/3*N] 
                   + 0.75*(r2[2/3*N]-r2[1/3*N])+0.5*(r1[2/3*N]-r1[1/3*N]) 
                   + 1.5*(r2[N]-r2[2/3*N])-(r1[N]-r1[2/3*N])+2/3*(r0[N]-r0[2/3*N])) / LB )
    ############################################################################################################################
    m.addConstr( B <= gamma*(alpha*sigma + 1.5*r1[1/3*N]) / LB 
                
                + (2/3*alpha*sigma + r1[1/3*N]-0.75*r2[1/3*N]) / LB 
                
                + 1 )
    ############################################################################################################################
    m.addConstr( (gamma*(r0[N]-r0[1/3*N]) + 0.5*(r1[N]-r1[1/3*N])) >= LB )
    ############################################################################################################################
    for i in dx:
        if i > 1:
            m.addConstr( r2[i]-r2[i-1] >= (i-1)/N * (r1[i]-r1[i-1]) ) 
            m.addConstr( r2[i]-r2[i-1] <= i/N * (r1[i]-r1[i-1]) )
            m.addConstr( r1[i]-r1[i-1] >= (i-1)/N * (r0[i]-r0[i-1]) )
            m.addConstr( r1[i]-r1[i-1] <= i/N * (r0[i]-r0[i-1]) )
        else:
            m.addConstr( r2[i] >= (i-1)/N * r1[i] ) 
            m.addConstr( r2[i] <= i/N * r1[i] )
            m.addConstr( r1[i] >= (i-1)/N * r0[i] )
            m.addConstr( r1[i] <= i/N * r0[i] )
    ############################################################################################################################        
    m.addConstr( r1[N] >= 1 ); m.addConstr( r1[N] <= 1 );
    ############################################################################################################################
    
    m.setObjective(y, GRB.MAXIMIZE); m.optimize()
    if m.status == GRB.OPTIMAL: return m.objVal
    else: return -1
    
### PART II ################################################################################################################################  
if 0:
    ### The following code can be used to obtained Figures 2, 3, and 6.
    x = np.linspace(0.1, 0.374, 50) ### \gamma\in (0,0.375)
    xx = np.linspace(0.375, 2, 150) ### \gamma\in (0.375, 2)
    
    plt.figure(figsize=(8, 6), dpi=600)
    ##########################
    y1 = np.array([R11(xi) for xi in x]) ### APPROX.1(\lambda, 0.5, p)
    y2 = np.array([R12(xi) for xi in xx]) ### APPROX.1(\lambda, 0.6677, p)
    y3 = np.array([R13(xi) for xi in xx]) ### APPROX.1(\lambda, \theta, p)
    y4 = np.array([R21(xi) for xi in xx]) ### APPROX.2
    y5 = np.array([R22(xi) for xi in xx]) ### APPROX.2(\sigma=1) the maximum value obtained from the linear programs
    
    plt.plot(x, y1, label='$APPROX.1(\lambda,0.5,p)$', color = 'black')
    plt.plot(xx, y2, label='$APPROX.1(1,0.6677,p)$', color = 'blue')
    plt.plot(xx, y3, label=r'$APPROX.1(1,\theta,p)$', color = 'green')
    plt.plot(xx, y4, label='$APPROX.2$', color = 'red')
    plt.plot(xx, y5, label='$APPROX.2(\sigma=1)$', color = 'orange')
    
    point_x = 0.375
    point_y_f = R11(point_x)
    plt.scatter([point_x], [point_y_f], color='black')
    plt.annotate(f'({point_x}, {point_y_f:.4f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,10), ha='center')
    
    point_x = 0.375
    point_y_f = R12(point_x)
    plt.scatter([point_x], [point_y_f], color='blue')
    plt.annotate(f'({point_x}, {point_y_f:.5f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,10), ha='center', color='blue')
    
    point_x = 1.444
    point_y_f = R12(point_x)
    plt.scatter([point_x], [point_y_f], color='red')
    plt.annotate(f'({point_x}, {point_y_f:.5f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,10), ha='left', color='red')
    
    point_x = 1.033
    point_y_f = R22(point_x)
    plt.scatter([point_x], [point_y_f], color='orange')
    plt.annotate(f'({point_x}, {point_y_f:.5f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,10), ha='left', color='orange')
    ##########################
    plt.legend(); plt.grid(True); 
    #plt.savefig('Cu-VRPSD.pdf', format='pdf')
    plt.show()
    
### PART III ###############################################################################################################################
if 0:
    ### The following code can be used to obtained Figure 4.
    gamma = 1.444
    x = np.linspace(1, 3, 300)
    
    plt.figure(figsize=(8, 6), dpi=600)

    y1 = np.array([lp1(gamma, xi) for xi in x])
    y2 = np.array([lp2(gamma, xi) for xi in x])


    plt.plot(x, y1, label='Case 1', color = 'red')
    plt.plot(x, y2, label=r'Case 2', color = 'blue')

    point_x = 1
    point_y_f = lp1(gamma, 1)
    plt.scatter([point_x], [point_y_f], color='blue')
    plt.annotate(f'({point_x}, {point_y_f:.5f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,5), ha='center', color='red')

    point_x = 1
    point_y_f = lp2(gamma, 1)
    plt.scatter([point_x], [point_y_f], color='blue')
    plt.annotate(f'({point_x}, {point_y_f:.5f})', (point_x, point_y_f), textcoords="offset points", xytext=(0,-13), ha='center', color='blue')
    ##########################
    plt.legend(); plt.grid(True); 
    #plt.savefig('Cu_VRPSD.pdf', format='pdf')
    plt.show()
    
### PART IV ################################################################################################################################
if 1:
    ### The following code can be used to obtained Figure 5.
    G = np.linspace(0.3, 2, 18)
    Gamma = np.array([ gi for gi in G])
    plt.figure(figsize=(8, 6), dpi=600)
    
    count=0
    
    for gi in Gamma:
        gamma = gi

        x = np.linspace(1, 3, 30)
        
        y1 = np.array([ max(lp1(gamma, xi), lp2(gamma, xi)) for xi in x])
        
        plt.plot(x, y1, color='black')
        
        R=count*0.1+0.3

        if count % 2 ==0: 
            point_x = 3
            point_y_f = max(lp1(gamma, 3), lp2(gamma, 3))
            plt.annotate(fr'{R:.1f}', (point_x, point_y_f), textcoords="offset points", xytext=(12,0), ha='center', color='black')
         
        count = count + 1
    ##########################
    plt.legend(); plt.grid(True); 
    #plt.savefig('Cu_VRPSD.pdf', format='pdf')
    plt.show()

































