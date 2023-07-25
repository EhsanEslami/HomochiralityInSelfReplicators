
import numpy as np
import gillespy2
import matplotlib.pyplot as plt


k_a_value = 1
k_n_value = 1

samples = 1000                                 # number of the samples

class chirality(gillespy2.Model):
    def __init__(self, parameter_values=None):
        
        # First call the gillespy2.Model initializer.
        
        super().__init__(self)

        # Define parameters for the rates 
        k_a = gillespy2.Parameter(name='k_a', expression= k_a_value)
        k_n = gillespy2.Parameter(name='k_n', expression= k_n_value)
        k_d = gillespy2.Parameter(name='k_d', expression= 0.01)
        
        self.add_parameter([k_a, k_n, k_d])

        # Define variables for the molecular species representing A & D & L.
        
        a = gillespy2.Species(name='A', initial_value = 1000)
        d = gillespy2.Species(name='D', initial_value = 0)
        l = gillespy2.Species(name='L',  initial_value = 0)
        
        self.add_species([a , d , l])

        # The list of reactants and products for a Reaction object are
        # each a Python dictionary in which the dictionary keys are
        # Species objects and the values are stoichiometries of the
        # species in the reaction.
        
        r_a1 = gillespy2.Reaction(name="autocatalysis1", rate=k_a,
                                 reactants={a:1 ,d:1}, products={d:2})
        r_a2 = gillespy2.Reaction(name="autocatalysis2", rate=k_a,
                                 reactants={a:1 ,l:1}, products={l:2})
        r_n1 = gillespy2.Reaction(name="non_autocatalysis1", rate=k_n,
                                 reactants={a:1}, products={d:1})
        r_n2 = gillespy2.Reaction(name="non_autocatalysis2", rate=k_n,
                                 reactants={a:1}, products={l:1})
        r_d1 = gillespy2.Reaction(name="decay1", rate=k_d,
                                 reactants={d:1}, products={a:1})
        r_d2 = gillespy2.Reaction(name="decay2", rate=k_d,
                                 reactants={l:1}, products={a:1})
        
        self.add_reaction([r_a1, r_a2 , r_n1, r_n2 , r_d1 , r_d2 ])

        # Set the timespan for the simulation.
        self.timespan(np.linspace(0, 100, 101))

w_array= np.array([])   # the array to store the values for w in

model = chirality()



for run in range(samples) :
    
    print(run)
       
    results = model.run(number_of_trajectories=1)
    
    d = results['D'][100]
    l = results['L'][100]
    
    w = (d - l)/ (d + l)
    
    
    w_array = np.append(w_array , np.array([w]))

#------------------- Figure Settings ------------------------------------------#

fig, ax1 = plt.subplots()
ax1.set(xlabel = r'$\omega$' , ylabel = r'P($\omega$)' , title = 'probability density from simulation v.s. analytic model \n ' r'$\alpha$ = 1')


#------------------ Analytical Curve -----------------------------------------#

w_axis = np.linspace(-0.99,1,200 , endpoint = False)
p_analytic = (1 - w_axis ** 2) ** (k_n_value/k_a_value - 1)

ax1.plot(w_axis , p_analytic/np.sum(p_analytic)/0.01 , '--' , label = 'analytic result' , color = 'cornflowerblue')
ax1.legend()

#------------------- Simulation ----------------------------------------------#

p_sim = np.histogram(w_array , bins = 20 , range = (-1 , 1))


ax1.plot(p_sim[1][:20] , p_sim[0] /np.sum(p_sim[0]) /0.1 , 'x' , ms = 5 , label = 'simulation result' , color = 'crimson')
ax1.legend()

plt.savefig("alpha ", dpi = 500)
plt.show()

