import numpy as np
import matplotlib.pyplot as plt
#from tqdm.notebook import tqdm
from scipy.optimize import curve_fit

# %%
class Simbox():
    
    def __init__(self, N_part, N_steps, filename, density, lambd, mean_sigma):
        self.N = N_part
        self.steps = N_steps
        self.file = filename
        self.dt = 1e-4
        self.mean_sigma = mean_sigma
        self.dens = density
        self.lambd = lambd
        
    def distance_MSD(self):
        delta_x = self.x1 - self.x2
        delta_y = self.y1 - self.y2
        delta_z = self.z1 - self.z2
        sqdist = delta_x**2 + delta_y**2 + delta_z**2 
        return sqdist
            
    def get_values(self, step):
        x, y, z, sigma = np.loadtxt(self.file, skiprows=step*521+9,
                                    usecols=(2, 3, 4, 5), max_rows=512, unpack=True)
        return x, y, z, sigma
        
    def MSD(self): # function that calculates MSD and ACF
        self.time = []
        self.msd = [] # mean squared displacement
        self.acf = [] # autocorrelation function
        for i in range (self.steps):
            self.time.append(i * self.dt * 1000)
            N_max = self.steps - i
            sqdist_mean1 = np.zeros(N_max)
            ACF_mean1 = np.zeros(N_max)
            if i%int(self.steps/100) == 0:
                print('{0}% of calculation'.format(int(100 * i / self.steps)))
            
            for j in range (N_max):
                self.x1, self.y1, self.z1, self.sigma1 = self.get_values(j) # get x, y, z, sigma at step t0 = j
                self.x2, self.y2, self.z2, self.sigma2 = self.get_values(j + i) # get values at step t0 + t = j+i
                sqdist_mean1[j] = np.mean(self.distance_MSD()) # mean over all particles of r^2
                ACF_mean1[j] = np.mean((self.sigma1 - self.mean_sigma) * (self.sigma2 - self.mean_sigma))
                
            self.msd.append(np.mean(sqdist_mean1)) # mean over all times t0
            self.acf.append(np.mean(ACF_mean1))
            self.ouput_prop('MSD', self.msd)
            self.ouput_prop('ACF', self.acf)
        return self.time, self.msd , self.acf  
    
    def ouput_prop(self, name, prop):
        filename = "{0}_density_{1}_Npart_{2}_lambda_{3}.txt".format(name, self.dens, self.N, self.lambd)
        file = open(filename,'a')
        for i in range (np.size(prop)):
            file.write('{0} {1} \n'.format(self.time[i], prop[i]))
        file.close()
        
# %%
def Gauss(x, a, b, c):
    appr = c * np.exp(- (x - a)**2 / (2 * b**2))
    return appr

# %%
def read_prop(prop_name, dens, N_part, alph, alph_or_lambd):
    filename = "{0}_density_{1}_Npart_{2}_{3}_{4}.txt".format(prop_name, dens, N_part, alph_or_lambd, alph)
    prop_val = np.loadtxt(filename)  
    return prop_val[:, 0], prop_val[:, 1]

#%% 
def fit_N(N_val, sigma_val):
    popt, pcov = curve_fit(Gauss, sigma_val, N_val)  # approximation
    appr = Gauss(sigma_val, *popt)
    return popt[0], popt[1], appr

# %%

lambd = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
N_002 = np.zeros((np.size(lambd), 100))
mean_002 = np.zeros(np.size(lambd))
var_002 = np.zeros(np.size(lambd)) 
appr_002 = np.zeros((np.size(lambd), 100))

N_095 = np.zeros((np.size(lambd), 100))
mean_095 = np.zeros(np.size(lambd))
var_095 = np.zeros(np.size(lambd))
appr_095 = np.zeros((np.size(lambd), 100))

for i in range (np.size(lambd)):
    sigmaval, N_002[i, :] = read_prop('N_sigma', 0.02, 512, lambd[i], 'lambda')
    mean_002[i], var_002[i], appr_002[i, :] = fit_N(N_002[i, :], sigmaval)
    
    sigmaval, N_095[i, :] = read_prop('N_sigma', 0.95, 512, lambd[i], 'lambda')
    mean_095[i], var_095[i], appr_095[i, :] = fit_N(N_095[i, :], sigmaval)

N_save = 9990 # number of saved steps be used
N_part = 512

msd1 = np.zeros((np.size(lambd), N_save))
acf1 = np.zeros((np.size(lambd), N_save))

for i in range (np.size(lambd)):
    Filename1 = "Prod_traj_unwrapped_density_0.02_Npart_512_lambda_{0}.dump".format(lambd[i])
    p1 = Simbox(N_part, N_save, Filename1, 0.02, lambd[i], mean_002[i])
    t, msd1[i, :], acf1[i, :] = p1.MSD()
    print('Simulation {0}\n'.format(i))
    
for i in range (np.size(lambd)):
    plt.plot(t, msd1[i, :],label=r'$T_\sigma$ = {0}'.format(lambd[i]))
    plt.xlabel('time, [s]')
    plt.ylabel('MSD')
    plt.legend()
plt.show()

for i in range (np.size(lambd)):
    plt.plot(t, acf1[i, :],label=r'$T_\sigma$ = {0}'.format(lambd[i]))
    plt.xlabel('time, [s]')
    plt.ylabel('ACF')
    plt.legend()
plt.show()

msd2 = np.zeros((np.size(lambd), N_save))
acf2 = np.zeros((np.size(lambd), N_save))

for i in range (np.size(lambd)):
    Filename1 = "Prod_traj_unwrapped_density_0.95_Npart_512_lambda_{0}.dump".format(lambd[i])
    p2 = Simbox(N_part, N_save, Filename1, 0.95, lambd[i], mean_095[i])
    t, msd2[i, :] , acf2[i, :] = p2.MSD()
    print('Simulation {0}\n'.format(i))
    
for i in range (np.size(lambd)):
    plt.plot(t, msd2[i, :],label=r'$T_\sigma$ = {0}'.format(lambd[i]))
    plt.xlabel('time, [s]')
    plt.ylabel('MSD')
    plt.legend()
plt.show()

for i in range (np.size(lambd)):
    plt.plot(t, acf2[i, :],label=r'$T_\sigma$ = {0}'.format(lambd[i]))
    plt.xlabel('time, [s]')
    plt.ylabel('ACF')
    plt.legend()
plt.show()