#!/usr/bin/env python
# coding: utf-8

# # Identify the Theodorsen unsteady lift
# 
# ## Data generation
# 
# As for the previous case the data is generated through a python script. Since the objective is to identify the Theodorsen function through its input-output response the classical model for a monochromatic oscillation cannot be used, since a wide spectrum input has to be prefered. The data-generation code uses then a linear approximation of the Theodorsen function obtained by balanced truncation in {cite}`brunton2013empirical`. 
# 
# <span style="color: red;">Here you might fill in with a very short overview of your code. If you want we might add a more complete explanation in an appendix, but I'd leave this part quite lean.</span> 

# In[1]:


import numpy as np
import control
import Theodorsen_control as theodorsen
from signals import prbs

a = 1/2  # pitch axis wrt to 1/2-chord
b = 1.  # half-chord length of the airfoil
# default values of C_1 and C_2 used
airfoil = theodorsen.AirfoilGeometry(a=a, b=b)

# THEDORSEN MODEL

# the balanced truncation Theodorsen function approximation
theodorsen_function_sys = theodorsen.theodorsen_function_balanced_truncation_ss()

# state-space system with both α" and h" as inputs
theodorsen_full_sys = theodorsen.unsteady_lift_ss(
    airfoil, theodorsen_function_sys)

# INPUT SIGNALS

NT = 20000 # number of time steps
t = np.linspace(0, 1, NT)
dt = 400 * t[-1] / NT

amplitude= 0.2

alpha_ddot = prbs(t, dt, min = -amplitude, max = amplitude)
h_ddot = prbs(t, 2 * dt, min = -amplitude, max = amplitude)
u_MISO = np.vstack((h_ddot.T, alpha_ddot.T))

# TIME RESPONSE

X0 = np.array([0,0,0,0,
              0.0,   # h_dot(0)
              0.05, # alpha(0)
              0.0])   # alpha_dot(0)

output = control.forced_response(
    theodorsen_full_sys, T=t, U=u_MISO, X0 = X0)

# TIME RESPONSE POSTPROCESSING

data = theodorsen.TheodorsenTimeResponse(
    output, inputs='both', sys=theodorsen_full_sys)


# ## Model regression
# 
# As stated in <span style="color: red;">section</span> we want to fit a model in the form:
# ```{math}
#     C_L(t) = \Theta(t)\xi
# ```
# where $\Theta$ is a library of linear and nonlinear functions and $\xi$ are the corresponding coefficients. The data obtained in the previous section is divided equally between a training set and a test set, such that

# In[2]:


# Training set
X_train = np.stack([data.alpha_e[:NT // 2],
                    data.alpha_dot[:NT // 2]] , axis = -1) # physical states of the system
CL_train = data.C_L[:NT // 2] # output of the system

# Test set
X_test = np.stack([data.alpha_e[NT // 2:],
                   data.alpha_dot[NT // 2:]], axis = -1)
CL_test = data.C_L[NT // 2:]


# Since we know that the possibly nonlinear part of the Theodorsen model lies in the circulatory lift only, and that in general the added mass term behaves linearly, we decided to only fit the first one, assumng as a first moment that the added mass coefficients are known.

# In[3]:


ttt = theodorsen_full_sys.D @ u_MISO[:,NT // 2:]


# In[4]:


CL_train -= (theodorsen_full_sys.D @ u_MISO[:,NT // 2:]).T
CL_test -= (theodorsen_full_sys.D @ u_MISO[:, :NT // 2]).T


# Since SINDy is based on the machine learning library `sklearn` it is possible to substitute the derivatives of the classical SINDy models with the $C_l$ without any modification of the underlying code. The set up for the regression is then:

# In[5]:


import pysindy as ps

expansionDegree = 4

# SETTING UP SINDy 
optimizer = ps.optimizers.stlsq.STLSQ(threshold = 0.1, alpha = 1e-06, max_iter = 50)
library = ps.feature_library.polynomial_library.PolynomialLibrary(degree = expansionDegree) # standard polynomial library

model = ps.SINDy(optimizer = optimizer, 
			     feature_library = library,
			     feature_names = ['alpha_e', 'alpha_dot']) 
t = data.t
classicalModel = model.fit(X_train, t = t[1] - t[0], x_dot = CL_train)
classicalModel.print()


# The results for the classical polynomial library can be compared with those for an orthogonal and othonormal base

# In[6]:


from PolynomialChaos import *

aPC = PolynomialChaos(
     X_train,
     expansionDegree = expansionDegree,
     numberOfInputs = 2)
aPC.ComputeCoefficients(threshold =0.01, normalize = False)
coefficients = aPC.coefficients
AlphaMatrix = aPC.AlphaMatrix
 
LibraryList = GenerateLibraryList(
     expansionDegree=expansionDegree,
     coefficients = coefficients,
     AlphaMatrix = AlphaMatrix)

aPC_norm = PolynomialChaos(
     X_train,
     expansionDegree = expansionDegree,
     numberOfInputs = 2)
aPC_norm.ComputeCoefficients(threshold =0.01, normalize = True)
coefficients_norm = aPC_norm.coefficients
 
LibraryList_norm = GenerateLibraryList(
     expansionDegree=expansionDegree,
     coefficients = coefficients_norm,
     AlphaMatrix = AlphaMatrix)

optimizer = ps.optimizers.stlsq.STLSQ(threshold = 0.1, alpha = 1e-06, max_iter = 50)
library_orth = ps.feature_library.custom_library.CustomLibrary(LibraryList)
library_orthnorm = ps.feature_library.custom_library.CustomLibrary(LibraryList_norm)

model = ps.SINDy(optimizer = optimizer, 
			     feature_library = library_orth,
			     feature_names = ['alpha_e', 'alpha_dot']) 
OrthModel = model.fit(X_train, t = t[1] - t[0], x_dot = CL_train) # Orthogonal model

model = ps.SINDy(optimizer = optimizer, 
			     feature_library = library_orthnorm,
			     feature_names = ['alpha_e', 'alpha_dot']) 
OrthNormModel = model.fit(X_train, t = t[1] - t[0], x_dot = CL_train) # Orthonormal model

print('-------- Orthogonal  model --------\n')
OrthModel.print()
print('\n-------- Orthonormal model --------\n')
OrthNormModel.print()


# In[ ]:




