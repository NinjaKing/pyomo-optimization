
# coding: utf-8

# In[1]:


from pyomo.environ import *


# In[2]:


model = AbstractModel()


# In[3]:


## Define sets ##
# N : buses
model.N = Set(doc='Buses')
# L : lines
model.L = Set(within=model.N*model.N, doc='Lines')
# R : set of lines to repair
model.R = Set(within=model.L, doc='Set of lines to repair')


# In[4]:


# K : set of [1,..,siz(R)]
model.K = RangeSet(len(model.R))


# In[5]:


## Define parameters ##
# p_l : desired active loads
model.p_l = Param(model.N)
# q_l : desired reactive loads
model.q_l = Param(model.N)

# g, b : admittance between bus n and bus m
model.g = Param(model.L) # conductance
model.b = Param(model.L) # susceptance

# S : line loading limit on line (n, m)
model.S = Param(model.L)

# interval of active generation
model.p_g_lb = Param(model.N) # lower bound
model.p_g_ub = Param(model.N) # upper bound

# interval of reactive generation
model.q_g_lb = Param(model.N) # lower bound
model.q_g_ub = Param(model.N) # upper bound

# interval of voltage magnitude
model.v_lb = Param(model.N) # lower bound
model.v_ub = Param(model.N) # upper bound


# In[6]:


## Define variables ##
# o : line (n,m) is repaired at step k
model.o = Var(model.R, model.K, domain=Binary, doc='line (n,m) is repaired at step k')

# z : line (n,m) is activated at step k
model.z = Var(model.L, model.K, domain=Binary, doc='line (n,m) is activated at step k')

# theta : phase angle of bus n (radians)
model.theta = Var(model.N, model.K)

# v : voltage magnitude at step k
def v_bound(model, n, k):
   return (model.v_lb[n], model.v_ub[n])
model.v = Var(model.N, model.K, bounds=v_bound)

# p_g : active generation at step k
def p_g_bound(model, n, k):
   return (model.p_g_lb[n], model.p_g_ub[n])
model.p_g = Var(model.N, model.K, bounds=p_g_bound)

# q_g : reactive generation at step k
def q_g_bound(model, n, k):
   return (model.q_g_lb[n], model.q_g_ub[n])
model.q_g = Var(model.N, model.K, bounds=q_g_bound)

# l : percentage load served  at step k
model.l = Var(model.N, model.K, bounds=(0, 1))

# p : active power on line (n, m) at step k
def p_bound(model, n, m, k):
   return (- model.S[n,m], model.S[n,m])
model.p = Var(model.L, model.K, bounds=p_bound)

# q : reactive power on line (n, m) at step k
def q_bound(model, n, m, k):
   return (- model.S[n,m], model.S[n,m])
model.q = Var(model.L, model.K, bounds=q_bound)


# In[7]:


## Define objective ##
# maximize the active load
def obj_rule(model):
    return sum(sum(model.p_l[n] * model.l[n,k] for n in model.N) for k in model.K)
model.obj = Objective(rule=obj_rule, sense=maximize)


# In[8]:


## Define constraint ##
def m2_2_rule(model, k):
    return sum(model.o[n,m,k] for (n,m) in model.R) == k
model.m2_2 = Constraint(model.K, rule=m2_2_rule)

def m2_3_rule(model, k, n, m):
    if k == 1:
        return Constraint.Skip
    return model.o[n,m,k-1] <= model.o[n,m,k]
model.m2_3 = Constraint(model.K, model.R, rule=m2_3_rule)

def m2_4_rule(model, k, n, m):
    return model.z[n,m,k] <= model.o[n,m,k]
model.m2_4 = Constraint(model.K, model.R, rule=m2_4_rule)

def m2_5_rule(model, k, n, m):
    return model.z[n,m,k] == 1
model.m2_5 = Constraint(model.K, model.L - model.R, rule=m2_5_rule)

def m2_6_rule(model, k, n):
    return model.p_g[n,k] - model.p_l[n] * model.l[n,k] <= sum(model.p[i,j,k] for (i,j) in model.L if i == n)
model.m2_6 = Constraint(model.K, model.N, rule=m2_6_rule)

def m2_7_rule(model, k, n):
    return model.q_g[n,k] - model.q_l[n] * model.l[n,k] <= sum(model.q[i,j,k] for (i,j) in model.L if i == n)
model.m2_7 = Constraint(model.K, model.N, rule=m2_7_rule)

def m2_8_rule(model, k, n, m):
    if model.z[n,m,k] == 0:
        return Constraint.Skip
    admit = model.v[n,k]**2 * model.g[n,m]         - model.v[n,k] * model.v[m,k] * model.g[n,m] * cos(model.theta[n,k] - model.theta[m,k])         - model.v[n,k] * model.v[m,k] * model.b[n,m] * sin(model.theta[n,k] - model.theta[m,k])
        
    return model.p[n,m,k] == admit 
model.m2_8 = Constraint(model.K, model.L, rule=m2_8_rule)

def m2_9_rule(model, k, n, m):
    if model.z[n,m,k] == 0:
        return Constraint.Skip
    admit = - model.v[n,k]**2 * model.b[n,m]         + model.v[n,k] * model.v[m,k] * model.b[n,m] * cos(model.theta[n,k] - model.theta[m,k])         - model.v[n,k] * model.v[m,k] * model.g[n,m] * sin(model.theta[n,k] - model.theta[m,k])
        
    return model.q[n,m] == admit 
model.m2_9 = Constraint(model.K, model.L, rule=m2_9_rule)

def m2_10_rule(model, k, n, m):
    if model.z[n,m,k] == 1:
        return Constraint.Skip
    return model.p[n,m,k] == 0 and model.q[n,m,k] == 0
model.m2_10 = Constraint(model.K, model.L, rule=m2_10_rule)

def m2_11_rule(model, k, n, m):
    return model.p[n,m,k]**2 + model.q[n,m,k]**2 <= model.S[n,m]**2
model.m2_11 = Constraint(model.K, model.L, rule=m2_11_rule)

