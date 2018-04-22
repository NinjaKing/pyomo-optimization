
# coding: utf-8

# In[1]:


from pyomo.environ import *


# In[2]:


model = AbstractModel()


# In[3]:


## Define sets ##
# H+ : departure depots
model.H_plus = Set(doc='Departure depots')
# H- : arrival depots
model.H_minus = Set(doc='Arrival depots')

# W+ : pickup sites for repairs
model.W_plus = Set(doc='Pickup sites for repairs')
# W+ : repair sites
model.W_minus = Set(doc='Repair sites')

# vehicles
model.V = Set(doc='Vehicles')

model.S_plus = model.H_plus | model.W_plus | model.W_minus
model.S_minus = model.H_minus | model.W_plus | model.W_minus

model.S = model.H_plus| model.H_minus | model.W_plus | model.W_minus

# C : the precedence constraints from the ROP solution
#     the odered (i, j) stand for: repair site j after site i
model.C = Set(within=model.W_minus*model.W_minus)


# In[4]:


## Define parameters ##
# h_plus : departure depots of vehicle v
model.h_plus = Param(model.V, within=model.H_plus)
# h_minus : arrival depots of vehicle v
model.h_minus = Param(model.V, within=model.H_minus)

# c : capacity of vehicle v
model.c = Param(model.V, within=NonNegativeReals)

# s : service time at site i in W-
model.s = Param(model.W_minus, within=NonNegativeReals)

# d : pickup load at site i in W+
model.d = Param(model.W_plus, within=NonNegativeReals)

# t : travel time between sites i and j
model.t = Param(model.S, model.S, within=NonNegativeReals)

# r : repair site associated with pickup site i
model.r = Param(model.W_plus, within=model.W_minus)


# In[5]:


## Define variables ##
# successor of site i
model.sigma = Var(model.S, within=model.S_minus)

# vehicle : vehicle visiting site i
model.vehicle = Var(model.S, within=model.V)

# weight : load of the vehicle visiting site i at i
model.weight = Var(model.S, within=NonNegativeReals)

# eat : earliest arrival time at site i
model.eat = Var(model.S, within=NonNegativeReals)


# In[6]:


## Define objective ##
# maximize the active load
def obj_rule(model):
    return sum(model.eat[i] + model.s[i] for i in model.W_minus)
model.obj = Objective(rule=obj_rule, sense=minimize)


# In[7]:


## Define constraint ##
def m3_2_rule(model, i, j):
    if i is j:
        return Constraint.Skip
    return model.sigma[i] != model.sigma[j]
model.m3_2 = Constraint(model.S_plus, model.S_plus, rule=m3_2_rule)

def m3_3_rule(model, v):
    return model.vehicle[model.h_plus[v]] == v
model.m3_3 = Constraint(model.V, rule=m3_3_rule)

def m3_4_rule(model, v):
    return model.vehicle[model.h_minus[v]] == v
model.m3_4 = Constraint(model.V, rule=m3_4_rule)

def m3_5_rule(model, i):
    return model.vehicle[model.sigma[i]] == model.vehicle[i]
model.m3_5 = Constraint(model.S_plus, rule=m3_5_rule)

def m3_6_rule(model, i):
    return model.vehicle[model.r[i]] == model.vehicle[i]
model.m3_6 = Constraint(model.W_plus, rule=m3_6_rule)

def m3_7_rule(model, i):
    return model.weight[i] == 0
model.m3_7 = Constraint(model.H_plus, rule=m3_7_rule)

def m3_8_rule(model, i):
    return model.weight[model.sigma[i]] == model.weight[i] + model.d[i]
model.m3_8 = Constraint(model.W_plus, rule=m3_8_rule)

#def m3_9_rule(model, i):
#    return model.weight[model.sigma[i]] == model.weight[i] - model.d[i]
#model.m3_9 = Constraint(model.W_minus, rule=m3_9_rule)
def m3_9_rule(model, i):
    return model.weight[model.sigma[model.r[i]]] == model.weight[model.r[i]] - model.d[i]
model.m3_9 = Constraint(model.W_plus, rule=m3_9_rule)

def m3_10_rule(model, i):
    return model.weight[i] <= model.c[model.vehicle[i]]
model.m3_10 = Constraint(model.S_minus, rule=m3_10_rule)

def m3_11_rule(model, i):
    return model.eat[i] == 0
model.m3_11 = Constraint(model.H_minus, rule=m3_11_rule)

def m3_12_rule(model, i):
    return model.eat[i] + model.t[i, model.sigma[i]] <= model.eat[model.sigma[i]]
model.m3_12 = Constraint(model.H_minus | model.W_plus, rule=m3_12_rule)

def m3_13_rule(model, i):
    return model.eat[i] + model.s[i] + model.t[i, model.sigma[i]] <= model.eat[model.sigma[i]]
model.m3_13 = Constraint(model.W_minus, rule=m3_13_rule)

def m3_14_rule(model, i, j):
    return model.eat[i] <= model.eat[j]
model.m3_14 = Constraint(model.C, rule=m3_14_rule)
