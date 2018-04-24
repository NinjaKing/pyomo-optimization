
from pyomo.environ import *

M = 1e100


model = AbstractModel()


## Define sets ## ========================================================
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


## Define parameters ## ========================================================
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

# extra param to define load at repair site
def d_r_init(model):
    val = {model.r[i]:model.d[i] for i in model.W_plus}
    return val
model.d_r = Param(model.W_minus, within=NonNegativeReals, initialize=d_r_init)


## Define variables ## ========================================================
# successor of site i
#model.sigma = Var(model.S, within=model.S_minus)
# binary variable: sigma_i_j = 1 when site j is successor of site i, = 0 otherwise
model.sigma = Var(model.S_plus, model.S_minus, domain=Binary)

# vehicle : vehicle visiting site i
#model.vehicle = Var(model.S, within=model.V)
# binary variable: vehicle_v_i = 1 if v visit i, = 0 otherwise
model.vehicle = Var(model.V, model.S, domain=Binary)

# weight : load of the vehicle visiting site i at i
model.weight = Var(model.S, within=NonNegativeReals)

# eat : earliest arrival time at site i
model.eat = Var(model.S, within=NonNegativeReals)


## Define objective ## ========================================================
# maximize the active load
def obj_rule(model):
    return sum(model.eat[i] + model.s[i] for i in model.W_minus)
model.obj = Objective(rule=obj_rule, sense=minimize)


## Define constraint ## ========================================================
# constrain for binary variable: all_different
def m3_2_c1_rule(model, j):
    return sum(model.sigma[i,j] for i in model.S_plus) == 1
model.m3_2_c1 = Constraint(model.S_minus, rule=m3_2_c1_rule)

def m3_2_c2_rule(model, i):
    return sum(model.sigma[i,j] for j in model.S_minus) == 1
model.m3_2_c2 = Constraint(model.S_plus, rule=m3_2_c2_rule)

# M3.3
def m3_3_rule(model, v):
    return model.vehicle[v, model.h_plus[v]] == 1
model.m3_3 = Constraint(model.V, rule=m3_3_rule)

# M3.4
def m3_4_rule(model, v):
    return model.vehicle[v, model.h_minus[v]] == 1
model.m3_4 = Constraint(model.V, rule=m3_4_rule)

# each site will only have one vehicle visit it
def vehicle_cons_rule(model, i):
    return sum(model.vehicle[v, i] for v in model.V) == 1
model.vehicle_cons = Constraint(model.W_plus | model.W_minus, rule=vehicle_cons_rule)
    
# M3.5
def m3_5_1_rule(model, i, j, v):
    return model.vehicle[v, i] - model.vehicle[v, j] <= M * (1 - model.sigma[i, j])
model.m3_5_1 = Constraint(model.S_plus, model.S_minus, model.V, rule=m3_5_1_rule)

def m3_5_2_rule(model, i, j, v):
    return - model.vehicle[v, i] + model.vehicle[v, j] <= M * (1 - model.sigma[i, j])
model.m3_5_2 = Constraint(model.S_plus, model.S_minus, model.V, rule=m3_5_2_rule)

# M3.6
def m3_6_rule(model, i, v):
    return model.vehicle[v, model.r[i]] == model.vehicle[v, i]
model.m3_6 = Constraint(model.W_plus, model.V, rule=m3_6_rule)

# M3.7
def m3_7_rule(model, i):
    return model.weight[i] == 0
model.m3_7 = Constraint(model.H_plus, rule=m3_7_rule)

# M3.8
def m3_8_1_rule(model, i, j):
    return model.weight[j] - model.weight[i] - model.d[i] <= M * (1 - model.sigma[i, j])
model.m3_8_1 = Constraint(model.W_plus, model.S_minus, rule=m3_8_1_rule)

def m3_8_2_rule(model, i, j):
    return - model.weight[j] + model.weight[i] + model.d[i] <= M * (1 - model.sigma[i, j])
model.m3_8_2 = Constraint(model.W_plus, model.S_minus, rule=m3_8_2_rule)

# M3.9
def m3_9_1_rule(model, i, j):
    return model.weight[j] - model.weight[i] + model.d_r[i] <= M * (1 - model.sigma[i, j])
model.m3_9_1 = Constraint(model.W_minus, model.S_minus, rule=m3_9_1_rule)

def m3_9_2_rule(model, i, j):
    return - model.weight[j] + model.weight[i] - model.d_r[i] <= M * (1 - model.sigma[i, j])
model.m3_9_2 = Constraint(model.W_minus, model.S_minus, rule=m3_9_2_rule)

# M3.10
def m3_10_rule(model, i, v):
    return model.weight[i] <= model.c[v] + M * (1 - model.vehicle[v, i])
model.m3_10 = Constraint(model.S_minus, model.V, rule=m3_10_rule)

# M3.11
def m3_11_rule(model, i):
    return model.eat[i] == 0
model.m3_11 = Constraint(model.H_plus, rule=m3_11_rule)

# M3.12
def m3_12_rule(model, i, j):
    return model.eat[i] + model.t[i, j] <= model.eat[j] + M * (1 - model.sigma[i, j])
model.m3_12 = Constraint(model.H_plus | model.W_plus, model.S_minus, rule=m3_12_rule)

# M3.13
def m3_13_rule(model, i, j):
    return model.eat[i] + model.s[i] + model.t[i, j] <= model.eat[j] + M * (1 - model.sigma[i, j])
model.m3_13 = Constraint(model.W_minus, model.S_minus, rule=m3_13_rule)

# M3.14
def m3_14_rule(model, i, j):
    return model.eat[i] <= model.eat[j]
model.m3_14 = Constraint(model.C, rule=m3_14_rule)