from pyomo.environ import *
import numpy as np


model = AbstractModel()


## Define sets ##
# N : buses
model.N = Set(doc='Buses')
# L : lines
model.L = Set(within=model.N*model.N, doc='Lines')


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


## Define variables ##
# theta : phase angle of bus n (radians)
model.theta = Var(model.N)

# v : voltage magnitude
model.v = Var(model.N)

# p_g : active generation
model.p_g = Var(model.N)

# q_g : reactive generation
model.q_g = Var(model.N)

# l : percentage load served 
model.l = Var(model.N, bounds=(0, 1))

# p : active power on line (n, m)
model.p = Var(model.L)

# q : reactive power on line (n, m)
model.q = Var(model.L)


## Define objective ##
# maximize the active load
def obj_rule(model):
    return sum(model.p_l[n] * model.l[n] for n in model.N)
model.obj = Objective(rule=obj_rule, sense=maximize)


## Define constraint ##
# active generation limit
def p_g_bound_rule(model, n):
    return model.p_g_lb[n] <= model.p_g[n] <= model.p_g_ub[n]
model.p_g_limit = Constraint(model.N, rule=p_g_bound_rule)

# reactive generation limit
def q_g_bound_rule(model, n):
    return model.q_g_lb[n] <= model.q_g[n] <= model.q_g_ub[n]
model.q_g_limit = Constraint(model.N, rule=q_g_bound_rule)

# voltage magnitude limit
def v_bound_rule(model, n):
    return model.v_lb[n] <= model.v[n] <= model.v_ub[n]
model.v_limit = Constraint(model.N, rule=v_bound_rule)

# active power bound
def p_bound_rule(model, n, m):
    return - model.S[n,m] <= model.p[n,m] <= model.S[n,m]
model.p_limit = Constraint(model.L, rule=p_bound_rule)

# reactive power bound
def q_bound_rule(model, n, m):
    return - model.S[n,m] <= model.q[n,m] <= model.S[n,m]
model.q_limit = Constraint(model.L, rule=q_bound_rule)


# flow conservation
def flow_p_rule(model, n):
    return model.p_g[n] - model.p_l[n] * model.l[n] <= sum(model.p[i,j] for (i,j) in model.L)
model.flow_p = Constraint(model.N, rule=flow_p_rule)

def flow_q_rule(model, n):
    return model.q_g[n] - model.q_l[n] * model.l[n] <= sum(model.q[i,j] for (i,j) in model.L)
model.flow_q = Constraint(model.N, rule=flow_q_rule)


# the real and reactive power flows on lines
def active_power_flow_rule(model, n, m):
    admit = model.v[n]**2 * model.g[n,m]         - model.v[n] * model.v[m] * model.g[n,m] * cos(model.theta[n] - model.theta[m])         - model.v[n] * model.v[m] * model.b[n,m] * sin(model.theta[n] - model.theta[m])
        
    return model.p[n,m] == admit
model.active_power_flow = Constraint(model.L, rule=active_power_flow_rule)

def reactive_power_flow_rule(model, n, m):
    admit = - model.v[n]**2 * model.b[n,m]         + model.v[n] * model.v[m] * model.b[n,m] * cos(model.theta[n] - model.theta[m])         - model.v[n] * model.v[m] * model.g[n,m] * sin(model.theta[n] - model.theta[m])
        
    return model.q[n,m] == admit
model.reactive_power_flow = Constraint(model.L, rule=reactive_power_flow_rule)


# the thermal limits for lines
def thermal_limit_rule(model, n, m):
    return model.p[n,m]**2 + model.q[n,m]**2 <= model.S[n,m]**2
model.thermal_limit = Constraint(model.L, rule=thermal_limit_rule)
