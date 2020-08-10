# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 16:38:51 2018

@author: lenovo
"""

'''
In diesem Code:

##### Vereinfachte Gleichungen für die Entropie und die Enthalpie sind aus MINLP_Optimierung Buch von Mark-Jüdes Anhang A 
'''

from pyomo.environ import *
from iapws import IAPWS97
from pyomo.gdp import *

m = ConcreteModel()

##################################################################

# Koeffizienten für die vereinfachte Stoffwerte für Wasser: 

##################################################################

a_superheated = [2012.24356, 1.70989, 0.00027, -0.00403, 9.37999, -0.00001, -5896.61864]
b_superheated = [5.61195, 0.00535, -1.53838*10**(-6),-0.46858, -0.00066]

a_liquid = [-1.754335*10**3, 9.20571709,-1.353812*10**(-2), 1.2044842403*10**(-5), 0.49891277, 0.010589532]
b_liquid = [1.085, -6.084*10**(-4), -7.302*10**(-4), 1.303, -0.5014]

a_dampfdruck = [119.990238, -156.368914, 76.930767, -16.95968, 1.415874]

a_Tau = [2363.159168, -1.723792, 0.011738, -0.000013]
b_Tau = [23.16267681942, -0.08701195933, 0.00015876507, -0.00000010484]

a_Siede = [-1033.6898, 3.5369, 0.001] 
b_Siede = [-4.655586398, 0.020460732, -0.000012036]


##################################################################

# The following rules will be used in more than one block. 
# To get a better overview, they are declared here at the beginning of the code.

################################################################## 

def Massenbilanz(m):
    return m.m_in == m.m_out 

def Isobar(m):
    return m.p_in == m.p_out 

def Waermestrom_zu(m):
    return m.Q == m.m_in * (m.h_out - m.h_in)

def Waermestrom_ab(m):
    return m.Q == m.m_in * (m.h_in - m.h_out)

# Formel aus dem Buch "MINLP_Optimierung.." von Mark-Jüdes Anhang A: 
def Enthalpie_ueberhitzt_vereinfacht(m):
    return m.h_out == (a_superheated[0] + a_superheated[1] * m.T_out + a_superheated[2] * m.T_out**2 + a_superheated[3] * m.T_out * m.p_out + a_superheated[4] * m.p_out + 
                       a_superheated[5] * m.p_out**2 + a_superheated[6] * m.p_out/m.T_out)

def Entropie_ueberhitzt_vereinfacht(m):
    return m.s_out == (b_superheated[0] + b_superheated[1] * m.T_out + b_superheated[2] * m.T_out**2 + b_superheated[3] * log(m.p_out) + b_superheated[4]  * m.p_out)


##################################################################

# Beginning the implementation of the blocks! 

##################################################################

'''
Block für die Pumpe 
- Annahmen: isentrop
'''
m.P = Block()

# Variables 
m.P.m_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(10.0,200.0))
m.P.m_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(10.0,200.0))

m.P.p_in = Param(default=0.03) # bar 
m.P.p_out = Var(initialize=210.0, within=NonNegativeReals, bounds=(10.0,270.0)) 
m.P.r_p = Var(initialize=9000.0, within=NonNegativeReals, bounds=(1.0,9000.0)) # Druckerhöhung 

m.P.T_in = Param(default=293.15) 
m.P.T_out = Var(initialize=300.0, within=NonNegativeReals, bounds=(292.0,300.0)) 

m.P.s_in = Param(default=IAPWS97(P=value(m.P.p_in)*0.1, T=value(m.P.T_in)).s)
m.P.s_out = Var(initialize=0.354, within=NonNegativeReals, bounds=(0.0, 11.0)) # in kJ/kg

m.P.h_in = Param(default=IAPWS97(P=value(m.P.p_in)*0.1, T=value(m.P.T_in)).h)
m.P.h_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(10.0,1000.0)) # in kJ/kg 

m.P.W  = Var(initialize=1000.0, within=NonNegativeReals, bounds=(0.0,2000.0))  # in kW --> Villeicht im MW umrechnen


# Connector 
m.P.IN = Connector(initialize=
                        {'Druck':       m.P.p_in,
                         'Temperatur':  m.P.T_in,
                         'Massenstrom': m.P.m_in,
                         'Enthalpie':   m.P.h_in,
                         'Entropie':    m.P.s_in})
    
m.P.OUT = Connector(initialize=
                        {'Druck':       m.P.p_out,
                         'Temperatur':  m.P.T_out,
                         'Massenstrom': m.P.m_out,
                         'Enthalpie':   m.P.h_out,
                         'Entropie':    m.P.s_out})

# Rules for the constraints 
def Druckerhoehung(m):
    return m.p_out == m.p_in * m.r_p
def Pumpenleistung(m):
    return m.W == m.m_in * (m.h_out - m.h_in)
def Isentrop(m):
    return m.s_in == m.s_out

# Formel aus dem Buch "MINLP_Optimierung.." von Mark-Jüdes Anhang A: 
def Enthalpie_fluessig_vereinfacht(m): 
    return (m.h_out == a_liquid[0] + a_liquid[1] * m.T_out + a_liquid[2] * m.T_out**2 +
                       a_liquid[3] * m.T_out**3 + a_liquid[4] * m.p_out/10 + a_liquid[5] * (m.p_out/10)**2) 
def Entropie_SiedendeFluessigkeit(m):
    return m.s_out==b_Siede[0] + b_Siede[1] * m.T_out + b_Siede[2] * m.T_out**2
# Constraints
m.P.Massenbilanz = Constraint(rule=Massenbilanz)
m.P.Austrittsdruck = Constraint(rule=Druckerhoehung)
m.P.Leistung = Constraint(rule=Pumpenleistung)
m.P.Isentrop = Constraint(rule=Isentrop)
m.P.Austrittsenthalpie = Constraint(rule=Enthalpie_fluessig_vereinfacht) 

# Folgende Constraint [s_out = f(T_out, p_out)] dient der Berechnung der Austrittstemperatur 
# da isentrop und s_in bekannt, p_out Optimierungsvariable:  
m.P.Austrittsentropie = Constraint(rule=Entropie_SiedendeFluessigkeit) 

'''
Block für den Verdampfer 
- Annahme: isobar
'''
m.EV = Block()

# Variables 
m.EV.m_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))
m.EV.m_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))

m.EV.p_in = Var(initialize=210.0, within=NonNegativeReals, bounds=(1.0,270.0))
m.EV.p_out = Var(initialize=210.0, within=NonNegativeReals, bounds=(1.0,270.0))

m.EV.T_in = Var(initialize=300.0, within=NonNegativeReals, bounds=(292.0, 300.0)) 
m.EV.T_out = Var(initialize=774.0, within=NonNegativeReals, bounds=(293.15,893.15)) 

m.EV.h_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.EV.h_out = Var(initialize=4000.0, within=NonNegativeReals, bounds=(10.0,10000.0))

m.EV.Q = Var(initialize=400000.0, within=NonNegativeReals, bounds=(50000.0,500000.0))

m.EV.s_in = Var(initialize=0.354, within=NonNegativeReals, bounds=(0.0, 11.0))
m.EV.s_out = Var(initialize=2.0, within=NonNegativeReals, bounds=(0.0, 11.0))

# Connector 
m.EV.IN = Connector(initialize=
                        {'Druck':       m.EV.p_in,
                         'Temperatur':  m.EV.T_in,
                         'Massenstrom': m.EV.m_in,
                         'Enthalpie':   m.EV.h_in,
                         'Entropie':    m.EV.s_in})
    
m.EV.OUT = Connector(initialize=
                        {'Druck':       m.EV.p_out,
                         'Temperatur':  m.EV.T_out,
                         'Massenstrom': m.EV.m_out,
                         'Enthalpie':   m.EV.h_out,
                         'Entropie':    m.EV.s_out})

# Connector constraint 
m.EV.connectIN = Constraint(expr=m.EV.IN==m.P.OUT)

# Constraints
m.EV.Massenbilanz = Constraint(rule=Massenbilanz)
m.EV.Druck = Constraint(rule=Isobar)
m.EV.Waermestrom = Constraint(rule=Waermestrom_zu)
m.EV.Austrittstemperatur = Constraint(expr=m.EV.T_in <= m.EV.T_out)

# Der Zustand 3, also HDT-Eintritt liegt auf jeden Fall im überhitzten Bereich, daher die Regel für die Berechnung der Enthalpie: 
m.EV.Austrittsenthalpie = Constraint(rule=Enthalpie_ueberhitzt_vereinfacht)
m.EV.Austrittsentropie = Constraint(rule=Entropie_ueberhitzt_vereinfacht)


'''
Blocks für die Dampfturbinen
'''

##################################################################

# The following rules will be used for both turbines. 
# To get a better overview, they are defined here.

################################################################## 

def Bruttoarbeit(m):
    return m.W == m.m_in * (m.h_in - m.h_out)

def IsentroperWirkungsgrad(m):
    return m.h_out == m.eta * (m.h_s - m.h_in) + m.h_in

##################################################################
    
def _d_real(disjunct, flag):
    m = disjunct.parent_block()
    if flag:
        disjunct.c1 = Constraint(expr=m.h_out == (a_superheated[0] + a_superheated[1] * m.T_out + a_superheated[2] * m.T_out**2 + a_superheated[3] * m.T_out * m.p_out + 
           a_superheated[4] * m.p_out + a_superheated[5] * m.p_out**2 + a_superheated[6] * m.p_out/m.T_out))
        disjunct.c2 = Constraint(expr=m.s_out == (b_superheated[0] + b_superheated[1] * m.T_out + b_superheated[2] * m.T_out**2 + b_superheated[3] * log(m.p_out) +
           b_superheated[4]  * m.p_out))
    else:
        disjunct.c1 = Constraint(expr=0 == m.s_Tau - (b_Tau[0] + b_Tau[1] * m.T_out + b_Tau[2] * m.T_out**2 + b_Tau[3] * m.T_out**3))
        disjunct.c2 = Constraint(expr=0 == m.s_Siede - (b_Siede[0] + b_Siede[1] * m.T_out + b_Siede[2] * m.T_out**2))
        disjunct.c3 = Constraint(expr=0 == m.h_Tau - (a_Tau[0] + a_Tau[1] * m.T_out + a_Tau[2] * m.T_out**2 + a_Tau[3] * m.T_out**3))
        disjunct.c4 = Constraint(expr=0 == m.h_Siede - (a_Siede[0] + a_Siede[1] * m.T_out + a_Siede[2] * m.T_out**2)   )
        disjunct.c5 = Constraint(expr=0 == m.x_Dampf - (m.h_out - m.h_Tau)/(m.h_Siede - m.h_Tau))
        disjunct.c6 = Constraint(expr=0 == m.x_Dampf - (m.s_out - m.s_Tau)/(m.s_Siede - m.s_Tau))
        disjunct.c7 = Constraint(expr=0 == m.p_out - (a_dampfdruck[0] + a_dampfdruck[1] * m.T_out/100 + a_dampfdruck[2] * (m.T_out/100)**2 + 
           a_dampfdruck[3] * (m.T_out/100)**3 + a_dampfdruck[4] * (m.T_out/100)**4))

def _d_isentrop(disjunct, flag):
    m = disjunct.parent_block()
    if flag:
        disjunct.c1 = Constraint(expr=m.h_s == (a_superheated[0] + a_superheated[1] * m.T_s + a_superheated[2] * m.T_s**2 + a_superheated[3] * m.T_s * m.p_out + 
           a_superheated[4] * m.p_out + a_superheated[5] * m.p_out**2 + a_superheated[6] * m.p_out/m.T_s))
        disjunct.c2 = Constraint(expr=m.s_s == (b_superheated[0] + b_superheated[1] * m.T_s + b_superheated[2] * m.T_s**2 + b_superheated[3] * log(m.p_out) +
           b_superheated[4]  * m.p_out))
        disjunct.c3 = Constraint(expr=m.s_in == m.s_s)
    else:
        disjunct.c1 = Constraint(expr=0 == m.s_Tau_s - (b_Tau[0] + b_Tau[1] * m.T_s + b_Tau[2] * m.T_s**2 + b_Tau[3] * m.T_s**3))
        disjunct.c2 = Constraint(expr=0 == m.s_Siede_s - (b_Siede[0] + b_Siede[1] * m.T_s + b_Siede[2] * m.T_s**2))
        disjunct.c3 = Constraint(expr=0 == m.h_Tau_s - (a_Tau[0] + a_Tau[1] * m.T_s + a_Tau[2] * m.T_s**2 + a_Tau[3] * m.T_s**3))
        disjunct.c4 = Constraint(expr=0 == m.h_Siede_s - (a_Siede[0] + a_Siede[1] * m.T_s + a_Siede[2] * m.T_s**2)  )
        disjunct.c5 = Constraint(expr=0 == m.x_Dampf_s - (m.h_s - m.h_Tau_s)/(m.h_Siede_s - m.h_Tau_s))
        disjunct.c6 = Constraint(expr=0 == m.x_Dampf_s - (m.s_s - m.s_Tau_s)/(m.s_Siede_s - m.s_Tau_s))
        disjunct.c7 = Constraint(expr=0 == m.p_out - (a_dampfdruck[0] + a_dampfdruck[1] * m.T_s/100 + a_dampfdruck[2] * (m.T_s/100)**2 + 
          a_dampfdruck[3] * (m.T_s/100)**3 + a_dampfdruck[4] * (m.T_s/100)**4))
        disjunct.c8 = Constraint(expr=m.s_s == m.s_in)

# Define the disjunction
def _c1(model):
    return [model.d1[0], model.d1[1]]

def _c2(model):
    return [model.d2[0], model.d2[1]]


'''
Block für die Hochdruckturbine
'''

m.HDT = Block()

# Variables 
m.HDT.W = Var(initialize=50600.0, within=NonNegativeReals, bounds=(1000.0,100000.0))   # Als Parameter 

m.HDT.eta = Param(default=0.9)

m.HDT.m_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))
m.HDT.m_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))

m.HDT.p_in = Var(initialize=210.0, within=NonNegativeReals, bounds=(0.03,270.0)) 
m.HDT.p_out = Var(initialize=38.0, within=NonNegativeReals, bounds=(0.03,270.0)) 
#m.HDT.p_s = Var(initialize=0.03, within=NonNegativeReals, bounds=(0.03,270.0))       # Druck bei einer isentropen Entspannung ist gleich den Druck bei der realen Entspannung

m.HDT.T_in = Var(initialize=774.0, within=NonNegativeReals, bounds=(293.15,893.15)) 
m.HDT.T_out = Var(initialize=600.0, within=NonNegativeReals, bounds=(500.0, 893.15))
m.HDT.T_s = Var(initialize=550.0)#, within=NonNegativeReals, bounds=(292.0, 893.15))       # Austrittstemperatur bei einer isentropen Entspannung = T_out wenn das im NDG ist 

m.HDT.h_in = Var(initialize=4000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.HDT.h_out = Var(initialize=3000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.HDT.h_s = Var(initialize=2500.0, within=NonNegativeReals, bounds=(10.0,10000.0))       # Austrittsenthalpie bei einer isentropen Entspannung


m.HDT.h_Tau = Var(initialize=2900.0)#, within=NonNegativeReals, bounds=(10.0,10000.0))     # Enthalpie gesättigter Dampf
m.HDT.h_Siede = Var(initialize=1000.0)#, within=NonNegativeReals, bounds=(10.0,10000.0))   # Enthalpie siedender Flüssigkeit 
m.HDT.h_Siede_s = Var(initialize=2900.0)#, within=NonNegativeReals, bounds=(10.0,10000.0)) 
m.HDT.h_Tau_s = Var(initialize=1000.0)#, within=NonNegativeReals, bounds=(10.0,10000.0)) 

m.HDT.s_in = Var(initialize=2.0, within=NonNegativeReals, bounds=(1.0, 11.0))
m.HDT.s_out = Var(initialize=5.0, within=NonNegativeReals, bounds=(1.0, 11.0))
m.HDT.s_s = Var(initialize=2.0, within=NonNegativeReals, bounds=(1.0, 11.0))


m.HDT.s_Siede = Var(initialize=7.0)#, within=NonNegativeReals, bounds=(1.0, 11.0))
m.HDT.s_Tau = Var(initialize=3.0)#, within=NonNegativeReals, bounds=(1.0, 11.0))
m.HDT.s_Siede_s = Var(initialize=7.0)#, within=NonNegativeReals, bounds=(1.0, 11.0))
m.HDT.s_Tau_s = Var(initialize=3.0)#, within=NonNegativeReals, bounds=(1.0, 11.0))

m.HDT.x_Dampf = Var(initialize=0.85, within=NonNegativeReals, bounds=(0.85,1.0))   # größer als 0.85, zu hoher Endnässe ist schlecht für die Turbinenschaufeln
m.HDT.x_Dampf_s = Var(initialize=0.85, within=NonNegativeReals, bounds=(0.3,1.0))

# Connector 
m.HDT.IN = Connector(initialize=
                        {'Druck':       m.HDT.p_in,
                         'Temperatur':  m.HDT.T_in,
                         'Massenstrom': m.HDT.m_in,
                         'Enthalpie':   m.HDT.h_in,
                         'Entropie':    m.HDT.s_in})
    
m.HDT.OUT = Connector(initialize=
                        {'Druck':       m.HDT.p_out,
                         'Temperatur':  m.HDT.T_out,
                         'Massenstrom': m.HDT.m_out,
                         'Enthalpie':   m.HDT.h_out,
                         'Entropie':    m.HDT.s_out})

# Connector constraint 
m.HDT.connectIN = Constraint(expr=m.HDT.IN==m.EV.OUT)

# Constraints
m.HDT.Massenbilanz = Constraint(rule=Massenbilanz)
m.HDT.Arbeit = Constraint(rule=Bruttoarbeit)
m.HDT.Isentrop = Constraint(rule=IsentroperWirkungsgrad)

m.HDT.d1 = Disjunct([0,1], rule=_d_real)
m.HDT.d2 = Disjunct([0,1], rule=_d_isentrop)
m.HDT.c1 = Disjunction(rule=_c1)
m.HDT.c2 = Disjunction(rule=_c2)

'''
Zwischenüberhitzer
- Annahme: isobar
'''

m.ZU = Block()

# Variables 
m.ZU.m_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))
m.ZU.m_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))

m.ZU.p_in = Var(initialize=38.0, within=NonNegativeReals, bounds=(0.03,270.0))
m.ZU.p_out = Var(initialize=38.0, within=NonNegativeReals, bounds=(0.03,270.0))

m.ZU.T_in = Var(initialize=400.0, within=NonNegativeReals, bounds=(293.15,893.15)) 
m.ZU.T_out = Var(initialize=774.0, within=NonNegativeReals, bounds=(293.15,893.15)) 

m.ZU.h_in = Var(initialize=3000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.ZU.h_out = Var(initialize=5000.0, within=NonNegativeReals, bounds=(10.0,10000.0))

m.ZU.s_in = Var(initialize=5.0, within=NonNegativeReals, bounds=(1.0, 11.0))
m.ZU.s_out = Var(initialize=8.0, within=NonNegativeReals, bounds=(1.0, 11.0))

m.ZU.Q = Var(initialize=10000.0, within=NonNegativeReals, bounds=(10000.0, 500000.0))

# Connector 
m.ZU.IN = Connector(initialize=
                        {'Druck':       m.ZU.p_in,
                         'Temperatur':  m.ZU.T_in,
                         'Massenstrom': m.ZU.m_in,
                         'Enthalpie':   m.ZU.h_in,
                         'Entropie':    m.ZU.s_in})
    
m.ZU.OUT = Connector(initialize=
                        {'Druck':       m.ZU.p_out,
                         'Temperatur':  m.ZU.T_out,
                         'Massenstrom': m.ZU.m_out,
                         'Enthalpie':   m.ZU.h_out,
                         'Entropie':    m.ZU.s_out})
 
# Connector constraints 
m.ZU.connectIN = Constraint(expr=m.ZU.IN==m.HDT.OUT)

# Constraint
m.ZU.Massenbilanz = Constraint(rule=Massenbilanz)
m.ZU.Druck = Constraint(rule=Isobar)
m.ZU.Waerme = Constraint(rule=Waermestrom_zu)

# NDT-Eintrittspunkt ist auf jeden Fall im überhitzten Bereich! Also hier die Gleichung für die Enthalpie im überhittzten Bereich 
# (Funktion wurde ganz am Anfang des Codes deklariert)
m.ZU.Austrittsenthalpie = Constraint(rule=Enthalpie_ueberhitzt_vereinfacht)
m.ZU.Austrittsentropie = Constraint(rule=Entropie_ueberhitzt_vereinfacht)

'''
Block für die Niederdruckturbine
'''

m.NDT = Block()

# Variables 

m.NDT.W = Var(initialize=1000.0, within=NonNegativeReals, bounds=(40000.0,1000000.0))

m.NDT.eta = Param(default=0.9)

m.NDT.m_in = Var(initialize=10.0, within=NonNegativeReals, bounds=(1.0,200.0))
m.NDT.m_out = Var(initialize=10.0, within=NonNegativeReals, bounds=(1.0,200.0))

m.NDT.p_in = Var(initialize=38.0, within=NonNegativeReals, bounds=(0.03,270.0)) 
m.NDT.p_out = Var(initialize=0.03, within=NonNegativeReals, bounds=(0.03,270.0)) 
#m.NDT.p_s = Var(initialize=0.03, within=NonNegativeReals, bounds=(0.03,100.0))       # Druck bei einer isentropen Entspannung ist gleich den Druck bei der realen Entspannung

m.NDT.T_in = Var(initialize=774.0, within=NonNegativeReals, bounds=(293.15,893.15)) 
m.NDT.T_out = Var(initialize=600.0, within=NonNegativeReals, bounds=(120.0,893.15)) 
m.NDT.T_s = Var(initialize=400.0)#, within=NonNegativeReals, bounds=(292.0,893.15))       # Austrittstemperatur bei einer isentropen Entspannung

m.NDT.h_in = Var(initialize=5000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.NDT.h_out = Var(initialize=3000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.NDT.h_s = Var(initialize=2500.0, within=NonNegativeReals, bounds=(10.0,10000.0))       # Austrittsenthalpie bei einer isentropen Entspannung


m.NDT.h_Tau = Var(initialize=2900.0)#, within=NonNegativeReals, bounds=(10.0,10000.0))     # Enthalpie gesättigter Dampf
m.NDT.h_Siede = Var(initialize=1000.0)#, within=NonNegativeReals, bounds=(10.0,10000.0))   # Enthalpie siedender Flüssigkeit 
m.NDT.h_Siede_s = Var(initialize=2900.0)#, within=NonNegativeReals, bounds=(10.0,10000.0)) 
m.NDT.h_Tau_s = Var(initialize=1000.0)#, within=NonNegativeReals, bounds=(10.0,10000.0)) 

m.NDT.s_in = Var(initialize=8.0, within=NonNegativeReals, bounds=(0.0, 11.0))
m.NDT.s_out = Var(initialize=10.0, within=NonNegativeReals, bounds=(0.0, 11.0))
m.NDT.s_s = Var(initialize=8.0, within=NonNegativeReals, bounds=(0.0, 11.0))

m.NDT.s_Siede = Var(initialize=7.0)#, within=NonNegativeReals, bounds=(0.0, 11.0))
m.NDT.s_Tau = Var(initialize=1.0)#, within=NonNegativeReals, bounds=(0.0, 11.0))
m.NDT.s_Siede_s = Var(initialize=7.0)#, within=NonNegativeReals, bounds=(0.0, 11.0))
m.NDT.s_Tau_s = Var(initialize=1.0)#, within=NonNegativeReals, bounds=(0.0, 11.0))

m.NDT.x_Dampf = Var(initialize=0.85, within=NonNegativeReals, bounds=(0.85,1.0))
m.NDT.x_Dampf_s = Var(initialize=0.85, within=NonNegativeReals, bounds=(0.85,1.0))   

# Connector 
m.NDT.IN = Connector(initialize=
                        {'Druck':       m.NDT.p_in,
                         'Temperatur':  m.NDT.T_in,
                         'Massenstrom': m.NDT.m_in,
                         'Enthalpie':   m.NDT.h_in,
                         'Entropie':    m.NDT.s_in})
    
m.NDT.OUT = Connector(initialize=
                        {'Druck':       m.NDT.p_out,
                         'Temperatur':  m.NDT.T_out,
                         'Massenstrom': m.NDT.m_out,
                         'Enthalpie':   m.NDT.h_out,
                         'Entropie':    m.NDT.s_out})

# Connector constraint:
m.NDT.connectIN = Constraint(expr=m.NDT.IN==m.ZU.OUT)

# Constraints
m.NDT.Massenbilanz = Constraint(rule=Massenbilanz)
m.NDT.Arbeit = Constraint(rule=Bruttoarbeit)
m.NDT.eta_isentrop = Constraint(rule=IsentroperWirkungsgrad)

m.NDT.d1 = Disjunct([0,1], rule=_d_real)
m.NDT.d2 = Disjunct([0,1], rule=_d_isentrop)
m.NDT.c1 = Disjunction(rule=_c1)
m.NDT.c2 = Disjunction(rule=_c2)

'''
Kondensator
- Annahme: isobar 
'''

m.K = Block()

# Variables 

m.K.m_in = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))
m.K.m_out = Var(initialize=100.0, within=NonNegativeReals, bounds=(1.0,200.0))

m.K.p_in = Var(initialize=0.03, within=NonNegativeReals, bounds=(0.03,270.0)) 
m.K.p_out = Var(initialize=0.03, within=NonNegativeReals, bounds=(0.03,270.0))

m.K.T_in = Var(initialize=500.0, within=NonNegativeReals, bounds=(125.0,893.15)) 
m.K.T_out = Param(default=293.15)

m.K.h_in = Var(initialize=3000.0, within=NonNegativeReals, bounds=(10.0,10000.0))
m.K.h_out = Var(initialize=83.914, within=NonNegativeReals, bounds=(10.0,10000.0))

m.K.s_in = Var(initialize=10.0, within=NonNegativeReals, bounds=(1.0, 11.0))
m.K.s_out = Var(initialize=0.365, within=NonNegativeReals, bounds=(0.0, 11.0))

m.K.Q = Var(initialize=289570.8, within=NonNegativeReals, bounds=(10000.0,300000.0))

# Connector 
m.K.IN = Connector(initialize=
                        {'Druck':       m.K.p_in,
                         'Temperatur':  m.K.T_in,
                         'Massenstrom': m.K.m_in,
                         'Enthalpie':   m.K.h_in,
                         'Entropie':    m.K.s_in})
    
m.K.OUT = Connector(initialize=
                        {'Druck':       m.K.p_out,
                         'Temperatur':  m.K.T_out,
                         'Massenstrom': m.K.m_out,
                         'Enthalpie':   m.K.h_out,
                         'Entropie':    m.K.s_out})
# Connector constraints
m.K.connectIN = Constraint(expr=m.K.IN==m.NDT.OUT)
m.K.connectOUT = Constraint(expr=m.K.OUT==m.P.IN)

# Constraints
m.K.Massenbilanz = Constraint(rule=Massenbilanz)
m.K.Druck = Constraint(rule=Isobar)
m.K.Waerme = Constraint(rule=Waermestrom_ab)


##################################################################

# Back to Parent Block:

# Declaration of general variables and constraints 
# Declaration of the objective function

##################################################################

    
# Applying the connector and the GDP extensions:
TransformationFactory('core.expand_connectors').apply_to(m)
TransformationFactory('gdp.bigm').apply_to(m, bigM=(-10000.0, 10000.0))

m.W_netto = Param(default=100000.0) # Nettoarbeit in kW 


def Nettoarbeit(m):
    return m.W_netto == - m.P.W + (m.HDT.W + m.NDT.W)
def _limit(m):
    return m.EV.Q + m.ZU.Q >= 100000.0

m.Nettoleistung = Constraint(rule=Nettoarbeit)
m.Begrenzung = Constraint(rule=_limit)

m.obj = Objective(expr=m.W_netto / (m.EV.Q + m.ZU.Q), sense=maximize)


##################################################################

# Solving the model:

# Declaration of the solver 
# Solver options

##############################################################

from pyomo.opt import SolverFactory 

solver = SolverFactory('bonmin')

Results = solver.solve(m,tee=True)

# For solving on the Neos-Server
'''
solver_manager = SolverManagerFactory('neos')
Results = solver_manager.solve(m, opt=solver, tee=True) 
'''

assert str(Results.solver.status) == "ok"
assert str(Results.solver.termination_condition) == "optimal" 

# Standard Output of the solver 
with open('bonmin_Results.yml', 'w') as f:
    f.write("Results of optimization!\n")
    m.display(ostream=f)