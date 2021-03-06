# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:58:21 2018

@author: Nesrine Ouanes
"""

from pyomo.environ import *
from iapws import IAPWS97

# Allgemeine gegebene Werte 
T_0 =    298.15 # K
p_0 =    1.013  # bar

LHV  =   50000.0 # Heizwert CH4[kJ/kg]

# Isentropenexponent, ideale Gase:
kappa_a = 1.4    
kappa_g = 1.33

# Molare Massen 
# M = Summe(x_i * M_i)
M_Luft =  28.6384 # kg/kmol // Luftzusammensetzung aus TDO-Buch Seite 30 
M_CH4 =   16.0426   # kg/kmol // Annahme: Erdgas besteht aus reinem Methan (CH4)

'''
# Zusammenseztung im Abgas und Luft: 
# Aus TDO Buch Seite 30.                             
'''
'''
Referenzszustände für Luft und Verbrennungsgas (bei T_0 ubd p_0): 
>>> Berechnung siehe Excel "Refenrenzenthalpien" 
'''
h_Luft_ref = -164.580434    # kJ/kg  

# 
h_Abgas_ref =  -1207.992721 # kJ/kg 

comp = ConcreteModel()

comp.CGAM_W_net = Param(default=30000.0) #kW

'''
# Parameter der Kostenfunktionen
# Werte aus der 10.ET 1 -Übung, Spalte (2) "hohe Komponentenkosten": 
'''

init_c = {(1,1):71.1, (1,2):0.9,                            # AC Parameters
          (2,1):46.08, (2,2):0.995, (2,3):0.018, (2,4):26.4,# CC Parameters 
          (3,1):479.34, (3,2):0.92, (3,3):0.036, (3,4):54.4,# GT Parameters 
          (4,1):4122,                                      # APH Parameters 
          (5,1):6570, (5,2):21276, (5,3):1184.4 }          # HRSG Parameters 
                                         
comp.K = Set(initialize=[(1,1),(1,2),
                          (2,1),(2,2),(2,3),(2,4),
                          (3,1),(3,2),(3,3),(3,4),
                          (4,1),
                          (5,1),(5,2),(5,3)])

def c_init(model,i,j):
    return init_c[(i,j)]

comp.c = Param(comp.K,default=c_init)

# U-Wert für Luftvorwärmer
comp.U = Param(default=0.018)

'''
# Parameter der Stoffwertpolynome, je nach Stoff: (Polynome aus Buch Thermal Design and Optimization Appendix C))
# Der 6. Parameter M enstpricht die molare Masse des Stoffes in kg/kmol -> dient zur Umrechnung von molenspez. auf massenspez. Enthalpie 
'''

comp.W = Set(initialize=['H_plus','a','b','c','d','M'])

# Parameter für N2
init_w_N2 = {'H_plus': -7.069,
             'a':       24.229,
             'b':       10.521,
             'c':       0.180,
             'd':      -2.315,
             'M':       28.0}

def w_N2_init(model, a):
    return init_w_N2[a] 
    
comp.w_N2 = Param(comp.W, default=w_N2_init)   

# Parameter für O2
init_w_O2 = {'H_plus': -9.589,
             'a':       29.154,
             'b':       6.477,
             'c':      -0.184,
             'd':      -1.017,
             'M':       32.0} 

def w_O2_init(model, a):
    return init_w_O2[a] 
    
comp.w_O2 = Param(comp.W, default=w_O2_init)

# Parameter für CH4
init_w_CH4 = {'H_plus': -81.242,
              'a':       11.933,
              'b':       77.647,
              'c':       0.142,
              'd':      -18.414,
              'M':       16.0} 

def w_CH4_init(model, a):
    return init_w_CH4[a] 
    
comp.w_CH4 = Param(comp.W, default=w_CH4_init)

# Parameter für H2O gas
init_w_H2O = {'H_plus': -253.871,
              'a':       34.376,
              'b':       7.841,
              'c':      -0.423,
              'd':       0.0,
              'M':       18.0} 

def w_H2O_init(model, a):
    return init_w_H2O[a] 
    
comp.w_H2O = Param(comp.W, default=w_H2O_init)
    
# Parameter für CO
init_w_CO2 = {'H_plus': -413.886,
              'a':       51.128,
              'b':       4.368,
              'c':      -1.469,
              'd':       0.0,
              'M':       44.0} 

def w_CO2_init(model, a):
    return init_w_CO2[a] 
    
comp.w_CO2 = Param(comp.W, default=w_CO2_init)

'''
Block for AC component
'''
comp.AC = Block()

# Variables in AC Block

comp.AC.eta = Var(initialize=0.8,within=PositiveReals, bounds=(0.5,0.8999))             # Obere Grenze sollte model.c(1,2) sein!
comp.AC.rp = Var(initialize=9.0, within=PositiveReals, bounds=(1.0,16.0))               # bounds from TDO book page 85
comp.AC.W =   Var(initialize=34507.03, within=PositiveReals, bounds=(10000.0,50000.0))  # kW
comp.AC.Z = Var(within=PositiveReals,initialize=300000.0)                               # Kostenfunktion    

# Ein- und Austrittszustände
comp.AC.m_in  = Var(initialize=112.18, within=PositiveReals,bounds=(10,500.0))
comp.AC.m_out  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,500.0))

comp.AC.p_in = Param(default=p_0) 
comp.AC.p_out = Var(initialize=9.117,within=PositiveReals, bounds=(1.013,20.0))

comp.AC.T_in = Param(default=T_0) 
comp.AC.T_out = Var(initialize=604.52, within=PositiveReals, bounds=(298.15,800.0))

# Enthalpien in kJ/kg 
comp.AC.h_in = Param(default=h_Luft_ref) 
comp.AC.h_out = Var(initialize=1000.0, within=PositiveReals, bounds=(1.0,1000.0))

# Rules for AC constraints
def massenbilanz_AC(model):
    return model.m_in == model.m_out 
def mech_arbeit(model):
    return model.W == model.m_in * (model.h_out - model.h_in)

def austrittstemperatur(model):
    return model.T_out == model.T_in * (1 + 1/model.eta * ((model.p_out/model.p_in)**((kappa_a-1)/kappa_a) - 1)) 

def Austrittsdruck_AC(model):
    return model.rp == model.p_out*model.p_in**(-1)

def Austrittsenthalpie_AC(model, w, a, b, c):
    return model.h_out == ((w['H_plus']*10**3 + w['a']*model.T_out + w['b']/2*(model.T_out)**2 / 10**3 - w['c']*(model.T_out)**(-1) * 10**6 + w['d']/3 * (model.T_out/10**2)**3)/w['M']*0.7748 +
                           (a['H_plus']*10**3 + a['a']*model.T_out + a['b']/2*(model.T_out)**2 / 10**3 - a['c']*(model.T_out)**(-1) * 10**6 + a['d']/3 * (model.T_out/10**2)**3)/a['M']*0.2059 +
                           (b['H_plus']*10**3 + b['a']*model.T_out + b['b']/2*(model.T_out)**2 / 10**3 - b['c']*(model.T_out)**(-1) * 10**6 + b['d']/3 * (model.T_out/10**2)**3)/b['M']*0.0003 +
                           (c['H_plus']*10**3 + c['a']*model.T_out + c['b']/2*(model.T_out)**2 / 10**3 - c['c']*(model.T_out)**(-1) * 10**6 + c['d']/3 * (model.T_out/10**2)**3)/c['M']*0.019)
# Constraints
comp.AC.massenb = Constraint(rule=massenbilanz_AC) 
comp.AC.Enthalpie_out = Constraint(expr=Austrittsenthalpie_AC(comp.AC, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O)) 
comp.AC.arbeit = Constraint(rule=mech_arbeit)
comp.AC.austrittstemperatur = Constraint(rule=austrittstemperatur)    
comp.AC.austrittsdruck = Constraint(rule=Austrittsdruck_AC) # Druckverhältnis 

comp.AC.Kosten = Constraint(expr= comp.AC.Z == comp.c[1,1] * comp.AC.m_in*(comp.c[1,2] - comp.AC.eta)**(-1) *
                            comp.AC.p_out * comp.AC.p_in**(-1) * log(comp.AC.p_out/comp.AC.p_in))

# Connector 
comp.AC.OUT = Connector(initialize=
                        {'Druck':       comp.AC.p_out,
                         'Temperatur':  comp.AC.T_out,
                         'Massenstrom': comp.AC.m_out,
                         'Enthalpie':   comp.AC.h_out})
'''
Block for APH component
'''
comp.APH = Block()

# Variables 
comp.APH.deltaP_a = Param(default=0.05)
comp.APH.deltaP_g = Param(default=0.03)
comp.APH.deltaT = Var(initialize=39.741,within=PositiveReals) # logarithmische Temperaturdifferenz
comp.APH.Z = Var(within=PositiveReals,initialize=2000.0)  # Kostenfunktion

# Ein- und Austrittszustände
#    Kalte Seite (Luft)
comp.APH.m_in_c  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,500.0))
comp.APH.m_out_c  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,500.0))

comp.APH.p_in_c = Var(initialize=9.117,within=PositiveReals, bounds=(1.013,20.0))
comp.APH.p_out_c = Var(initialize=8.661,within=PositiveReals, bounds=(1.013,20.0))

comp.APH.T_in_c = Var(initialize=604.52, within=PositiveReals, bounds=(298.15,800.0))
comp.APH.T_out_c = Var(initialize=894.02, within=PositiveReals, bounds=(298.15,1200.0))

comp.APH.h_in_c = Var(initialize=1000.0)#, within=PositiveReals, bounds=(1.0,2000.0)) # kJ/kg
comp.APH.h_out_c = Var(initialize=1500.0)#, within=PositiveReals, bounds=(1.0,2500.0))

#    Connector für kalte Seite 
comp.APH.IN_cold = Connector(initialize=
                        {'Druck':       comp.APH.p_in_c,
                         'Temperatur':  comp.APH.T_in_c,
                         'Massenstrom': comp.APH.m_in_c,
                         'Enthalpie':   comp.APH.h_in_c})
comp.APH.OUT_cold = Connector(initialize=
                        {'Druck':       comp.APH.p_out_c,
                         'Temperatur':  comp.APH.T_out_c,
                         'Massenstrom': comp.APH.m_out_c,
                         'Enthalpie':   comp.APH.h_out_c})
#    Warme Seite (Turbinenabgas)
comp.APH.m_in_w  = Var(initialize=113.68, within=PositiveReals,bounds=(1.0,500.0))
comp.APH.m_out_w  = Var(initialize=113.68, within=PositiveReals,bounds=(1.0,500.0))

comp.APH.p_in_w = Var(initialize=1.099,within=PositiveReals, bounds=(1.013,20.0))
comp.APH.p_out_w = Var(initialize=1.066,within=PositiveReals, bounds=(1.013,20.0))

comp.APH.T_in_w = Var(initialize=915.54, within=PositiveReals, bounds=(298.15,1700.0))
comp.APH.T_out_w = Var(initialize=670.66, within=PositiveReals, bounds=(298.15,1700.0))

comp.APH.h_in_w = Var(initialize=1000.0) # kJ/kg
comp.APH.h_out_w = Var(initialize=800.0)

#    Connector für die warme Seite 
comp.APH.IN_warm = Connector(initialize=
                        {'Druck':       comp.APH.p_in_w,
                         'Temperatur':  comp.APH.T_in_w,
                         'Massenstrom': comp.APH.m_in_w,
                         'Enthalpie':   comp.APH.h_in_w})
comp.APH.OUT_warm = Connector(initialize=
                        {'Druck':       comp.APH.p_out_w,
                         'Temperatur':  comp.APH.T_out_w,
                         'Massenstrom': comp.APH.m_out_w,
                         'Enthalpie':   comp.APH.h_out_w})

# Rules for the constraints
def massenbilanz_APH_kalt(model):
    return model.m_in_c == model.m_out_c 
def massenbilanz_APH_warm(model):
    return model.m_in_w == model.m_out_w
def waermestrom_APH(model):
    return model.m_in_c * (model.h_out_c - model.h_in_c) == model.m_in_w * (model.h_in_w - model.h_out_w)
def austrittsdruck_kalt(model):
    return model.p_out_c == model.p_in_c * (1 - model.deltaP_a)
def austrittsdruck_warm(model):
    return model.p_out_w == model.p_in_w * (1 - model.deltaP_g)
def deltaT_log(model):
    return model.deltaT == (model.T_in_w - model.T_out_c - model.T_out_w + model.T_in_c)*log((model.T_in_w - model.T_out_c)/(model.T_out_w - model.T_in_c))**(-1)
def Austrittsenthalpie_kalt_APH(model, w, a, b, c):
    return model.h_out_c == ((w['H_plus']*10**3 + w['a']*model.T_out_c + w['b']/2*(model.T_out_c)**2 / 10**3 - w['c']*(model.T_out_c)**(-1) * 10**6 + w['d']/3 * (model.T_out_c/10**2)**3)/w['M']*0.7748 +
                           (a['H_plus']*10**3 + a['a']*model.T_out_c + a['b']/2*(model.T_out_c)**2 / 10**3 - a['c']*(model.T_out_c)**(-1) * 10**6 + a['d']/3 * (model.T_out_c/10**2)**3)/a['M']*0.2059 +
                           (b['H_plus']*10**3 + b['a']*model.T_out_c + b['b']/2*(model.T_out_c)**2 / 10**3 - b['c']*(model.T_out_c)**(-1) * 10**6 + b['d']/3 * (model.T_out_c/10**2)**3)/b['M']*0.0003 +
                           (c['H_plus']*10**3 + c['a']*model.T_out_c + c['b']/2*(model.T_out_c)**2 / 10**3 - c['c']*(model.T_out_c)**(-1) * 10**6 + c['d']/3 * (model.T_out_c/10**2)**3)/c['M']*0.019)
def Austrittsenthalpie_warm_APH(model, w, a, b, c):
    return model.h_out_w == ((w['H_plus']*10**3 + w['a']*model.T_out_w + w['b']/2*model.T_out_w**2 / 10**3 - w['c'] * (model.T_out_w)**(-1) * 10**6 + w['d']/3*(model.T_out_w/10**2)**3)/w['M']*0.7507  + 
                             (a['H_plus']*10**3 + a['a']*model.T_out_w + a['b']/2*model.T_out_w**2 / 10**3 - a['c'] * (model.T_out_w)**(-1) * 10**6 + a['d']/3*(model.T_out_w/10**2)**3)/a['M']*0.1372 + 
                             (b['H_plus']*10**3 + b['a']*model.T_out_w + b['b']/2*model.T_out_w**2 / 10**3 - b['c'] * (model.T_out_w)**(-1) * 10**6 + b['d']/3*(model.T_out_w/10**2)**3 )/b['M']*0.0314 +
                             (c['H_plus']*10**3 + c['a']*model.T_out_w + c['b']/2*model.T_out_w**2 / 10**3 - c['c'] * (model.T_out_w)**(-1) * 10**6 + c['d']/3*(model.T_out_w/10**2)**3 )/c['M']*0.0807)

# Constraints for APH

comp.APH.massenb_kalt = Constraint(rule=massenbilanz_APH_kalt)
comp.APH.massenb_warm = Constraint(rule=massenbilanz_APH_warm)
 
comp.APH.Enthalpie_out_c = Constraint(expr=Austrittsenthalpie_kalt_APH(comp.APH, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O)) 
comp.APH.Enthalpie_out_w = Constraint(expr=Austrittsenthalpie_warm_APH(comp.APH,comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O)) 

comp.APH.waerme = Constraint(rule=waermestrom_APH)
comp.APH.austrittP_kalt = Constraint(rule=austrittsdruck_kalt)
comp.APH.austrittP_warm = Constraint(rule=austrittsdruck_warm)
comp.APH.deltaT_log = Constraint(rule=deltaT_log)


comp.APH.Kosten = Constraint(expr= comp.APH.Z == comp.c[4,1]*(comp.APH.m_in_w * (comp.APH.h_in_w - 
                             comp.APH.h_out_w) / (comp.U * comp.APH.deltaT)) ** 0.6)

# Connector constraint: 
comp.APH_connectIN_cold = Constraint(expr=comp.APH.IN_cold==comp.AC.OUT)

'''
Block for CC Component
'''
# Variables
comp.CC = Block()
comp.CC.deltaP = Param(default=0.05)
comp.CC.eta = Param(default=0.98)
comp.CC.Z = Var(within=PositiveReals,initialize=100000.0) # Kostenfunktion 

# Ein- und Austrittszustände 
comp.CC.m_in = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,500.0))
comp.CC.m_out = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,500.0))
comp.CC.m_fuel = Var(initialize=1.73, within=PositiveReals,bounds=(1.0,20.0))

comp.CC.T_in = Var(initialize=894.02, within=PositiveReals, bounds=(298.15,1700.0))
comp.CC.T_out = Var(initialize=1400.0, within=PositiveReals, bounds=(298.15,1700.0))
comp.CC.T_fuel = Param(default=T_0)

comp.CC.p_in = Var(initialize=8.661,within=PositiveReals, bounds=(1.013,20.0))
comp.CC.p_out = Var(initialize=8.228,within=PositiveReals, bounds=(1.013,20.0))
# p_fuel soll gleich p_in sein!  

comp.CC.h_in = Var(initialize=1500.0) # kJ/kg
comp.CC.h_out = Var(initialize=2200.0)
#comp.CC.h_fuel = Param(default=909.841) # 909.841196424 in kJ/kg:  Wert mit Coolprop berechnet

# Connectors 
comp.CC.IN = Connector(initialize=
                        {'Druck':       comp.CC.p_in,
                         'Temperatur':  comp.CC.T_in,
                         'Massenstrom': comp.CC.m_in,
                         'Enthalpie':   comp.CC.h_in})
    
comp.CC.OUT = Connector(initialize=
                        {'Druck':       comp.CC.p_out,
                         'Temperatur':  comp.CC.T_out,
                         'Massenstrom': comp.CC.m_out,
                         'Enthalpie':   comp.CC.h_out})

# Rules for the constraints
def massenbilanz_CC(model):
    return model.m_out == model.m_in + model.m_fuel
def waermestrom_CC(model):
    return model.m_fuel * model.eta * LHV == model.m_out * (model.h_out - h_Abgas_ref)  - model.m_in * (model.h_in - h_Luft_ref)
def austrittsdruck(model):
    return model.p_out == model.p_in * (1-model.deltaP)
def Austrittsenthalpie_CC(model, w, a, b, c):
    return model.h_out == ((w['H_plus']*10**3 + w['a']*model.T_out + w['b']/2*model.T_out**2 / 10**3 - w['c'] * (model.T_out)**(-1) * 10**6 + w['d']/3*(model.T_out/10**2)**3)/w['M']*0.7507  + 
                             (a['H_plus']*10**3 + a['a']*model.T_out + a['b']/2*model.T_out**2 / 10**3 - a['c'] * (model.T_out)**(-1) * 10**6 + a['d']/3*(model.T_out/10**2)**3)/a['M']*0.1372 + 
                             (b['H_plus']*10**3 + b['a']*model.T_out + b['b']/2*model.T_out**2 / 10**3 - b['c'] * (model.T_out)**(-1) * 10**6 + b['d']/3*(model.T_out/10**2)**3 )/b['M']*0.0314 +
                             (c['H_plus']*10**3 + c['a']*model.T_out + c['b']/2*model.T_out**2 / 10**3 - c['c'] * (model.T_out)**(-1) * 10**6 + c['d']/3*(model.T_out/10**2)**3 )/c['M']*0.0807)
# Constraints 
comp.CC.massenb = Constraint(rule= massenbilanz_CC)
comp.CC.Enthalpie_out = Constraint(expr=Austrittsenthalpie_CC(comp.CC, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O))

comp.CC.waerme = Constraint(rule= waermestrom_CC)
comp.CC.austritt = Constraint(rule= austrittsdruck)

comp.CC.Kosten = Constraint(expr= comp.CC.Z == comp.c[2,1] * comp.CC.m_in / (comp.c[2,2] - (1 - comp.CC.deltaP)) *  
                                                (1 + exp(comp.c[2,3] * comp.CC.T_out - comp.c[2,4])))

# Connector constraint
comp.CC_connectIN = Constraint(expr= comp.CC.IN == comp.APH.OUT_cold)

'''
Block for GT component
'''
comp.GT = Block()

# Variables 
comp.GT.W = Var(initialize=64433.09, within=PositiveReals, bounds=(30000.0,80000.0))
comp.GT.eta = Var(initialize=0.88, within=PositiveReals, bounds=(0.5,0.919999))  # Beachte obere Grenze sollte model.c32 sein.
comp.GT.Z = Var(within=PositiveReals,initialize=100000.0)                        # Kostenfunktion

# Ein- und Austrittszustände
comp.GT.m_in  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,500.0))
comp.GT.m_out  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,500.0))

comp.GT.p_in = Var(initialize=8.228,within=PositiveReals, bounds=(1.013,20.0))
comp.GT.p_out = Var(initialize=1.099,within=PositiveReals, bounds=(1.013,20.0))

comp.GT.T_in = Var(initialize=1400.0, within=PositiveReals, bounds=(298.15,1700.0))
comp.GT.T_out = Var(initialize=915.54, within=PositiveReals, bounds=(298.15,1700.0))

comp.GT.h_in = Var(initialize=1200.0) # kJ/kg
comp.GT.h_out = Var(initialize=1000.0)

# Connector
comp.GT.IN = Connector(initialize=
                        {'Druck':       comp.GT.p_in,
                         'Temperatur':  comp.GT.T_in,
                         'Massenstrom': comp.GT.m_in,
                         'Enthalpie':   comp.GT.h_in})
    
comp.GT.OUT = Connector(initialize=
                        {'Druck':       comp.GT.p_out,
                         'Temperatur':  comp.GT.T_out,
                         'Massenstrom': comp.GT.m_out,
                         'Enthalpie':   comp.GT.h_out})
# Rules for constraints
def bruttoarbeit(model):
    return model.W == model.m_in * (model.h_in - model.h_out)
def austritttemp(model):
    return model.T_out == model.T_in*(1-model.eta*(1-(model.p_in/model.p_out)**((1-kappa_g)/kappa_g)))
def Austrittsenthalpie_GT(model, w, a, b, c):
    return model.h_out == ((w['H_plus']*10**3 + w['a']*model.T_out + w['b']/2*model.T_out**2 / 10**3 - w['c'] * (model.T_out)**(-1) * 10**6 + w['d']/3*(model.T_out/10**2)**3)/w['M']*0.7507  + 
                           (a['H_plus']*10**3 + a['a']*model.T_out + a['b']/2*model.T_out**2 / 10**3 - a['c'] * (model.T_out)**(-1) * 10**6 + a['d']/3*(model.T_out/10**2)**3)/a['M']*0.1372 + 
                           (b['H_plus']*10**3 + b['a']*model.T_out + b['b']/2*model.T_out**2 / 10**3 - b['c'] * (model.T_out)**(-1) * 10**6 + b['d']/3*(model.T_out/10**2)**3 )/b['M']*0.0314 +
                           (c['H_plus']*10**3 + c['a']*model.T_out + c['b']/2*model.T_out**2 / 10**3 - c['c'] * (model.T_out)**(-1) * 10**6 + c['d']/3*(model.T_out/10**2)**3 )/c['M']*0.0807)
# Constraints
#   Massenbilanz vom AC verwendet
comp.GT.massenb = Constraint(rule=massenbilanz_AC) 
comp.GT.Enthalpie_out = Constraint(expr=Austrittsenthalpie_GT(comp.GT, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O))
comp.GT.austritttemp = Constraint(rule=austritttemp)
comp.GT.bruttoarbeit = Constraint(rule=bruttoarbeit)

comp.GT.Kosten = Constraint(expr= comp.GT.Z == comp.c[3,1] * comp.GT.m_in * (comp.c[3,2] - comp.GT.eta)**(-1) *
                            log(comp.GT.p_in / comp.GT.p_out) * (1 + exp(comp.c[3,3] * comp.GT.T_in - comp.c[3,4])))

# Connector constraints:
comp.GT_connectIN = Constraint(expr=comp.GT.IN==comp.CC.OUT)
comp.GT_connectOUT = Constraint(expr=comp.GT.OUT==comp.APH.IN_warm)

'''
Block for HRSG component
'''

comp.HRSG = Block()
# Variables 
comp.HRSG.deltaP = Param(default=0.05)
comp.HRSG.deltaT_A = Param(default=15.0) # gegebene Approach-Temperaturdifferenz
comp.HRSG.deltaT_P = Var(initialize=5.0, within=PositiveReals, bounds=(1.0,60.0))
comp.HRSG.Z = Var(within=PositiveReals,initialize=2285850.112) # Kostenfunktion

# PH
comp.HRSG.PH_Q = Var(initialize=15696.1496, within=PositiveReals)
comp.HRSG.PH_deltaT = Var(initialize=58.0489, within=PositiveReals) # logarithmische Temeraturdifferenz

# EV 
comp.HRSG.EV_Q = Var(initialize=121987.6194)
comp.HRSG.EV_deltaT = Var(initialize=90.1297, within=PositiveReals) # logarithmische Temperaturdifferenz

# Ein- und Austrittszustände
#    Kalte Seite (Wasserdampf)
comp.HRSG.m_dampf = Param(default=14.0)
comp.HRSG.p_dampf = Param(default=20.0) # Annahme: keine Druckverluste auf der Wasserseite des Abhitzekessels
#    Dampftemperaturen
T_Sattdampf = IAPWS97(P=value(comp.HRSG.p_dampf)*0.1,x=1).T
comp.HRSG.T_in_c = Param(default=T_0)
comp.HRSG.T_mitte_c = Param(default=T_Sattdampf - comp.HRSG.deltaT_A)  
comp.HRSG.T_out_c = Param(default=T_Sattdampf)

comp.HRSG.h_in_c = Param(default=106.468) # Wert mit Coolprop berechnet: 106.468191318 
comp.HRSG.h_mitte_c = Var(initialize=100.0, within=NonNegativeReals, bounds=(0.0,3000.0)) 
comp.HRSG.h_out_c = Param(default=2798.293) # Wert mit Coolprop berechnet: 2798.2926032 
  
#    Connector für kalte Seite 
comp.HRSG.IN_cold = Connector(initialize=
                              {'Druck':       comp.HRSG.p_dampf,
                               'Temperatur':  comp.HRSG.T_in_c,
                               'Massenstrom': comp.HRSG.m_dampf,
                               'Enthalpie':   comp.HRSG.h_in_c})
comp.HRSG.OUT_cold = Connector(initialize=
                        {'Druck':       comp.HRSG.p_dampf,
                         'Temperatur':  comp.HRSG.T_out_c,
                         'Massenstrom': comp.HRSG.m_dampf,
                         'Enthalpie':   comp.HRSG.h_out_c})
#    Warme Seite (Turbinenabgas)
comp.HRSG.m_in_w  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,500.0))
comp.HRSG.m_out_w  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,500.0))

comp.HRSG.p_in_w = Var(initialize=1.066,within=PositiveReals, bounds=(1.013,20.0))
comp.HRSG.p_mitte_w = Var(initialize=1.039,within=PositiveReals, bounds=(1.013,20.0))
comp.HRSG.p_out_w = Param(default=p_0)

comp.HRSG.T_in_w = Var(initialize=670.66, within=PositiveReals, bounds=(298.15,900.0))
comp.HRSG.T_mitte_w = Var(initialize=505.53, within=PositiveReals, bounds=(298.15,800.0))
comp.HRSG.T_out_w = Var(initialize=387.65, within=PositiveReals, bounds=(298.15,600.0))

comp.HRSG.h_in_w = Var(initialize=400.0) #kJ/kg
comp.HRSG.h_mitte_w = Var(initialize=300.0)
comp.HRSG.h_out_w = Var(initialize=200.0)

#   Connector für die warme Seite 
comp.HRSG.IN_warm = Connector(initialize=
                        {'Druck':       comp.HRSG.p_in_w,
                         'Temperatur':  comp.HRSG.T_in_w,
                         'Massenstrom': comp.HRSG.m_in_w,
                         'Enthalpie':   comp.HRSG.h_in_w})
comp.HRSG.OUT_warm = Connector(initialize=
                        {'Druck':       comp.HRSG.p_out_w,
                         'Temperatur':  comp.HRSG.T_out_w,
                         'Massenstrom': comp.HRSG.m_out_w,
                         'Enthalpie':   comp.HRSG.h_out_w})

# Rules for HRSG constraints
def massenbilanz_HRSG(model):
    return model.m_in_w == model.m_out_w
def pinch_deltaT(model):
    return model.T_mitte_w == model.T_out_c + model.deltaT_P
def waermestrom_HRSG(model):
    return (model.m_in_w * (model.h_in_w-model.h_out_w)==model.m_dampf*
            (IAPWS97(P=value(model.p_dampf)*0.1,x=1).h-IAPWS97(P=value(model.p_dampf)*0.1,
                    T=value(model.T_in_c)).h))
def waermestrom_EV(model):
    return model.EV_Q == model.m_dampf*(IAPWS97(P=value(model.p_dampf)*0.1,x=1).h-
                                        IAPWS97(P=value(model.p_dampf)*0.1,
                                               T=value(model.T_mitte_c)).h) 
def druckverluste_PH(model):
    return model.p_out_w == model.p_mitte_w * (1-model.deltaP/2)
def waermestrom_PH(model):
    return model.PH_Q == model.m_in_w * (model.h_mitte_w - model.h_out_w) 
def waerme_EV(model):
    return model.EV_Q == model.m_in_w * (model.h_in_w - model.h_mitte_w)

# Rechenregeln für die logarithmischen Temperaturen
def logT_PH(model):
    return model.PH_deltaT == ((model.T_mitte_w - model.T_mitte_c - model.T_out_w + model.T_in_c)*log((model.T_mitte_w - model.T_mitte_c)/(model.T_out_w - model.T_in_c))**(-1))

def logT_EV(model):
    return model.EV_deltaT == ((model.T_mitte_w - model.T_mitte_c - model.T_in_w + model.T_out_c)*log((model.T_mitte_w - model.T_mitte_c)/(model.T_in_w - model.T_out_c))**(-1))

def druckverluste_EV(model):
    return model.p_mitte_w == model.p_in_w*(1-model.deltaP/2)
    
def Austrittsenthalpie_warm_HRSG(model, w, a, b, c):
    return model.h_out_w == ((w['H_plus']*10**3 + w['a']*model.T_out_w + w['b']/2*model.T_out_w**2 / 10**3 - w['c'] * (model.T_out_w)**(-1) * 10**6 + w['d']/3*(model.T_out_w/10**2)**3)/w['M']*0.7507  + 
                           (a['H_plus']*10**3 + a['a']*model.T_out_w + a['b']/2*model.T_out_w**2 / 10**3 - a['c'] * (model.T_out_w)**(-1) * 10**6 + a['d']/3*(model.T_out_w/10**2)**3)/a['M']*0.1372 + 
                           (b['H_plus']*10**3 + b['a']*model.T_out_w + b['b']/2*model.T_out_w**2 / 10**3 - b['c'] * (model.T_out_w)**(-1) * 10**6 + b['d']/3*(model.T_out_w/10**2)**3 )/b['M']*0.0314 +
                           (c['H_plus']*10**3 + c['a']*model.T_out_w + c['b']/2*model.T_out_w**2 / 10**3 - c['c'] * (model.T_out_w)**(-1) * 10**6 + c['d']/3*(model.T_out_w/10**2)**3 )/c['M']*0.0807)

def Austrittsenthalpie_mitte_HRSG(model, w, a, b, c):
    return model.h_mitte_w == ((w['H_plus']*10**3 + w['a']*model.T_mitte_w + w['b']/2*model.T_mitte_w**2 / 10**3 - w['c'] * (model.T_mitte_w)**(-1) * 10**6 + w['d']/3*(model.T_mitte_w/10**2)**3)/w['M']*0.7507  + 
                               (a['H_plus']*10**3 + a['a']*model.T_mitte_w + a['b']/2*model.T_mitte_w**2 / 10**3 - a['c'] * (model.T_mitte_w)**(-1) * 10**6 + a['d']/3*(model.T_mitte_w/10**2)**3)/a['M']*0.1372 + 
                               (b['H_plus']*10**3 + b['a']*model.T_mitte_w + b['b']/2*model.T_mitte_w**2 / 10**3 - b['c'] * (model.T_mitte_w)**(-1) * 10**6 + b['d']/3*(model.T_mitte_w/10**2)**3 )/b['M']*0.0314 +
                               (c['H_plus']*10**3 + c['a']*model.T_mitte_w + c['b']/2*model.T_mitte_w**2 / 10**3 - c['c'] * (model.T_mitte_w)**(-1) * 10**6 + c['d']/3*(model.T_mitte_w/10**2)**3 )/c['M']*0.0807)
# HRSG constraints
comp.HRSG.massenbilanz_warm = Constraint(rule=massenbilanz_HRSG)
comp.HRSG.PinchT = Constraint(rule=pinch_deltaT)


comp.HRSG.Enthalpie_out_w = Constraint(expr=Austrittsenthalpie_warm_HRSG(comp.HRSG, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O)) 
comp.HRSG.Enthalpie_mitte_w = Constraint(expr=Austrittsenthalpie_mitte_HRSG(comp.HRSG, comp.w_N2, comp.w_O2, comp.w_CO2, comp.w_H2O)) 

comp.HRSG.waermestrom = Constraint(rule=waermestrom_HRSG) # Enthalpien für Wasser/Dampf mit IAPWS97 
comp.HRSG.waermestrom_EV = Constraint(rule=waermestrom_EV)
comp.HRSG.druckverluste_PH = Constraint(rule=druckverluste_PH)
comp.HRSG.PH_waermestrom = Constraint(rule=waermestrom_PH)
comp.HRSG.EV_waermestrom = Constraint(rule=waerme_EV) 
comp.HRSG.EV_logT = Constraint(rule=logT_EV)
comp.HRSG.PH_logT = Constraint(rule=logT_PH)
comp.HRSG.druckverluste_EV = Constraint(rule=druckverluste_EV)

comp.HRSG.Kosten = Constraint(expr= comp.HRSG.Z == comp.c[5,1] * ((comp.HRSG.PH_Q / comp.HRSG.PH_deltaT) ** 0.8 +
                              (comp.HRSG.EV_Q / comp.HRSG.EV_deltaT) ** 0.8) + comp.c[5,2] * comp.HRSG.m_dampf + comp.c[5,3] *
                               (comp.HRSG.m_in_w) ** 1.2)

# Connector constraint
comp.HRSG.connectIN = Constraint(expr=comp.HRSG.IN_warm==comp.APH.OUT_warm)


'''      
Applying the connector extention
'''
TransformationFactory('core.expand_connectors').apply_to(comp)

'''
Back to Parent block
'''
def nettoarbeit(model):
    return model.GT.W - model.AC.W == model.CGAM_W_net 

comp.Nettoarbeit = Constraint(rule=nettoarbeit)

# Objective Declaration
comp.obj = Objective(expr=comp.CGAM_W_net/(comp.CC.m_fuel * LHV), sense=maximize)

# Optimization: Solver definition
#from pyomo.opt import SolverFactory

opt = SolverFactory("ipopt")

# Increasing the maximum number of iterations for Ipopt or Bonmin 
#opt.options['max_iter'] = 10000 
#opt.options['linear_solver'] = 'ma27'

Results = opt.solve(comp, tee=True) #, keepfiles=True) 


'''

# Calling Neos Server 
#from pyomo.opt import SolverManagerFactory

solver_manager = SolverManagerFactory('neos') 
Results = solver_manager.solve(comp, opt=opt, tee=True)  #keepfiles=True
'''



assert str(Results.solver.status) == "ok"
assert str(Results.solver.termination_condition) == "optimal" 


'''
Print results to a file 
'''

# Standard Output of the solver 
with open('ipopt_Thermodyn_Results.yml', 'w') as f:
    f.write("Results of optimization!\n")
    comp.display(ostream=f)

#######################  PrettyDIsplay function  #############################
   
#Formatierte Ausgabe
def display(f):
    
    CGAM_eta=comp.CGAM_W_net.value/(comp.CC.m_fuel.value*LHV)
    print (" ")    
    # Formatierte Ausgabe
    print(" #  m[kg/s]    p[bar]    h[kJ/kg]    T[K]    T[°C]", file=f)
    print("---------------------------------------------------", file=f)
    print("Verdichter (AC): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(1,comp.AC.m_in.value,
          comp.AC.p_in.value, comp.AC.T_in.value, comp.AC.T_in.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(2,comp.AC.m_out.value,
          comp.AC.p_out.value, comp.AC.T_out.value, comp.AC.T_out.value-273.15),'\n', file=f)
    
    print("Luftvorwärmer (APH): ",'\n', file=f)
    print("    Kalte Seite (Luft): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(2,comp.APH.m_in_c.value,
          comp.APH.p_in_c.value, comp.APH.T_in_c.value, comp.APH.T_in_c.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(3,comp.APH.m_out_c.value,
          comp.APH.p_out_c.value, comp.APH.T_out_c.value, comp.APH.T_out_c.value-273.15), '\n', file=f)
    print("    Warme Seite (Abgas): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(5,comp.APH.m_in_w.value,
          comp.APH.p_in_w.value, comp.APH.T_in_w.value, comp.APH.T_in_w.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(6,comp.APH.m_out_w.value,
          comp.APH.p_out_w.value, comp.APH.T_out_w.value, comp.APH.T_out_w.value-273.15),'\n', file=f)
    
    print("Brennkammer (CC): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(3,comp.CC.m_in.value,
          comp.CC.p_in.value, comp.CC.T_in.value, comp.CC.T_in.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(4,comp.CC.m_out.value,
          comp.CC.p_out.value, comp.CC.T_out.value, comp.CC.T_out.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(12,comp.CC.m_fuel.value,
          comp.CC.p_in.value, comp.CC.T_fuel.value, comp.CC.T_fuel.value-273.15),'\n', file=f)
    
    print("Gasturbine (GT): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(4,comp.GT.m_in.value,
          comp.GT.p_in.value, comp.GT.T_in.value, comp.GT.T_in.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(5,comp.GT.m_out.value,
          comp.GT.p_out.value, comp.GT.T_out.value, comp.GT.T_out.value-273.15),'\n', file=f)
    
    print("Abhiltzekessel (HRSG): ", '\n', file=f)
    print("    Kalte Seite (Wasserdampf): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(9,comp.HRSG.m_dampf.value,
          comp.HRSG.p_dampf.value, comp.HRSG.T_in_c.value, comp.HRSG.T_in_c.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(10,comp.HRSG.m_dampf.value,
          comp.HRSG.p_dampf.value, comp.HRSG.T_mitte_c.value, comp.HRSG.T_mitte_c.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(11,comp.HRSG.m_dampf.value,
          comp.HRSG.p_dampf.value, comp.HRSG.T_out_c.value, comp.HRSG.T_out_c.value-273.15),'\n', file=f)
    print("    Warme Seite (Abgas): ", '\n', file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(6,comp.HRSG.m_in_w.value,
          comp.HRSG.p_in_w.value, comp.HRSG.T_in_w.value, comp.HRSG.T_in_w.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(7,comp.HRSG.m_in_w.value,
          comp.HRSG.p_mitte_w.value, comp.HRSG.T_mitte_w.value, comp.HRSG.T_mitte_w.value-273.15), file=f)
    print("{:2.0f}  {:7.2f}  {:8.3f} {:7.2f}  {:7.2f}".format(8,comp.HRSG.m_out_w.value,
          comp.HRSG.p_out_w.value, comp.HRSG.T_out_w.value, comp.HRSG.T_out_w.value-273.15), '\n', file=f)

    print (" ", file=f)
    print("W_AC: {:7.2f}".format(comp.AC.W.value), 'kW', file=f)
    print("W_GT: {:7.2f}".format(comp.GT.W.value), 'kW', file=f)
    print("eta: {:6.3f}".format(CGAM_eta*100), '%', file=f)
    print(" ", file=f)
    print("## Entscheidungsvariablen:",'\n', file=f)
    print("     Druckverhältnis im Verdichter: {:3.3f}".format(comp.AC.rp.value), file=f)
    print("     Isentroper Wirkungsgrad des Verdichters: {:6.3f} % und der Turbine {:6.3f} %".format(comp.AC.eta.value*100, 
          comp.GT.eta.value*100), file=f)
    print("     Temperaturdifferenz am Pinch: {:3.3f}".format(comp.HRSG.deltaT_P.value), 'K', file=f)
    print("     Temperatur am Turbineneintritt: {:3.3f}".format(comp.GT.T_in.value), 'K', '\n', file=f)



# Formatted Output with the display function defined above
with open('ipopt_PrettyResults.yml', 'w') as f:
    f.write("Results of optimization!\n")
    print('  ', file=f)
    display(f)