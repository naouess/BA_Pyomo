# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:58:21 2018

@author: Nesrine Ouanes 
"""

from pyomo.environ import *
from iapws import IAPWS97

# Allgemeine gegebene Werte 
T_0 =    298.15 #K
p_0 =    1.013  #bar

LHV  =   50000.0 # Heizwert CH4[kJ/kg]

# Spez. Wärmekapazität [kJ/kgK]:
cp_a =   1.004   
cp_g =   1.17   
# Isentropenexponent, ideale Gase:
kappa_a = 1.4    
kappa_g = 1.33

comp = ConcreteModel()

comp.CGAM_W_net = Param(default=30000.0) # kW

# Parameter der Kostenfunktionen
# Werte aus der Übung, Spalte (1): (CGAM Paper Valero 1994)

init_c = {(1,1):71.1, (1,2):0.9,                            # AC Parameter
          (2,1):46.08, (2,2):0.995, (2,3):0.018, (2,4):26.4,# CC Parameter
          (3,1):479.34, (3,2):0.92, (3,3):0.036, (3,4):54.4,# GT Parameter
          (4,1):4122,                                      # APH Parameter
          (5,1):6570, (5,2):21276, (5,3):1184.4 }          # HRSG Parameter
                                        
comp.K = Set(initialize=[(1,1),(1,2),
                          (2,1),(2,2),(2,3),(2,4),
                          (3,1),(3,2),(3,3),(3,4),
                          (4,1),
                          (5,1),(5,2),(5,3)])

def c_init(model,i,j):
    return init_c[(i,j)]

comp.c = Param(comp.K,default=c_init)

# U-Wert für Luftvorwärmer:
comp.U = Param(default=0.018)

'''
Block for AC component
'''
comp.AC = Block()

# Variablen im AC-Block

comp.AC.eta = Var(initialize=0.85,within=PositiveReals, bounds=(0.5,0.899999))           # Obere Grenze sollte model.c(1,2) sein!
comp.AC.rp =  Var(initialize=9.0, within=PositiveReals, bounds=(1.0,60.0))
comp.AC.W =   Var(initialize=34507.03, within=PositiveReals, bounds=(10000.0,100000.0))  # kW
comp.AC.Z = Var(within=PositiveReals,initialize=3154564.6245)                            # Kostenfunktion   

# Ein- und Austrittszustände
comp.AC.m_in  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,1000.0))
comp.AC.m_out  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,1000.0))


comp.AC.p_in = Param(default=p_0) 
comp.AC.p_out = Var(initialize=9.117,within=PositiveReals, bounds=(1.013,50.0))

comp.AC.T_in = Param(default=T_0) 
comp.AC.T_out = Var(initialize=604.52, within=PositiveReals, bounds=(298.15,1700.0))


# Rules for AC-constraints
def massenbilanz_AC(model):
    return model.m_in == model.m_out     
def mech_arbeit(model):
    return model.W == model.m_in * cp_a * (model.T_out - model.T_in)
def austrittstemperatur(model):
    return model.T_out == model.T_in * (1 + 1/model.eta * ((model.p_out/model.p_in)**((kappa_a-1)/kappa_a) - 1)) 
def Austrittsdruck_AC(model):
    return model.rp == model.p_out*model.p_in**(-1)

# Constraints
comp.AC.massenb = Constraint(rule=massenbilanz_AC)  
comp.AC.arbeit = Constraint(rule=mech_arbeit)
comp.AC.austritt = Constraint(rule=austrittstemperatur)    
comp.AC.druckverhaeltnis = Constraint(rule=Austrittsdruck_AC)

comp.AC.Kosten = Constraint(expr= comp.AC.Z == (comp.c[1,1] * comp.AC.m_in*(comp.c[1,2] - comp.AC.eta)**(-1)) *
                           (comp.AC.p_out * comp.AC.p_in**(-1)) * log(comp.AC.p_out/comp.AC.p_in))

# Connector 
comp.AC.OUT = Connector(initialize=
                        {'Druck':       comp.AC.p_out,
                         'Temperatur':  comp.AC.T_out,
                         'Massenstrom': comp.AC.m_out})
'''
Block for APH component
'''
comp.APH = Block()

# Variables 
comp.APH.deltaP_a = Param(default=0.05)
comp.APH.deltaP_g = Param(default=0.03)
comp.APH.deltaT = Var(initialize=39.741,within=PositiveReals) # logarithmische Temperaturdifferenz [K]
comp.APH.Z = Var(within=PositiveReals,initialize=2572.673)    # Kostenfunktion

# Ein- und Austrittszustände
#     Kalte Seite (Luft)
comp.APH.m_in_c  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,1000.0))
comp.APH.m_out_c  = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,1000.0))

comp.APH.p_in_c = Var(initialize=9.117,within=PositiveReals, bounds=(1.013,50.0))
comp.APH.p_out_c = Var(initialize=8.661,within=PositiveReals, bounds=(1.013,50.0))

comp.APH.T_in_c = Var(initialize=604.52, within=PositiveReals, bounds=(298.15,1700.0))
comp.APH.T_out_c = Var(initialize=894.02, within=PositiveReals, bounds=(298.15,1700.0))

#     Connector für kalte Seite 
comp.APH.IN_cold = Connector(initialize=
                        {'Druck':       comp.APH.p_in_c,
                         'Temperatur':  comp.APH.T_in_c,
                         'Massenstrom': comp.APH.m_in_c})
comp.APH.OUT_cold = Connector(initialize=
                        {'Druck':       comp.APH.p_out_c,
                         'Temperatur':  comp.APH.T_out_c,
                         'Massenstrom': comp.APH.m_out_c})

#     Warme Seite (Turbinenabgas)
comp.APH.m_in_w  = Var(initialize=113.68, within=PositiveReals,bounds=(1.0,1000.0))
comp.APH.m_out_w  = Var(initialize=113.68, within=PositiveReals,bounds=(1.0,1000.0))

comp.APH.p_in_w = Var(initialize=1.099,within=PositiveReals, bounds=(1.013,50.0))
comp.APH.p_out_w = Var(initialize=1.066,within=PositiveReals, bounds=(1.013,50.0))

comp.APH.T_in_w = Var(initialize=915.54, within=PositiveReals, bounds=(298.15,1700.0))
comp.APH.T_out_w = Var(initialize=670.66, within=PositiveReals, bounds=(298.15,1700.0))

#     Connector für die warme Seite 
comp.APH.IN_warm = Connector(initialize=
                        {'Druck':       comp.APH.p_in_w,
                         'Temperatur':  comp.APH.T_in_w,
                         'Massenstrom': comp.APH.m_in_w})
comp.APH.OUT_warm = Connector(initialize=
                        {'Druck':       comp.APH.p_out_w,
                         'Temperatur':  comp.APH.T_out_w,
                         'Massenstrom': comp.APH.m_out_w})

# Rules for the constraints
def massenbilanz_APH_kalt(model):
    return model.m_in_c == model.m_out_c 
def massenbilanz_APH_warm(model):
    return model.m_in_w == model.m_out_w
def waermestrom_APH(model):
    return model.m_in_c * cp_a * (model.T_out_c - model.T_in_c) == model.m_in_w * cp_g * (model.T_in_w - model.T_out_w)
def austrittsdruck_kalt(model):
    return model.p_out_c == model.p_in_c * (1 - model.deltaP_a)
def austrittsdruck_warm(model):
    return model.p_out_w == model.p_in_w * (1 - model.deltaP_g)
def deltaT_log(model):
    return model.deltaT == (model.T_in_w - model.T_out_c - model.T_out_w + model.T_in_c)*(log((model.T_in_w - model.T_out_c)/(model.T_out_w - model.T_in_c)))**(-1)

# Constraints for APH

comp.APH.massenb_kalt = Constraint(rule=massenbilanz_APH_kalt)
comp.APH.massenb_warm = Constraint(rule=massenbilanz_APH_warm)
comp.APH.waerme = Constraint(rule=waermestrom_APH)
comp.APH.austrittP_kalt = Constraint(rule=austrittsdruck_kalt)
comp.APH.austrittP_warm = Constraint(rule=austrittsdruck_warm)
comp.APH.deltaT_log = Constraint(rule=deltaT_log)

comp.APH.Kosten = Constraint(expr= comp.APH.Z == comp.c[4,1]*(comp.APH.m_in_w * cp_g * (comp.APH.T_in_w - 
                             comp.APH.T_out_w) / (comp.U * comp.APH.deltaT)) ** 0.6)

# Connector constraint: 
comp.APH_connectIN_cold = Constraint(expr=comp.APH.IN_cold==comp.AC.OUT)

'''
Block for CC Component
'''
# Variablen
comp.CC = Block()
comp.CC.deltaP = Param(default=0.05)                        # Druckverlust
comp.CC.eta = Param(default=0.98)                           # Wirkungsgrad
comp.CC.Z = Var(within=PositiveReals,initialize=149473.758) # Kostenfunktion 

# Ein- und Austrittszustände 
comp.CC.m_in = Var(initialize=112.18, within=PositiveReals,bounds=(1.0,1000.0))
comp.CC.m_out = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,1000.0))
comp.CC.m_fuel = Var(initialize=1.73, within=PositiveReals,bounds=(1.0,50.0))

comp.CC.T_in = Var(initialize=894.02, within=PositiveReals, bounds=(298.15,1700.0))
comp.CC.T_out = Var(initialize=1400.0, within=PositiveReals, bounds=(298.15,1700.0))
comp.CC.T_fuel = Param(default=T_0)

comp.CC.p_in = Var(initialize=8.661,within=PositiveReals, bounds=(1.013,50.0))
comp.CC.p_out = Var(initialize=8.228,within=PositiveReals, bounds=(1.013,50.0))

# Connectors 
comp.CC.IN = Connector(initialize=
                        {'Druck':       comp.CC.p_in,
                         'Temperatur':  comp.CC.T_in,
                         'Massenstrom': comp.CC.m_in})
    
comp.CC.OUT = Connector(initialize=
                        {'Druck':       comp.CC.p_out,
                         'Temperatur':  comp.CC.T_out,
                         'Massenstrom': comp.CC.m_out})

# Rules for the constraints
def massenbilanz_CC(model):
    return model.m_out == model.m_in + model.m_fuel
def waermestrom_CC(model):
    return model.m_fuel * model.eta * LHV == model.m_out * cp_g * (model.T_out - T_0) - model.m_in * cp_a * (model.T_in - T_0)
def austrittsdruck(model):
    return model.p_out == model.p_in * (1-model.deltaP)

# Constraints 
comp.CC.massenb = Constraint(rule= massenbilanz_CC)
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

# Variablen 
comp.GT.W = Var(initialize=64433.09, within=PositiveReals, bounds=(10000.0,100000.0))
comp.GT.eta = Var(initialize=0.88, within=PositiveReals, bounds=(0.5,0.919999))    # Obere Grenze sollte model.c32 sein.
comp.GT.Z = Var(within=PositiveReals,initialize=2796365.220)                       # Kostenfunktion

# Ein- und Austrittszustände
comp.GT.m_in  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,1000.0))
comp.GT.m_out  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,1000.0))

comp.GT.p_in = Var(initialize=8.228, within=PositiveReals, bounds=(1.013,50.0))
comp.GT.p_out = Var(initialize=1.099,within=PositiveReals, bounds=(1.013,50.0))

comp.GT.T_in = Var(initialize=1400.0, within=PositiveReals, bounds=(298.15,1700.0))
comp.GT.T_out = Var(initialize=915.54, within=PositiveReals, bounds=(298.15,1700.0))

# Connector
comp.GT.IN = Connector(initialize=
                        {'Druck':       comp.GT.p_in,
                         'Temperatur':  comp.GT.T_in,
                         'Massenstrom': comp.GT.m_in})
    
comp.GT.OUT = Connector(initialize=
                        {'Druck':       comp.GT.p_out,
                         'Temperatur':  comp.GT.T_out,
                         'Massenstrom': comp.GT.m_out})
# Rules for constraints
def bruttoarbeit(model):
    return model.W == model.m_in * cp_g * (model.T_in - model.T_out)
def austritttemp(model):
    return model.T_out == model.T_in*(1-model.eta*(1-(model.p_in/model.p_out)**((1-kappa_g)/kappa_g)))

# Constraints

comp.GT.massenb = Constraint(rule=massenbilanz_AC)
comp.GT.bruttoarbeit = Constraint(rule=bruttoarbeit)
comp.GT.austritttemp = Constraint(rule=austritttemp)

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
comp.HRSG.deltaT_A = Param(default=15.0) # Gegebene Approach-Temperaturdifferenz
comp.HRSG.deltaT_P =  within=PositiveReals)ar(initialize=20.0, within=PositiveReals, bounds=(1.0,50.0))
comp.HRSG.Z = Var(within=PositiveReals,initialize=2285850.112) #Kostenfunktion

# PH
comp.HRSG.PH_Q = Var(initialize=15696.1496, within=PositiveReals)
comp.HRSG.PH_deltaT = Var(initialize=58.0489, within=PositiveReals) # Logarithmische Temeraturdifferenz

# EV 
comp.HRSG.EV_Q = Var(initialize=121987.6194, within=PositiveReals)
comp.HRSG.EV_deltaT = Var(initialize=90.1297, within=PositiveReals) # Logarithmische Temperaturdifferenz

# Ein- und Austrittszustände
#     Kalte Seite (Wasserdampf)
comp.HRSG.m_dampf = Param(default=14.0)
comp.HRSG.p_dampf = Param(default=20.0) # Annahme: keine Druckverluste auf der Wasserseite des Abhitzekessels
#     Dampftemperaturen
T_Sattdampf = IAPWS97(P=value(comp.HRSG.p_dampf)*0.1,x=1).T
comp.HRSG.T_in_c = Param(default=T_0)
comp.HRSG.T_mitte_c = Param(default=T_Sattdampf - comp.HRSG.deltaT_A)  
comp.HRSG.T_out_c = Param(default=T_Sattdampf)    
#     Connector für kalte Seite 
comp.HRSG.IN_cold = Connector(initialize=
                              {'Druck':       comp.HRSG.p_dampf,
                               'Temperatur':  comp.HRSG.T_in_c,
                               'Massenstrom': comp.HRSG.m_dampf})
comp.HRSG.OUT_cold = Connector(initialize=
                        {'Druck':       comp.HRSG.p_dampf,
                         'Temperatur':  comp.HRSG.T_out_c,
                         'Massenstrom': comp.HRSG.m_dampf})
#     Warme Seite (Turbinenabgas)
comp.HRSG.m_in_w  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,1000.0))
comp.HRSG.m_out_w  = Var(initialize=113.81, within=PositiveReals,bounds=(1.0,1000.0))

comp.HRSG.p_in_w = Var(initialize=1.066,within=PositiveReals, bounds=(1.013,50.0))
comp.HRSG.p_mitte_w = Var(initialize=1.039,within=PositiveReals, bounds=(1.013,50.0))
comp.HRSG.p_out_w = Param(default=p_0)

comp.HRSG.T_in_w = Var(initialize=670.66, within=PositiveReals, bounds=(298.15,1700.0))
comp.HRSG.T_mitte_w = Var(initialize=505.53, within=PositiveReals, bounds=(298.15,1700.0))
comp.HRSG.T_out_w = Var(initialize=387.65, within=PositiveReals, bounds=(400.0,1700.0))

#     Connector für die warme Seite 
comp.HRSG.IN_warm = Connector(initialize=
                        {'Druck':       comp.HRSG.p_in_w,
                         'Temperatur':  comp.HRSG.T_in_w,
                         'Massenstrom': comp.HRSG.m_in_w})
comp.HRSG.OUT_warm = Connector(initialize=
                        {'Druck':       comp.HRSG.p_out_w,
                         'Temperatur':  comp.HRSG.T_out_w,
                         'Massenstrom': comp.HRSG.m_out_w})

# Rules for HRSG constraints
def massenbilanz_HRSG(model):
    return model.m_in_w == model.m_out_w
def pinch_deltaT(model):
    return model.T_mitte_w == model.deltaT_P + model.T_out_c
def waermestrom_HRSG(model):
    return (model.m_in_w*cp_g*(model.T_in_w-model.T_out_w)==model.m_dampf*
            (IAPWS97(P=value(model.p_dampf)*0.1,x=1).h-IAPWS97(P=value(model.p_dampf)*0.1,
                    T=value(model.T_in_c)).h))
def waermestrom_EV(model):
    return model.EV_Q == model.m_dampf*(IAPWS97(P=value(model.p_dampf)*0.1,x=1).h-
                                        IAPWS97(P=value(model.p_dampf)*0.1,
                                               T=value(model.T_mitte_c)).h)  
    
def druckverluste_PH(model):
    return model.p_out_w == model.p_mitte_w*(1-model.deltaP/2)
def Eintrittsdruck_EV(model):
    return model.p_mitte_w == model.p_in_w*(1-model.deltaP/2)
def waermestrom_PH(model):
    return model.PH_Q==model.m_in_w*cp_g*(model.T_mitte_w-model.T_out_w) 
def waerme_EV(model):
    return model.EV_Q==model.m_in_w*cp_g*(model.T_in_w-model.T_mitte_w)

# Rechenregeln für die logarithmischen Temperaturen
def logT_PH(model):
    return model.PH_deltaT == ((model.T_mitte_w - model.T_mitte_c - model.T_out_w + model.T_in_c)*(log((model.T_mitte_w - model.T_mitte_c)/(model.T_out_w - model.T_in_c)))**(-1))

def logT_EV(model):
    return model.EV_deltaT == ((model.T_mitte_w - model.T_mitte_c - model.T_in_w + model.T_out_c)*(log((model.T_mitte_w - model.T_mitte_c)/(model.T_in_w - model.T_out_c)))**(-1))

# HRSG constraints
comp.HRSG.massenbilanz_warm = Constraint(rule=massenbilanz_HRSG)
comp.HRSG.PinchT_1 = Constraint(rule=pinch_deltaT)
comp.HRSG.waermestrom = Constraint(rule=waermestrom_HRSG)  # Verbrennungsgas als ideales Gas, Enthalpien für Wasser/Dampf mit IAPWS97 
comp.HRSG.waermestrom_EV = Constraint(rule=waermestrom_EV) # Constraint declared under the HRSG block instead of EV, since it has variables that are only available in the HRSG block
comp.HRSG.druckverluste_PH = Constraint(rule=druckverluste_PH)
comp.HRSG.PH_waermestrom = Constraint(rule=waermestrom_PH)
comp.HRSG.EV_waermestrom = Constraint(rule=waerme_EV) 
comp.HRSG.EV_logT = Constraint(rule=logT_EV)
comp.HRSG.PH_logT = Constraint(rule=logT_PH)
comp.HRSG.Eintrittsdruck_EV = Constraint(rule=Eintrittsdruck_EV)

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

# Parameters for the objective function 
comp.CRF = Param(default=0.182)         # CRF
comp.phi = Param(default=1.06)          # Betriebs- und Wartungsfaktor (maintenance factor)
comp.N = Param(default=8000)            # Jahresvolllastbenutzung in Stunden
comp.C_f = Param(default=0.000004)      # Brennstoffkosten in dollar/kJ mit der Annahme Steigerungsrate = 0
comp.SummZ = Var(within=PositiveReals, initialize=8388826.3868)  # Summe der nivellierten kapitalgebundene Kosten und Betrieb- und Wartungskosten

def _SummeKosten(model):
    return model.SummZ == model.AC.Z + model.CC.Z + model.APH.Z + model.GT.Z + model.HRSG.Z

comp.SummeKosten = Constraint(rule=_SummeKosten)

# Objective Declaration
comp.obj = Objective(expr= comp.C_f * LHV * comp.CC.m_fuel + comp.SummZ * comp.CRF * comp.phi/ 
                      (comp.N * 3600.0),sense=minimize)

# Optimization: Solver definition
from pyomo.opt import SolverFactory
opt = SolverFactory("ipopt")
# opt.options['maxit'] = 80000 #-> increasing the maximum number of iterations for ipopt


# The following formulation is used, if local solvers are to be used. 

Results = opt.solve(comp, tee=True) #keepfiles=True)
'''

# Following Formulation is used, if the Neos server is to be used.  
# Calling Neos Server 

solver_manager = SolverManagerFactory('neos') 
Results = solver_manager.solve(comp, opt=opt, tee=True) #keepfiles=True

'''
# Making sure that the server status is okay and the solution is optimal, before
# loading the results to the model

assert str(Results.solver.status) == "ok"
assert str(Results.solver.termination_condition) == "optimal"    

comp.solutions.store_to(Results)

#Solution display function
def display(f):
    
    CGAM_eta=comp.CGAM_W_net.value/(comp.CC.m_fuel.value*LHV)
    print (" ")    
    #Formatierte Ausgabe
    print(" #  m[kg/s]    p[bar]    T[K]    T[°C]", file=f)
    print("--------------------------------------", file=f)
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
    print("Brennstoffkosten: {:3.4f}".format(comp.C_f.value*LHV*comp.CC.m_fuel.value), '$/s', file=f)
    print("Kapitalgebundene Kosten und Betrieb und Wartungskosten: {:3.4f}".format(comp.SummZ.value * 
          comp.CRF.value * comp.phi.value/(comp.N.value * 3600)), '$/s', file=f)
    print("Gesamtkosten: {:3.4f}".format( comp.C_f.value*LHV*comp.CC.m_fuel.value + comp.SummZ.value * 
          comp.CRF.value * comp.phi.value/(comp.N.value * 3600)), '$/s','\n', file=f)

'''
Print results to a file 
'''

# Standard Output of the solver 
with open('ipopt_Results.yml', 'w') as f:
    f.write("Results of optimization!\n")
    comp.display(ostream=f)

# Formatted Output with the display function defined above
with open('ipopt_PrettyResults.yml', 'w') as f:
    f.write("Results of optimization!\n")
    print('  ', file=f)
    display(f)

