# -*- coding: utf-8 -*-
"""
Created on Thur July 14 09:42:04 2022

@author: zhangx
"""

from casadi import*
import pandas as pd
import numpy as np
import math

# %% import data

df = pd.read_excel('Final_result.xlsx', header = 0)
df = np.array(df)

z_min = 0
z_max = 1

gap = 2

# %% steps
def pressurization(n_space, yE_feed, T_feed, P_high):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_feed = SX.sym('PA_feed', 1)                  # feed state
    PE_feed = SX.sym('PE_feed', 1)                  # feed state
    E_input = SX.sym('E_input', 1)                  # energy input

    # left boundary conditions
    AE_L1 = yE[0] - yE_feed
    AE_L2 = P[0] - P_high
    AE_L3 = T[0] - T_feed
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_half[n_space]
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - 0
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE_feed
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i]-0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4)) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 1,...,n_space
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance, index 1,...,n_space
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance, n_space equations, index 1,...,n_space
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1-Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1-Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i]-Tw[i]) / R_inner) / EB[i])

    # wall energy balance, n_space equations, index 1, ..., n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 1,...,n_space
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space equations, index 1,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i])

    # mass balances
    ODE_PA_feed = P[0] * u[0] * yA[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_feed = P[0] * u[0] * yE[0] * CSA * Porosity_bed / R / T[0]
    ODE_E_input = P[0] * u[0]
    
    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_feed,     PE_feed,     E_input)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_feed, ODE_PE_feed, ODE_E_input)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]

def adsorption(n_space, yE_feed, T_feed, P_high, P_abs_end):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K 
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_feed = SX.sym('PA_feed', 1)                  # feed state
    PE_feed = SX.sym('PE_feed', 1)                  # feed state
    PA_out  = SX.sym('PA_out',  1)                  # output state
    PE_out  = SX.sym('PE_out',  1)                  # output state
    E_input = SX.sym('E_input', 1)                  # energy for input
    E_out   = SX.sym('E_out',   1)                  # energy for output

    # left boundary conditions
    AE_L1 = yE[0] - yE_feed
    AE_L2 = P[0] - P_high
    AE_L3 = T[0] - T_feed
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_abs_end
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - u_half[n_space]
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE_feed
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i] - 0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4) ) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4) )

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4) )
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 1,...,n_space
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])

    # gas and adsorbent energy balance, n_space equations, index 0,...,n_space-1
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g - Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g - Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance, n_space equations, index 1, ..., n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 1,...,n_space
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space equations, index 1,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balances
    ODE_PA_feed = P[0] * u[0] * yA[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_feed = P[0] * u[0] * yE[0] * CSA * Porosity_bed / R / T[0]
    ODE_PA_out = P[n_space+1] * u[n_space+1] * yA[n_space+1] * CSA * Porosity_bed / R / T[0]
    ODE_PE_out = P[n_space+1] * u[n_space+1] * yE[n_space+1] * CSA * Porosity_bed / R / T[0]
    ODE_E_input = P[0] * u[0]
    ODE_E_out = P[n_space+1] * u[n_space+1]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],     qE[1:n_space+1],     PA_feed,     PE_feed,     PA_out,     PE_out,     E_input,     E_out)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1], ODE_qE[1:n_space+1], ODE_PA_feed, ODE_PE_feed, ODE_PA_out, ODE_PE_out, ODE_E_input, ODE_E_out)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]

def equalization_one(n_space, P_equ_end):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_out  = SX.sym('PA_out',  1)                  # output state
    PE_out  = SX.sym('PE_out',  1)                  # output state
    E_out   = SX.sym('E_out',   1)                  # energy for output

    # left boundary conditions
    AE_L1 = yE[0] - yE[1]
    AE_L2 = P[0] - P[1]
    AE_L3 = T[0] - T[1]
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_equ_end
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - u_half[n_space]
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE[1]
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i]-0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4) ) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P[1]
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4) )

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4) )
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance,        n_space equations, index 1,...,n_space
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance,        n_space equation, index 1,...,n_space
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance,     n_space equations, index 1,...,n_space
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance,                  n_space equations, index 1,...,n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance,                   n_space equations, index 1,...,n_space
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation,                       n_space equations, index 1,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balances
    ODE_PA_out = u[n_space+1] * yA[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PE_out = u[n_space+1] * yE[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_E_out = u[n_space+1] * P[n_space+1]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_out,       PE_out,     E_out)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_out,   ODE_PE_out, ODE_E_out)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]

def rinse(n_space, yE_rinse, T_feed, P_high, P_rin_end):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_feed = SX.sym('PA_feed', 1)                  # feed state
    PE_feed = SX.sym('PE_feed', 1)                  # feed state
    PA_out  = SX.sym('PA_out',  1)                  # output state
    PE_out  = SX.sym('PE_out',  1)                  # output state
    E_input = SX.sym('E_input', 1)                  # energy for input
    E_out   = SX.sym('E_out',   1)                  # energy for output

    # left boundary conditions
    AE_L1 = yE[0] - yE_rinse
    AE_L2 = P[0] - P_high
    AE_L3 = T[0] - T_feed
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_rin_end
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - u_half[n_space]
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE_rinse
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i]-0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4) ) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4) )

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4) )
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space+1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 0,...,n_space-1
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance, n_space equations, index 0,...,n_space-1
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance,                  n_space equations, index 1,...,n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 0,...,n_space-1
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space+1 equations, index 0,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balances
    ODE_PA_feed = P[0] * u[0] * yA[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_feed = P[0] * u[0] * yE[0] * CSA * Porosity_bed / R / T[0]
    ODE_PA_out = u[n_space+1] * yA[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PE_out = u[n_space+1] * yE[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_E_input = P[0] * u[0]
    ODE_E_out = P[n_space+1] * u[n_space+1]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_feed,     PE_feed,     PA_out,     PE_out,     E_input,     E_out)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_feed, ODE_PE_feed, ODE_PA_out, ODE_PE_out, ODE_E_input, ODE_E_out)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]    

def desorption(n_space, P_low):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_out  = SX.sym('PA_out',  1)                  # output state
    PE_out  = SX.sym('PE_out',  1)                  # output state
    E_out   = SX.sym('E_out',   1)                  # energy for output

    # left boundary conditions
    AE_L1 = yE[0] - yE_half[0]
    AE_L2 = P[0] - P_low
    AE_L3 = T[0] - T_half[0]
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_half[n_space]
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - u_half[n_space]
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions        
    AE_yE_half = yE_half[0] - yE[1]
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i] + 0.5 * yE[i+1]) * 0.6667 / (Sigma + yE[i] - yE[i+1])**4 + (1.5 * yE[i+1] - 0.5 * yE[i+2]) * 0.3333 / (Sigma + yE[i+1] - yE[i+2])**4) / (0.6667 / (Sigma + yE[i] - yE[i+1])**4 + 0.3333 / (Sigma + yE[i+1] - yE[i+2])**4))
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_low
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i] + 0.5 * P[i+1]) * 0.6667 / (Sigma + P[i] - P[i+1])**4 + (1.5 * P[i+1] - 0.5 * P[i+2]) * 0.3333 / (Sigma + P[i+1] - P[i+2])**4) / (0.6667 / (Sigma + P[i] - P[i+1])**4 + 0.3333 / (Sigma + P[i+1] - P[i+2])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i] + 0.5 * T[i+1]) * 0.6667 / (Sigma + T[i] - T[i+1])**4 + (1.5 * T[i+1] - 0.5 * T[i+2]) * 0.3333 / (Sigma + T[i+1] - T[i+2])**4) / (0.6667 / (Sigma + T[i] - T[i+1])**4 + 0.3333 / (Sigma + T[i+1] - T[i+2])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space+1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 0,...,n_space-1
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance, n_space equations, index 0,...,n_space-1
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance,                  n_space equations, index 1,...,n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 0,...,n_space-1
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space+1 equations, index 0,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balance
    ODE_PA_out = - u[0] * P[0] * yA[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_out = - u[0] * P[0] * yE[0] * CSA * Porosity_bed / R / T[0]
    ODE_E_out = - u[0] * P[0]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_out,     PE_out,     E_out)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_out, ODE_PE_out, ODE_E_out)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]        

def equalization_two(n_space, yE_equ_out, T_equ_out, P_equ_end):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_feed = SX.sym('PA_feed', 1)                  # feed state
    PE_feed = SX.sym('PE_feed', 1)                  # feed state

    # left boundary conditions
    AE_L1 = yE[0] - yE_equ_out
    AE_L2 = P[0] - P_equ_end
    AE_L3 = T[0] - T_equ_out
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_half[n_space]
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - 0
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE_equ_out
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i]-0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4) ) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_equ_end
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4) )
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T_equ_out
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4) )
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 0,...,n_space-1
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance, n_space equations, index 0,...,n_space-1
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance,                  n_space equations, index 1,...,n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 0,...,n_space-1
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space+1 equations, index 0,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balances
    ODE_PA_feed = u[0] * (1 - yE[0]) * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_feed = u[0] * yE[0] * P[0] * CSA * Porosity_bed / R / T[0]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_feed,     PE_feed)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_feed, ODE_PE_feed)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0],  yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1],  yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6,  AE_R1,          AE_R2,          AE_R3,          AE_R6,          AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]

def reflux(n_space, yE_rin_out, T_rin_out, P_rin_end):

    # variables
    yE      = SX.sym('yE',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # state,             unit m/s
    yE_half = SX.sym('yE_half', n_space+1)          # wall state variable
    T_half  = SX.sym('T_half',  n_space+1)          # wall state variable
    P_half  = SX.sym('P_half',  n_space+1)          # wall state variable
    u_half  = SX.sym('u_half',  n_space+1)          # wall state variable
    Tw      = SX.sym('Tw',      n_space+2)          # state,             unit K
    qE      = SX.sym('qE',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qE_Star = SX.sym('qE_Star', n_space+2)          # state,             unit mol/kg 
    qA_Star = SX.sym('qA_Star', n_space+2)          # state,             unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # time differential
    EB      = SX.sym('EB',      n_space+2)          # state
    PA_feed = SX.sym('PA_feed', 1)                  # feed state
    PE_feed = SX.sym('PE_feed', 1)                  # feed state

    # left boundary conditions
    AE_L1 = yE[0] - yE_rin_out
    AE_L2 = P[0] - P_rin_end
    AE_L3 = T[0] - T_rin_out
    AE_L4 = u[0] - u_half[0]
    AE_L5 = EB[0] - 0
    AE_L6 = Tw[0] - T_amb

    # right boundary conditions
    AE_R1 = yE[n_space+1] - yE_half[n_space]
    AE_R2 = P[n_space+1] - P_half[n_space]
    AE_R3 = T[n_space+1] - T_half[n_space]
    AE_R4 = u[n_space+1] - 0
    AE_R5 = EB[n_space+1] - 0
    AE_R6 = Tw[n_space+1] - T_amb

    # WENO scheme for walls + boundary conditions
    AE_yE_half = yE_half[0] - yE_rin_out
    for i in range(1, n_space):
        AE_yE_half = vertcat(AE_yE_half, yE_half[i] - ((0.5 * yE[i+1] + 0.5 * yE[i]) * 0.6667 / (Sigma + yE[i+1] - yE[i])**4 + (1.5 * yE[i]-0.5 * yE[i-1]) * 0.3333 / (Sigma + yE[i] - yE[i-1])**4) / (0.6667 / (Sigma + yE[i+1] - yE[i])**4 + 0.3333 / (Sigma + yE[i] - yE[i-1])**4) ) 
    AE_yE_half = vertcat(AE_yE_half, yE_half[n_space] - yE[n_space])

    AE_P_half = P_half[0] - P_rin_end
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4) )
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T_rin_out
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4) )
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yE_half[0]) + MW_PE * yE_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yE_half[i]) + MW_PE * yE_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yE + yA - 1

    # DSL isotherm & LDF equation
    AE_qA_Star = qA_Star - Qsat_1_PA * exp(k_1_PA + k_2_PA / T) * P * yA / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PA * exp(k_3_PA + k_4_PA / T) * P * yA / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    AE_qE_Star = qE_Star - Qsat_1_PE * exp(k_1_PE + k_2_PE / T) * P * yE / (1 + exp(k_1_PA + k_2_PA / T) * P * yA + exp(k_1_PE + k_2_PE / T) * P * yE) - Qsat_2_PE * exp(k_3_PE + k_4_PE / T) * P * yE / (1 + exp(k_3_PA + k_4_PA / T) * P * yA + exp(k_3_PE + k_4_PE / T) * P * yE)
    ODE_qA = kA * (qA_Star - qA)
    ODE_qE = kE * (qE_Star - qE)

    # Single component mass balance, n_space equations, index 0,...,n_space-1
    ODE_yE = DxZS * ((yE[2] + yE[0] - 2 * yE[1]) - (yE_half[1] - yE_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yE_half[1] - yE_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yE_half[1] - yE_half[0]) / delta_z - M * T[1] * (ODE_qE[1] - yE[1] * (ODE_qA[1] + ODE_qE[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yE = vertcat(ODE_yE, DxZS * ((yE[i+1] + yE[i-1] - 2 * yE[i]) - (yE_half[i] - yE_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yE_half[i] - yE_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yE_half[i] - yE_half[i-1]) / delta_z - M * T[i] * (ODE_qE[i] - yE[i] * (ODE_qA[i] + ODE_qE[i])) / P[i])

    # denominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qE[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qE[i])) - Cp_g * Porosity_bed * P[i] / R / T[i] )

    # gas and adsorbent energy balance, n_space equations, index 0,...,n_space-1
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qE[1]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[1] + delta_H_E * ODE_qE[1]) - 2 * h_W_in * (T[1] - Tw[1]) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1-Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qE[i]) + (1 - Porosity_bed) * (delta_H_A * ODE_qA[i] + delta_H_E * ODE_qE[i]) - 2 * h_W_in * (T[i] - Tw[i]) / R_inner) / EB[i])

    # wall energy balance,                  n_space equations, index 1,...,n_space
    ODE_Tw = EB_1 * (Tw[2] + Tw[0] - 2 * Tw[1]) / delta_zS + EB_2 * (T[1] - Tw[1]) - EB_3 * (Tw[1] - T_amb)
    for i in range(2, n_space + 1):
        ODE_Tw = vertcat(ODE_Tw, EB_1 * (Tw[i+1] + Tw[i-1] - 2 * Tw[i]) / delta_zS + EB_2 * (T[i] - Tw[i]) - EB_3 * (Tw[i] - T_amb)) 

    # Total mass balance, n_space equations, index 0,...,n_space-1
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qE[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qE[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])   

    # Ergun equation, n_space+1 equations, index 0,...,n_space 
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PE * yE[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PE * yE[i]) / R / T[i] )

    # mass balances
    ODE_PA_feed = u[0] * (1 - yE[0]) * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_PE_feed = u[0] * yE[0] * P[0] * CSA * Porosity_bed / R / T[0]

    x = vertcat(    yE[1:n_space+1], P[1:n_space+1], T[1:n_space+1], Tw[1:n_space+1], qA[1:n_space+1],       qE[1:n_space+1],     PA_feed,     PE_feed)
    ode = vertcat(  ODE_yE,          ODE_P,          ODE_T,          ODE_Tw,          ODE_qA[1:n_space+1],   ODE_qE[1:n_space+1], ODE_PA_feed, ODE_PE_feed)

    z = vertcat(    yE[0],  P[0],   T[0],   Tw[0], yE[n_space+1],  P[n_space+1],   T[n_space+1],   Tw[n_space+1], yE_half,    P_half,     T_half,     u_half,     u,                  yA,     qA_Star,    qE_Star,    EB)
    alg = vertcat(  AE_L1,  AE_L2,  AE_L3,  AE_L6, AE_R1,          AE_R2,          AE_R3,          AE_R6,         AE_yE_half, AE_P_half,  AE_T_half,  AE_u_half,  AE_L4,AE_u,AE_R4,   AE_yA,  AE_qA_Star, AE_qE_Star, AE_L5,AE_EB,AE_R5)

    return [x, ode, z, alg]

# %% bed and gas parameters

n_space = 30                        # number of space interval
L = 10                              # bed length: m
delta_z = L / n_space
delta_zS = delta_z * delta_z           
Diameter = L / 5                    # column inner diameter: m
R_inner = Diameter / 2              # Column inner radius: m
R_out = Diameter / 2 + 0.06         # Column outer radius: m

yE_feed = 0.85                      # feed gas composition
yA_feed = 1 - yE_feed
yE_rinse = 0.99                     # input composition of rinse
T_feed = 323                        # feed gas temperature: K
T_amb = 298                         # ambient temperature: K
P_atm = 101325                      # Pa

K_z = 0.0903                        # thermal diffusivity of gas: J/m/s/K
Dx = 5e-5                           # axial dispersion: m2/s from adsorption kinetics of propane and propylene in zeolite 4A
DxZS = Dx / delta_zS
Porosity_bed = 0.45                 # assumed
rho_s = 1210                        # particle density: kg/m3 from Kim 2019
viscosity = 8e-6                    # gas viscosity: kg/m/s
Radius_SP = 8e-4                    # pellet particle radius: m from Kim 2019
Cp_g = 119                          # averaged heat capacity of gas: J/mol/K
Cp_s = 920                          # heat capacity of adsorbent: J/kg/K
Cp_a = 30                           # averaged heat capacity of adsorbed phase: J/mol/K

K_W = 16                            # thermal diffusivity of wall: J/m/s/K
rho_W = 8238                        # density of wall: kg/m3
Cp_W = 500                          # heat capacity of wall: J/kg/K
h_W_in = 8.6                        # heat transfer coefficient of inner wall: W/m2/K
h_W_out = 2.5                       # heat transfer coefficient of outer wall: W/m2/K

R = 8.3145                          # gas constant
IC = 7.67                           # isentropic coefficient: k/k-1
Eff = 1.333                         # 1/efficiency of vacuum and compressor
Sigma = 1e-5                        # small number in WENO

EB_1 = K_W / rho_W / Cp_W
EB_2 = 2 * h_W_in * R_inner / (R_out * R_out - R_inner * R_inner) / Cp_W / rho_W
EB_3 = 2 * h_W_out * R_out / (R_out * R_out - R_inner * R_inner) / Cp_W / rho_W
EM = 150 * viscosity * (1 - Porosity_bed) * (1 - Porosity_bed) / (4 * Radius_SP * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed)
FM = 1.75 * (1 - Porosity_bed) / (2 * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed)
M = (1 - Porosity_bed) * R * rho_s / Porosity_bed
CSA = 3.1416 * Diameter * Diameter / 4
rho_bed = (1 - Porosity_bed) * rho_s
ERTI = Eff * Porosity_bed * CSA * IC / 1e6

# %% parameters of 4A zeolites

Qsat_1_PA = 1.7527              # mol/kg
Qsat_2_PA = 0                   # mol/kg
k_1_PA = -16.9856               # exp(k_1_PA) = 4.2e-8: 1/Pa
k_2_PA = 16230 / R              # K
k_3_PA = -50                    # exp(k_3_PA) ~= 0: 1/Pa
k_4_PA = 0                      # K

Qsat_1_PE = 1.1866              # mol/kg
Qsat_2_PE = 0.7656              # mol/kg
k_1_PE = -16.9188               # exp(k_1_PE) = 4.55e-12: 1/Pa
k_2_PE = 20390 / R              # K
k_3_PE = -26.1159               # exp(k_3_PE) = 4.49e-8: 1/Pa
k_4_PE = 58010 / R              # K

kA = 7.56e-6                    # 1/s
kE = 2.24e-2                    # 1/s

MW_PA = 0.0441                  # kg/mol
MW_PE = 0.04208                 # kg/mol

delta_H_A = 15605               # J/mol
delta_H_E = 28265               # J/mol

# %% decision variables

PE_purity_result = []
PE_recovery_result = []
PE_productivity_result = []
Energy_result = []

# %% sampling

for z in range(z_min, z_max):
    # print(z)
    # try:
        P_high = df[z][0]
        P_equ_end = df[z][1]
        P_low = df[z][2]
        delta_t1 = df[z][3]         # pressurization
        delta_t2 = df[z][4]         # adsorption
        delta_t3 = df[z][5]         # equalization
        delta_t4 = df[z][6]         # rinse
        delta_t5 = df[z][7]         # desorption

        P_high = P_high * P_atm
        P_low = P_low * P_atm
        P_abs_end = P_equ_end
        P_rin_end = P_equ_end
        P_abs_end = P_abs_end * P_atm
        P_equ_end = P_equ_end * P_atm
        P_rin_end = P_rin_end * P_atm

        # %% initial
        yE_ini = 0.99
        yA_ini = 1 - yE_ini
        P_ini = P_equ_end           
        T_ini = T_feed
        Tw_ini = T_amb                 
        
        rho_g_ini = P_ini * (MW_PA * yA_ini + MW_PE * yE_ini) / R / T_ini
        qA_ini = Qsat_1_PA * math.exp(k_1_PA + k_2_PA / T_ini) * P_ini * yA_ini / (1 + math.exp(k_1_PA + k_2_PA / T_ini) * P_ini * yA_ini + math.exp(k_1_PE + k_2_PE / T_ini) * P_ini * yE_ini)                + Qsat_2_PA * math.exp(k_3_PA + k_4_PA / T_ini) * P_ini * yA_ini / (1 + math.exp(k_3_PA + k_4_PA / T_ini) * P_ini * yA_ini + math.exp(k_3_PE + k_4_PE / T_ini) * P_ini * yE_ini)
        qE_ini = Qsat_1_PE * math.exp(k_1_PE + k_2_PE / T_ini) * P_ini * yE_ini / (1 + math.exp(k_1_PA + k_2_PA / T_ini) * P_ini * yA_ini + math.exp(k_1_PE + k_2_PE / T_ini) * P_ini * yE_ini)                + Qsat_2_PE * math.exp(k_3_PE + k_4_PE / T_ini) * P_ini * yE_ini / (1 + math.exp(k_3_PA + k_4_PA / T_ini) * P_ini * yA_ini + math.exp(k_3_PE + k_4_PE / T_ini) * P_ini * yE_ini)
        
        u_ini = 0
        qA_Star_ini = qA_ini
        qE_Star_ini = qE_ini
        EB_ini = (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA_ini + qE_ini)) - Cp_g * Porosity_bed * P_ini / R / T_ini
        
        x_ini = [yE_ini, P_ini, T_ini, Tw_ini, qA_ini, qE_ini]
        x_ini = [item for item in x_ini for i in range(n_space)]
        
        z_ini_1 = [yE_feed, P_high, T_feed, Tw_ini, yE_ini, P_ini, T_ini, Tw_ini]
        z_ini_2 = [yE_ini, P_ini, T_ini, u_ini]
        z_ini_2 = [item for item in z_ini_2 for i in range(n_space+1)]
        z_ini_3 = [u_ini, yA_ini, qA_Star_ini, qE_Star_ini, EB_ini]
        z_ini_3 = [item for item in z_ini_3 for i in range(n_space+2)]
        z_ini = z_ini_1 + z_ini_2 + z_ini_3

        # %% iterations
        for s in range(0, 300):

            x_ini = x_ini[0:6*n_space] + [0, 0, 0]
            [x_1, ode_1, z_1, alg_1] = pressurization(n_space, yE_feed, T_feed, P_high)
            dae_1 = {'x':x_1, 'ode':ode_1, 'z':z_1, 'alg':alg_1}
            F = integrator('F', 'idas', dae_1, {'t0':0, 'tf': delta_t1})
            F_sol_1 = F(x0 = x_ini, z0 = z_ini)
            xf_1 = F_sol_1['xf'].elements()
            zf_1 = F_sol_1['zf'].elements()
            PRE_PA_feed = np.float(xf_1[6*n_space])
            PRE_PE_feed = np.float(xf_1[6*n_space + 1])
            PRE_E_input = np.float(xf_1[6*n_space + 2])

            xf_1 = xf_1[0:6*n_space] + [0, 0, 0, 0, 0, 0]
            [x_2, ode_2, z_2, alg_2] = adsorption(n_space, yE_feed, T_feed, P_high, P_abs_end)
            dae_2 = {'x':x_2, 'ode':ode_2, 'z':z_2, 'alg':alg_2}
            F = integrator('F', 'idas', dae_2, {'t0':0, 'tf': delta_t2})
            F_sol_2 = F(x0 = xf_1, z0 = zf_1)
            xf_2 = F_sol_2['xf'].elements()
            zf_2 = F_sol_2['zf'].elements()
            ADS_PA_feed = np.float(xf_2[6*n_space])
            ADS_PE_feed = np.float(xf_2[6*n_space + 1])
            ADS_PA_out = np.float(xf_2[6*n_space + 2])
            ADS_PE_out = np.float(xf_2[6*n_space + 3])
            ADS_E_input = np.float(xf_2[6*n_space + 4])
            ADS_E_out = np.float(xf_2[6*n_space + 5])

            xf_2 = xf_2[0:6*n_space] + [0, 0, 0]
            [x_3, ode_3, z_3, alg_3] = equalization_one(n_space, P_equ_end)
            dae_3 = {'x':x_3, 'ode':ode_3, 'z':z_3, 'alg':alg_3}
            F = integrator('F', 'idas', dae_3, {'t0':0, 'tf': delta_t3})
            F_sol_3 = F(x0 = xf_2, z0 = zf_2)
            xf_3 = F_sol_3['xf'].elements()
            zf_3 = F_sol_3['zf'].elements()
            EQU_PA_out = np.float(xf_3[6 * n_space])
            EQU_PE_out = np.float(xf_3[6 * n_space + 1])
            EQU_E_out = np.float(xf_3[6 * n_space + 2])
            yE_equ_out = np.abs(EQU_PE_out) / np.abs(EQU_PA_out + EQU_PE_out)

            xf_3 = xf_3[0:6*n_space] + [0, 0, 0, 0, 0, 0]
            [x_4, ode_4, z_4, alg_4] = rinse(n_space, yE_rinse, T_feed, P_high, P_rin_end)
            dae_4 = {'x':x_4, 'ode':ode_4, 'z':z_4, 'alg':alg_4}
            F = integrator('F', 'idas', dae_4, {'t0':0, 'tf': delta_t4})
            F_sol_4 = F(x0 = xf_3, z0 = zf_3)
            xf_4 = F_sol_4['xf'].elements()
            zf_4 = F_sol_4['zf'].elements()
            RIN_PA_feed = np.float(xf_4[6 * n_space])
            RIN_PE_feed = np.float(xf_4[6 * n_space + 1]) 
            RIN_PA_out = np.float(xf_4[6 * n_space + 2])
            RIN_PE_out = np.float(xf_4[6 * n_space + 3])
            RIN_E_input = np.float(xf_4[6 * n_space + 4])
            RIN_E_out = np.float(xf_4[6 * n_space + 5])
            yE_rin_out = RIN_PE_out / (RIN_PA_out + RIN_PE_out)

            xf_4 = xf_4[0:6*n_space] + [0, 0, 0]
            [x_5, ode_5, z_5, alg_5] = desorption(n_space, P_low)
            dae_5 = {'x':x_5, 'ode':ode_5, 'z':z_5, 'alg':alg_5}
            F = integrator('F', 'idas', dae_5, {'t0':0, 'tf': delta_t5})
            F_sol_5 = F(x0 = xf_4, z0 = zf_4)
            xf_5 = F_sol_5['xf'].elements()
            zf_5 = F_sol_5['zf'].elements()
            DES_PA_out = np.float(xf_5[6 * n_space])
            DES_PE_out = np.float(xf_5[6 * n_space + 1])
            DES_E_out = np.float(xf_5[6 * n_space + 2])
            
            xf_5 = xf_5[0:6*n_space] + [0, 0]
            for i in range(0, 1000):
                delta_eqn = 0.05 + 0.05 * i
                [x_6, ode_6, z_6, alg_6] = equalization_two(n_space, yE_equ_out, T_feed, P_equ_end)
                dae_6 = {'x':x_6, 'ode':ode_6, 'z':z_6, 'alg':alg_6}
                F = integrator('F', 'idas', dae_6, {'t0':0, 'tf': delta_eqn})
                F_sol_6 = F(x0 = xf_5, z0 = zf_5)
                xf_6 = F_sol_6['xf'].elements()
                EQN_PA_feed = np.float(xf_6[6 * n_space])
                EQN_PE_feed = np.float(xf_6[6 * n_space + 1])
                if EQN_PA_feed + EQN_PE_feed >= np.abs(EQU_PA_out + EQU_PE_out):
                    break

            xf_6 = xf_6[0:6*n_space] + [0, 0]
            zf_6 = F_sol_6['zf'].elements()    
            for i in range(0, 1000):
                delta_ref = 0.05 + 0.05 * i
                [x_7, ode_7, z_7, alg_7] = reflux(n_space, yE_rin_out, T_feed, P_rin_end)
                dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_ref})
                F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                xf_7 = F_sol_7['xf'].elements()
                REF_PA_feed = np.float(xf_7[6 * n_space])
                REF_PE_feed = np.float(xf_7[6 * n_space + 1])
                if REF_PA_feed + REF_PE_feed >= np.abs(RIN_PA_out + RIN_PE_out):
                    break
            xf_7 = xf_7[0:6*n_space]
            zf_7 = F_sol_7['zf'].elements()

            xstart = np.asarray(x_ini[0:6*n_space])
            xend = np.asarray(xf_7)
            error = np.max(np.abs((xstart - xend) / (xstart + 1e-10)))
            abs_error = np.max(np.abs(xstart - xend))
            print("%d th iter, max error: %f, %f" %(s, error, abs_error))
            if error <= 0.01:
                break
            else:        
                x_ini = xf_7
                z_ini = zf_7
                z_ini[0] = yE_feed
                z_ini[1] = P_high
                z_ini[2] = T_feed
                z_ini[3] = T_amb               

# %% performance
        if P_high >= 2 * P_atm:
            E_1_feed = ERTI * PRE_E_input * ((P_high / (2 * P_atm))**(1/IC) - 1) 
            E_2_feed = ERTI * ADS_E_input * ((P_high / (2 * P_atm))**(1/IC) - 1)
        else:
            E_1_feed = 0
            E_2_feed = 0

        if P_equ_end >= P_atm:
            E_2_out = 0
            E_3_out = 0
            E_4_out = 0
        else:
            E_2_out = ERTI * ADS_E_out * ((P_atm / P_equ_end)**(1/IC) - 1)
            E_3_out = ERTI * EQU_E_out * ((P_atm / P_equ_end)**(1/IC) - 1)
            E_4_out = ERTI * RIN_E_out * ((P_atm / P_equ_end)**(1/IC) - 1)
    
        E_4_feed = ERTI * RIN_E_input * ((P_high / P_low)**(1/IC) - 1)
        E_5_out = ERTI * np.abs(DES_E_out) * ((P_atm / P_low)**(1/IC) - 1)

        PE_purity = (DES_PE_out) / (DES_PA_out + DES_PE_out)
        PE_recovery = (DES_PE_out - RIN_PE_feed) / (PRE_PE_feed + ADS_PE_feed)
        PE_productivity = (DES_PE_out - RIN_PE_feed) / (delta_t1 + delta_t2 + delta_t3 + delta_t4 + delta_t5 + delta_eqn + delta_ref + 1e-10)
        Energy = (E_1_feed + E_2_feed + E_2_out + E_3_out + E_4_out + E_4_feed + E_5_out) / (DES_PE_out - RIN_PE_feed + 1e-10)
        
        PE_purity_result.append(PE_purity)
        PE_recovery_result.append(PE_recovery)
        PE_productivity_result.append(PE_productivity)
        Energy_result.append(Energy)
    # except:
    #     PE_purity_result.append(np.nan)
    #     PE_recovery_result.append(np.nan)
    #     PE_productivity_result.append(np.nan)
    #     Energy_result.append(np.nan)
        
# # %% save

# PE_purity_result = np.asarray(PE_purity_result).reshape(-1, 1)
# PE_recovery_result = np.asarray(PE_recovery_result).reshape(-1, 1)
# PE_productivity_result = np.asarray(PE_productivity_result).reshape(-1, 1)
# Energy_result = np.asarray(Energy_result).reshape(-1, 1)    

# df = df[z_min:z_max, :]
# Result = np.hstack([PE_purity_result, PE_recovery_result, PE_productivity_result, Energy_result])
# Result = np.hstack([df, Result])

# Header = np.asarray(['P_high', 'P_equ', 'P_low', 'delta_t1',
#                       'delta_t2', 'delta_t3', 'delta_t4', 'delta_t5', 'PE Purity',
#                       'PE Recovery', 'PE Productivity', 'Energy']).reshape([1, -1])

# Result = np.vstack([Header, Result])

# # %% output to excel

# Result = pd.DataFrame(Result)

# save_name = int(z_min / gap)
# save_name = 'Result_' + str(save_name) + '.xlsx'
# Result.to_excel(save_name, header = False, index = False)                