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

rank = 2

# z_min = 30 * (rank - 1) + 3
# z_max = 30 * (rank - 1) + 4
z_min = 0
z_max = 1

df = pd.read_excel('Test_2_CS2.xlsx', header = 0)
df = np.array(df)

# %% steps
def pressurization(n_space, yN_feed, yO_feed, T_feed, P_high):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    Input   = SX.sym('Input',   1)                  # state

    # left boundary conditions
    AE_L1 = yN[0] - yN_feed
    AE_L2 = yO[0] - yO_feed
    AE_L3 = P[0] - P_high
    AE_L4 = T[0] - T_feed
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_half[n_space]
    AE_R2 = yO[n_space+1] - yO_half[n_space]
    AE_R3 = P[n_space+1] - P_half[n_space]
    AE_R4 = T[n_space+1] - T_half[n_space]
    AE_R5 = u[n_space+1] - 0
    AE_R6 = EB[n_space+1] - 0

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN_feed
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i+1] + 0.5 * yN[i]) * 0.6667 / (Sigma + yN[i+1] - yN[i])**4 + (1.5 * yN[i]-0.5 * yN[i-1]) * 0.3333 / (Sigma + yN[i] - yN[i-1])**4) / (0.6667 / (Sigma + yN[i+1] - yN[i])**4 + 0.3333 / (Sigma + yN[i] - yN[i-1])**4)) 
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN[n_space])

    AE_yO_half = yO_half[0] - yO_feed
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i+1] + 0.5 * yO[i]) * 0.6667 / (Sigma + yO[i+1] - yO[i])**4 + (1.5 * yO[i]-0.5 * yO[i-1]) * 0.3333 / (Sigma + yO[i] - yO[i-1])**4) / (0.6667 / (Sigma + yO[i+1] - yO[i])**4 + 0.3333 / (Sigma + yO[i] - yO[i-1])**4)) 
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # mass balance
    ODE_Input = u[0] * P[0] * CSA * Porosity_bed / R / T[0]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     Input)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_Input)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def adsorption_one(n_space, yN_feed, yO_feed, T_feed, P_high, P_ads_end):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    Input   = SX.sym('Input',   1)                  # state
    Output  = SX.sym('Output',  1)                  # output state
    PN_out  = SX.sym('PN_out',  1)                  # output state
    PO_out  = SX.sym('PO_out',  1)                  # output state

    # left boundary conditions
    AE_L1 = yN[0] - yN_feed
    AE_L2 = yO[0] - yO_feed
    AE_L3 = P[0] - P_high
    AE_L4 = T[0] - T_feed
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_half[n_space]
    AE_R2 = yO[n_space+1] - yO_half[n_space]
    AE_R3 = P[n_space+1] - P_ads_end
    AE_R4 = T[n_space+1] - T_half[n_space]
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN_feed
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i+1] + 0.5 * yN[i]) * 0.6667 / (Sigma + yN[i+1] - yN[i])**4 + (1.5 * yN[i]-0.5 * yN[i-1]) * 0.3333 / (Sigma + yN[i] - yN[i-1])**4) / (0.6667 / (Sigma + yN[i+1] - yN[i])**4 + 0.3333 / (Sigma + yN[i] - yN[i-1])**4)) 
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN[n_space])

    AE_yO_half = yO_half[0] - yO_feed
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i+1] + 0.5 * yO[i]) * 0.6667 / (Sigma + yO[i+1] - yO[i])**4 + (1.5 * yO[i]-0.5 * yO[i-1]) * 0.3333 / (Sigma + yO[i] - yO[i-1])**4) / (0.6667 / (Sigma + yO[i+1] - yO[i])**4 + 0.3333 / (Sigma + yO[i] - yO[i-1])**4)) 
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4))

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # mass balance
    ODE_Input = u[0] * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_Output = u[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PN_out = u[n_space+1] * P[n_space+1] * yN[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PO_out = u[n_space+1] * P[n_space+1] * yO[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     Input,     Output,     PN_out,     PO_out)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_Input, ODE_Output, ODE_PN_out, ODE_PO_out)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def adsorption_two(n_space, yN_feed, yO_feed, T_feed, P_high, P_adp_end):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    Input   = SX.sym('Input',   1)                  # state
    Output  = SX.sym('Output',  1)                  # state
    PN_out  = SX.sym('PN_out',  1)                  # state
    PO_out  = SX.sym('PO_out',  1)                  # state

    # left boundary conditions
    AE_L1 = yN[0] - yN_feed
    AE_L2 = yO[0] - yO_feed
    AE_L3 = P[0] - P_high
    AE_L4 = T[0] - T_feed
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_half[n_space]
    AE_R2 = yO[n_space+1] - yO_half[n_space]
    AE_R3 = P[n_space+1] - P_adp_end
    AE_R4 = T[n_space+1] - T_half[n_space]
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN_feed
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i+1] + 0.5 * yN[i]) * 0.6667 / (Sigma + yN[i+1] - yN[i])**4 + (1.5 * yN[i]-0.5 * yN[i-1]) * 0.3333 / (Sigma + yN[i] - yN[i-1])**4) / (0.6667 / (Sigma + yN[i+1] - yN[i])**4 + 0.3333 / (Sigma + yN[i] - yN[i-1])**4)) 
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN[n_space])

    AE_yO_half = yO_half[0] - yO_feed
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i+1] + 0.5 * yO[i]) * 0.6667 / (Sigma + yO[i+1] - yO[i])**4 + (1.5 * yO[i]-0.5 * yO[i-1]) * 0.3333 / (Sigma + yO[i] - yO[i-1])**4) / (0.6667 / (Sigma + yO[i+1] - yO[i])**4 + 0.3333 / (Sigma + yO[i] - yO[i-1])**4)) 
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO[n_space])

    AE_P_half = P_half[0] - P_high
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4))

    AE_T_half = T_half[0] - T_feed
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # material balances
    ODE_Input = u[0] * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_Output = u[n_space+1] * P[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PN_out = u[n_space+1] * P[n_space+1] * yN[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PO_out = u[n_space+1] * P[n_space+1] * yO[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     Input,     Output,     PN_out,     PO_out)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_Input, ODE_Output, ODE_PN_out, ODE_PO_out)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def equalization_one(n_space, P_equ_end):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    PN_out  = SX.sym('PN_out',  1)                  # state
    PO_out  = SX.sym('PO_out',  1)                  # state
    PA_out  = SX.sym('PA_out',  1)                  # state
    
    # left boundary conditions
    AE_L1 = yN[0] - yN[1]
    AE_L2 = yO[0] - yO[1]
    AE_L3 = P[0] - P[1]
    AE_L4 = T[0] - T[1]
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_half[n_space]
    AE_R2 = yO[n_space+1] - yO_half[n_space]
    AE_R3 = P[n_space+1] - P_equ_end
    AE_R4 = T[n_space+1] - T_half[n_space]
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN[1]
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i+1] + 0.5 * yN[i]) * 0.6667 / (Sigma + yN[i+1] - yN[i])**4 + (1.5 * yN[i]-0.5 * yN[i-1]) * 0.3333 / (Sigma + yN[i] - yN[i-1])**4) / (0.6667 / (Sigma + yN[i+1] - yN[i])**4 + 0.3333 / (Sigma + yN[i] - yN[i-1])**4)) 
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN[n_space])

    AE_yO_half = yO_half[0] - yO[1]
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i+1] + 0.5 * yO[i]) * 0.6667 / (Sigma + yO[i+1] - yO[i])**4 + (1.5 * yO[i]-0.5 * yO[i-1]) * 0.3333 / (Sigma + yO[i] - yO[i-1])**4) / (0.6667 / (Sigma + yO[i+1] - yO[i])**4 + 0.3333 / (Sigma + yO[i] - yO[i-1])**4)) 
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO[n_space])

    AE_P_half = P_half[0] - P[1]
    for i in range(1, n_space + 1):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i+1] + 0.5 * P[i]) * 0.6667 / (Sigma + P[i+1] - P[i])**4 + (1.5 * P[i] - 0.5 * P[i-1]) * 0.3333 / (Sigma + P[i] - P[i-1])**4) / (0.6667 / (Sigma + P[i+1] - P[i])**4 + 0.3333 / (Sigma + P[i] - P[i-1])**4))

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i+1] + 0.5 * T[i]) * 0.6667 / (Sigma + T[i+1] - T[i])**4 + (1.5 * T[i] - 0.5 * T[i-1]) * 0.3333 / (Sigma + T[i] - T[i-1])**4) / (0.6667 / (Sigma + T[i+1] - T[i])**4 + 0.3333 / (Sigma + T[i] - T[i-1])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # material balances
    ODE_PN_out = u[n_space+1] * P[n_space+1] * yN[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PO_out = u[n_space+1] * P[n_space+1] * yO[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PA_out = u[n_space+1] * P[n_space+1] * yA[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     PN_out,     PO_out,     PA_out)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_PN_out, ODE_PO_out, ODE_PA_out)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def desorption_one(n_space, P_low):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    Output  = SX.sym('Output',  1)                  # state
    PN_out  = SX.sym('PN_out',  1)                  # state
    PO_out  = SX.sym('PO_out',  1)                  # state

    # left boundary conditions
    AE_L1 = yN[0] - yN_half[0]
    AE_L2 = yO[0] - yO_half[0]
    AE_L3 = P[0] - P_low
    AE_L4 = T[0] - T_half[0]
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_half[n_space]
    AE_R2 = yO[n_space+1] - yO_half[n_space]
    AE_R3 = P[n_space+1] - P_half[n_space]
    AE_R4 = T[n_space+1] - T_half[n_space]
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0    

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN[1]
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i] + 0.5 * yN[i+1]) * 0.6667 / (Sigma + yN[i] - yN[i+1])**4 + (1.5 * yN[i+1] - 0.5 * yN[i+2]) * 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4) / (0.6667 / (Sigma + yN[i] - yN[i+1])**4 + 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4))
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN[n_space])

    AE_yO_half = yO_half[0] - yO[1]
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i] + 0.5 * yO[i+1]) * 0.6667 / (Sigma + yO[i] - yO[i+1])**4 + (1.5 * yO[i+1] - 0.5 * yO[i+2]) * 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4) / (0.6667 / (Sigma + yO[i] - yO[i+1])**4 + 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4))
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO[n_space])

    AE_P_half = P_half[0] - P_low
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i] + 0.5 * P[i+1]) * 0.6667 / (Sigma + P[i] - P[i+1])**4 + (1.5 * P[i+1] - 0.5 * P[i+2]) * 0.3333 / (Sigma + P[i+1] - P[i+2])**4) / (0.6667 / (Sigma + P[i] - P[i+1])**4 + 0.3333 / (Sigma + P[i+1] - P[i+2])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P[n_space])

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i] + 0.5 * T[i+1]) * 0.6667 / (Sigma + T[i] - T[i+1])**4 + (1.5 * T[i+1] - 0.5 * T[i+2]) * 0.3333 / (Sigma + T[i+1] - T[i+2])**4) / (0.6667 / (Sigma + T[i] - T[i+1])**4 + 0.3333 / (Sigma + T[i+1] - T[i+2])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T[n_space])

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])
    
    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # mass balance
    ODE_Output = - u[0] * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_PN_out = - u[0] * P[0] * yN[0] * CSA * Porosity_bed / R / T[0]
    ODE_PO_out = - u[0] * P[0] * yO[0] * CSA * Porosity_bed / R / T[0]
    
    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     Output,     PN_out,     PO_out)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_Output, ODE_PN_out, ODE_PO_out)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def purge(n_space, yN_pur_feed, yO_pur_feed, T_pur_feed, P_low, P_pur_feed):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    PN_feed = SX.sym('PN_feed', 1)                  # state
    PO_feed = SX.sym('PO_feed', 1)                  # state
    PA_feed = SX.sym('PA_feed', 1)                  # state
    Output  = SX.sym('Output',  1)                  # state
    PN_out  = SX.sym('PN_out',  1)                  # state
    PO_out  = SX.sym('PO_out',  1)                  # state
   
    # left boundary conditions
    AE_L1 = yN[0] - yN_half[0]
    AE_L2 = yO[0] - yO_half[0]
    AE_L3 = P[0] - P_low
    AE_L4 = T[0] - T_half[0]
    AE_L5 = u[0] - u_half[0]
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_pur_feed
    AE_R2 = yO[n_space+1] - yO_pur_feed
    AE_R3 = P[n_space+1] - P_pur_feed
    AE_R4 = T[n_space+1] - T_pur_feed
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0    

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN[1]
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i] + 0.5 * yN[i+1]) * 0.6667 / (Sigma + yN[i] - yN[i+1])**4 + (1.5 * yN[i+1] - 0.5 * yN[i+2]) * 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4) / (0.6667 / (Sigma + yN[i] - yN[i+1])**4 + 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4))
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN_pur_feed)

    AE_yO_half = yO_half[0] - yO[1]
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i] + 0.5 * yO[i+1]) * 0.6667 / (Sigma + yO[i] - yO[i+1])**4 + (1.5 * yO[i+1] - 0.5 * yO[i+2]) * 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4) / (0.6667 / (Sigma + yO[i] - yO[i+1])**4 + 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4))
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO_pur_feed)

    AE_P_half = P_half[0] - P_low
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i] + 0.5 * P[i+1]) * 0.6667 / (Sigma + P[i] - P[i+1])**4 + (1.5 * P[i+1] - 0.5 * P[i+2]) * 0.3333 / (Sigma + P[i+1] - P[i+2])**4) / (0.6667 / (Sigma + P[i] - P[i+1])**4 + 0.3333 / (Sigma + P[i+1] - P[i+2])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P_pur_feed)

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i] + 0.5 * T[i+1]) * 0.6667 / (Sigma + T[i] - T[i+1])**4 + (1.5 * T[i+1] - 0.5 * T[i+2]) * 0.3333 / (Sigma + T[i+1] - T[i+2])**4) / (0.6667 / (Sigma + T[i] - T[i+1])**4 + 0.3333 / (Sigma + T[i+1] - T[i+2])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T_pur_feed)

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # material balance
    ODE_PN_feed = - u[n_space+1] * P[n_space+1] * yN[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PO_feed = - u[n_space+1] * P[n_space+1] * yO[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PA_feed = - u[n_space+1] * P[n_space+1] * yA[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_Output  = - u[0] * P[0] * CSA * Porosity_bed / R / T[0]
    ODE_PN_out = - u[0] * P[0] * yN[0] * CSA * Porosity_bed / R / T[0]
    ODE_PO_out = - u[0] * P[0] * yO[0] * CSA * Porosity_bed / R / T[0]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     PN_feed,     PO_feed,     PA_feed,     Output,     PN_out,     PO_out)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_PN_feed, ODE_PO_feed, ODE_PA_feed, ODE_Output, ODE_PN_out, ODE_PO_out)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

def equalization_two(n_space, yN_eqn_feed, yO_eqn_feed, T_eqn_feed, P_eqn_feed):

    # variables
    yN      = SX.sym('yN',      n_space+2)          # time differential
    yO      = SX.sym('yO',      n_space+2)          # time differential
    P       = SX.sym('P',       n_space+2)          # time differential, unit Pa
    T       = SX.sym('T',       n_space+2)          # time differential, unit K
    u       = SX.sym('u',       n_space+2)          # time differential, unit K
    yN_half = SX.sym('yN_half', n_space+1)          # wall
    yO_half = SX.sym('yO_half', n_space+1)          # wall
    P_half  = SX.sym('P_half',  n_space+1)          # wall
    T_half  = SX.sym('T_half',  n_space+1)          # wall
    u_half  = SX.sym('u_half',  n_space+1)          # wall
    qN      = SX.sym('qN',      n_space+2)          # time differential, unit mol/kg
    qO      = SX.sym('qO',      n_space+2)          # time differential, unit mol/kg
    qA      = SX.sym('qA',      n_space+2)          # time differential, unit mol/kg
    qN_Star = SX.sym('qN_Star', n_space+2)          # state, unit mol/kg
    qO_Star = SX.sym('qO_Star', n_space+2)          # state, unit mol/kg
    qA_Star = SX.sym('qA_Star', n_space+2)          # state, unit mol/kg
    yA      = SX.sym('yA',      n_space+2)          # state
    EB      = SX.sym('EB',      n_space+2)          # state
    PN_feed = SX.sym('PN_feed', 1)                  # state
    PO_feed = SX.sym('PO_feed', 1)                  # state
    PA_feed = SX.sym('PA_feed', 1)                  # state

    # left boundary conditions
    AE_L1 = yN[0] - yN_half[0]
    AE_L2 = yO[0] - yO_half[0]
    AE_L3 = P[0] - P_half[0]
    AE_L4 = T[0] - T_half[0]
    AE_L5 = u[0] - 0
    AE_L6 = EB[0] - 0

    # right boundary conditions
    AE_R1 = yN[n_space+1] - yN_eqn_feed
    AE_R2 = yO[n_space+1] - yO_eqn_feed
    AE_R3 = P[n_space+1] - P_eqn_feed
    AE_R4 = T[n_space+1] - T_eqn_feed
    AE_R5 = u[n_space+1] - u_half[n_space]
    AE_R6 = EB[n_space+1] - 0    

    # WENO wall + boundary conditions
    AE_yN_half = yN_half[0] - yN[1]
    for i in range(1, n_space):
        AE_yN_half = vertcat(AE_yN_half, yN_half[i] - ((0.5 * yN[i] + 0.5 * yN[i+1]) * 0.6667 / (Sigma + yN[i] - yN[i+1])**4 + (1.5 * yN[i+1] - 0.5 * yN[i+2]) * 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4) / (0.6667 / (Sigma + yN[i] - yN[i+1])**4 + 0.3333 / (Sigma + yN[i+1] - yN[i+2])**4))
    AE_yN_half = vertcat(AE_yN_half, yN_half[n_space] - yN_eqn_feed)

    AE_yO_half = yO_half[0] - yO[1]
    for i in range(1, n_space):
        AE_yO_half = vertcat(AE_yO_half, yO_half[i] - ((0.5 * yO[i] + 0.5 * yO[i+1]) * 0.6667 / (Sigma + yO[i] - yO[i+1])**4 + (1.5 * yO[i+1] - 0.5 * yO[i+2]) * 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4) / (0.6667 / (Sigma + yO[i] - yO[i+1])**4 + 0.3333 / (Sigma + yO[i+1] - yO[i+2])**4))
    AE_yO_half = vertcat(AE_yO_half, yO_half[n_space] - yO_eqn_feed)

    AE_P_half = P_half[0] - P[1]
    for i in range(1, n_space):
        AE_P_half = vertcat(AE_P_half, P_half[i] - ((0.5 * P[i] + 0.5 * P[i+1]) * 0.6667 / (Sigma + P[i] - P[i+1])**4 + (1.5 * P[i+1] - 0.5 * P[i+2]) * 0.3333 / (Sigma + P[i+1] - P[i+2])**4) / (0.6667 / (Sigma + P[i] - P[i+1])**4 + 0.3333 / (Sigma + P[i+1] - P[i+2])**4))
    AE_P_half = vertcat(AE_P_half, P_half[n_space] - P_eqn_feed)

    AE_T_half = T_half[0] - T[1]
    for i in range(1, n_space):
        AE_T_half = vertcat(AE_T_half, T_half[i] - ((0.5 * T[i] + 0.5 * T[i+1]) * 0.6667 / (Sigma + T[i] - T[i+1])**4 + (1.5 * T[i+1] - 0.5 * T[i+2]) * 0.3333 / (Sigma + T[i+1] - T[i+2])**4) / (0.6667 / (Sigma + T[i] - T[i+1])**4 + 0.3333 / (Sigma + T[i+1] - T[i+2])**4))
    AE_T_half = vertcat(AE_T_half, T_half[n_space] - T_eqn_feed)

    AE_u_half = (P[1] - P[0]) / delta_z + EM * u_half[0] + FM * u_half[0] * fabs(u_half[0]) * P_half[0] * (MW_PA * (1 - yN_half[0] - yO_half[0]) + MW_PN * yN_half[0] + MW_PO * yO_half[0]) / R / T_half[0]
    for i in range(1, n_space + 1):
        AE_u_half = vertcat(AE_u_half, (P[i+1] - P[i]) / delta_z + EM * u_half[i] + FM * u_half[i] * fabs(u_half[i]) * P_half[i] * (MW_PA * (1 - yN_half[i] - yO_half[i]) + MW_PN * yN_half[i] + MW_PO * yO_half[i]) / R / T_half[i])

    # sum of mol fraction
    AE_yA = yN + yO + yA - 1

    # Isotherm and LDF equations
    AE_qN_Star = qN_Star - P * yN * exp(Q_PN + b_PN / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qO_Star = qO_Star - P * yO * exp(Q_PO + b_PO / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    AE_qA_Star = qA_Star - P * yA * exp(Q_PA + b_PA / T) / (1e5 + P * yN * exp(Q_PN_b + b_PN_b / T) + P * yO * exp(Q_PO_b + b_PO_b / T) + P * yA * exp(Q_PA_b + b_PA_b / T))
    ODE_qN = kN * (qN_Star - qN)
    ODE_qO = kO * (qO_Star - qO)
    ODE_qA = kA * (qA_Star - qA)

    # single component mass balance
    ODE_yN = DxZS * ((yN[2] + yN[0] - 2 * yN[1]) - (yN_half[1] - yN_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yN_half[1] - yN_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yN_half[1] - yN_half[0]) / delta_z - M * T[1] * (ODE_qN[1] - yN[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yN = vertcat(ODE_yN, DxZS * ((yN[i+1] + yN[i-1] - 2 * yN[i]) - (yN_half[i] - yN_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yN_half[i] - yN_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yN_half[i] - yN_half[i-1]) / delta_z - M * T[i] * (ODE_qN[i] - yN[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    ODE_yO = DxZS * ((yO[2] + yO[0] - 2 * yO[1]) - (yO_half[1] - yO_half[0]) * (T_half[1] - T_half[0]) / T[1] + (yO_half[1] - yO_half[0]) * (P_half[1] - P_half[0]) / P[1]) - u[1] * (yO_half[1] - yO_half[0]) / delta_z - M * T[1] * (ODE_qO[1] - yO[1] * (ODE_qN[1] + ODE_qO[1] + ODE_qA[1])) / P[1]
    for i in range(2, n_space + 1):
        ODE_yO = vertcat(ODE_yO, DxZS * ((yO[i+1] + yO[i-1] - 2 * yO[i]) - (yO_half[i] - yO_half[i-1]) * (T_half[i] - T_half[i-1]) / T[i] + (yO_half[i] - yO_half[i-1]) * (P_half[i] - P_half[i-1]) / P[i]) - u[i] * (yO_half[i] - yO_half[i-1]) / delta_z - M * T[i] * (ODE_qO[i] - yO[i] * (ODE_qN[i] + ODE_qO[i] + ODE_qA[i])) / P[i])

    # dominator in energy balance
    AE_EB = EB[1] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[1] + qN[1] + qO[1])) - Cp_g * Porosity_bed * P[1] / R / T[1]
    for i in range(2, n_space + 1):
        AE_EB = vertcat(AE_EB, EB[i] - (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qA[i] + qN[i] + qO[i])) - Cp_g * Porosity_bed * P[i] / R / T[i])
    
    # column energy balance
    ODE_T = (K_z * (T[2] + T[0] - 2 * T[1]) / delta_zS - Cp_g * Porosity_bed * u[1] * P[1] * (T_half[1] - T_half[0]) / delta_z / R / T[1] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[1] + deltaH_PN * ODE_qN[1] + deltaH_PO * ODE_qO[1]) - 2 * h_W_in * (T[1] - T_amb) / R_inner) / EB[1]
    for i in range(2, n_space + 1):
        ODE_T = vertcat(ODE_T, (K_z * (T[i+1] + T[i-1] - 2 * T[i]) / delta_zS - Cp_g * Porosity_bed * u[i] * P[i] * (T_half[i] - T_half[i-1]) / delta_z / R / T[i] + (1 - Porosity_bed) * (Cp_g-Cp_a) * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + (1 - Porosity_bed) * (deltaH_PA * ODE_qA[i] + deltaH_PN * ODE_qN[i] + deltaH_PO * ODE_qO[i]) - 2 * h_W_in * (T[i] - T_amb) / R_inner) / EB[i])

    # Total mass balance
    ODE_P = - (P_half[1] * u_half[1] - P_half[0] * u_half[0]) / delta_z - M * T[1] * (ODE_qA[1] + ODE_qN[1] + ODE_qO[1]) + P[1] * (ODE_T[0] + u[1] * (T_half[1] - T_half[0]) / delta_z) / T[1]
    for i in range(2, n_space + 1):
        ODE_P = vertcat(ODE_P, - (P_half[i] * u_half[i] - P_half[i-1] * u_half[i-1]) / delta_z - M * T[i] * (ODE_qA[i] + ODE_qN[i] + ODE_qO[i]) + P[i] * (ODE_T[i-1] + u[i] * (T_half[i] - T_half[i-1]) / delta_z) / T[i])

    # Ergun equation
    AE_u = (P_half[1] - P_half[0]) / delta_z + EM * u[1] + FM * u[1] * fabs(u[1]) * P[1] * (MW_PA * yA[1] + MW_PN * yN[1] + MW_PO * yO[1]) / R / T[1]
    for i in range(2, n_space + 1):
        AE_u = vertcat(AE_u, (P_half[i] - P_half[i-1]) / delta_z + EM * u[i] + FM * u[i] * fabs(u[i]) * P[i] * (MW_PA * yA[i] + MW_PN * yN[i] + MW_PO * yO[i]) / R / T[i])

    # material balances
    ODE_PN_feed = - u[n_space+1] * P[n_space+1] * yN[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PO_feed = - u[n_space+1] * P[n_space+1] * yO[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]
    ODE_PA_feed = - u[n_space+1] * P[n_space+1] * yA[n_space+1] * CSA * Porosity_bed / R / T[n_space+1]

    x = vertcat(    yN[1:n_space+1], yO[1:n_space+1], P[1:n_space+1], T[1:n_space+1], qN[1:n_space+1],     qO[1:n_space+1],     qA[1:n_space+1],     PN_feed,     PO_feed,     PA_feed)
    ode = vertcat(  ODE_yN,          ODE_yO,          ODE_P,          ODE_T,          ODE_qN[1:n_space+1], ODE_qO[1:n_space+1], ODE_qA[1:n_space+1], ODE_PN_feed, ODE_PO_feed, ODE_PA_feed)

    z = vertcat(    yN[0], yO[0], P[0],  T[0],  yN[n_space+1], yO[n_space+1], P[n_space+1], T[n_space+1], u,                yN_half,    yO_half,    P_half,    T_half,    u_half,    yA,    qN_Star,    qO_Star,    qA_Star,    EB)
    alg = vertcat(  AE_L1, AE_L2, AE_L3, AE_L4, AE_R1,         AE_R2,         AE_R3,        AE_R4,        AE_L5,AE_u,AE_R5, AE_yN_half, AE_yO_half, AE_P_half, AE_T_half, AE_u_half, AE_yA, AE_qN_Star, AE_qO_Star, AE_qA_Star, AE_L6,AE_EB,AE_R6)

    return [x, ode, z, alg]

# %% parameters

n_space = 30                        # number of space interval
L = 2                               # bed lengthL: m
delta_z = L / n_space               # space interval: m                            
delta_zS = delta_z * delta_z
Diameter = 3.5                      # bed diameter: m
R_inner = Diameter / 2              # bed radius: m
R = 8.3145                          # gas constant
P_bar = 1e5                         # 1bar  = 1e5 Pa

yN_feed = 0.78                      # air composition
yO_feed = 0.21                      # air composition
yA_feed = 0.01                      # air composition
T_feed = 298.15                     # feed temperature: K
T_amb = 298.15                      # ambinent temperature: K

Dx = 5e-5                           # axial dispersion coefficient: m2/s
DxZS = Dx / delta_zS
K_z = 0.044                         # thermal conductivity: W/m/K
Cp_g = 30                           # heat capacity of gas: J/mol/K
Cp_a = 28.85                        # heat capacity of air: J/mol/K
Cp_s = 1210                         # heat capacity of solid: J/kg/K
rho_s = 630                         # particle density: kg/m3
viscosity = 8e-6
Radius_SP = 8.5e-4                  # particle radius: m
Porosity_bed = 0.36                 # bed porosity
h_W_in = 0.3                        # heat transfer coefficient: W/m2/K
Sigma = 1e-6                    

EM = 150 * viscosity * (1 - Porosity_bed) * (1 - Porosity_bed) / (4 * Radius_SP * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed)
FM = 1.75 * (1 - Porosity_bed) / (2 * Radius_SP * Porosity_bed * Porosity_bed * Porosity_bed)
M = (1 - Porosity_bed) * R * rho_s / Porosity_bed
CSA = 3.1416 * R_inner * R_inner

# %% isotherm parameters

Q_PN = -9.5518
b_PN = 2910
Q_PO = -7.2845
b_PO = 1567
Q_PA = -7.3613
b_PA = 1634

Q_PN_b = -5.9666
b_PN_b = 1612
Q_PO_b = -5.3763
b_PO_b = 441.3
Q_PA_b = -5.4433
b_PA_b = 450.6

kN = 1
kO = 2.5
kA = 2.5

deltaH_PN = 23430       
deltaH_PO = 13220
deltaH_PA = 12650

MW_PN = 0.028
MW_PO = 0.032
MW_PA = 0.04

# %% decision variables

P_pur_feed_final_result = []
P_eqn_feed_final_result = []
PO_purity_result = []
PO_productivity_result = []

PRE_Input_result = []
ADS_Input_result = []
ADS_Output_result = []
ADS_PO_out_result = []
ADP_Input_result = []
ADP_Output_result = []
EQU_Output_result = []
DES_Output_result = []
PUR_Input_result = []
PUR_Output_result = []
EQN_Input_result = []

# %% sampling

for z in range(z_min, z_max):
    # try:
        P_high    = df[z][0]
        P_ads_end = df[z][1]
        P_adp_end = df[z][2]
        P_equ_end = df[z][3]
        P_low     = df[z][4]
        delta_t1  = df[z][5]     # pressurization, desorption_one
        delta_t2  = df[z][6]     # adsorption_one, desorption_two
        delta_t3  = df[z][7]     # adsorption_two
        delta_t4  = df[z][8]     # equalization_one
        
        P_high    = P_high * P_bar
        P_ads_end = P_ads_end * P_bar
        P_adp_end = P_adp_end * P_bar
        P_equ_end = P_equ_end * P_bar
        P_low     = P_low * P_bar
    
        # initial
        yN_ini = 0.19
        yO_ini = 0.8
        yA_ini = 0.01
        P_ini  = P_low * 1.5
        T_ini  = T_amb
    
        rho_g_ini = P_ini * (MW_PN * yN_ini + MW_PO * yO_ini + MW_PA * yA_ini) / R / T_ini
        qN_ini = P_ini * yN_ini * math.exp(Q_PN + b_PN / T_ini) / (1e5 + P_ini * yN_ini * math.exp(Q_PN_b + b_PN_b / T_ini) + P_ini * yO_ini * math.exp(Q_PO_b + b_PO_b / T_ini) + P_ini * yA_ini * math.exp(Q_PA_b + b_PA_b / T_ini))
        qO_ini = P_ini * yO_ini * math.exp(Q_PO + b_PO / T_ini) / (1e5 + P_ini * yN_ini * math.exp(Q_PN_b + b_PN_b / T_ini) + P_ini * yO_ini * math.exp(Q_PO_b + b_PO_b / T_ini) + P_ini * yA_ini * math.exp(Q_PA_b + b_PA_b / T_ini))
        qA_ini = P_ini * yA_ini * math.exp(Q_PA + b_PA / T_ini) / (1e5 + P_ini * yN_ini * math.exp(Q_PN_b + b_PN_b / T_ini) + P_ini * yO_ini * math.exp(Q_PO_b + b_PO_b / T_ini) + P_ini * yA_ini * math.exp(Q_PA_b + b_PA_b / T_ini))
    
        u_ini = 0
        qN_Star_ini = qN_ini
        qO_Star_ini = qO_ini
        qA_Star_ini = qA_ini
        EB_ini = (1 - Porosity_bed) * (rho_s * Cp_s + Cp_a * (qN_ini + qO_ini + qA_ini)) + Cp_g * Porosity_bed * P_ini / R / T_ini
    
        x_ini = [yN_ini, yO_ini, P_ini, T_ini, qN_ini, qO_ini, qA_ini]
        x_ini = [item for item in x_ini for i in range(n_space)]
    
        z_ini_1 = [yN_feed, yO_feed, P_high, T_feed, yN_ini, yO_ini, P_ini, T_ini]
        z_ini_2 = [u_ini]
        z_ini_2 = [item for item in z_ini_2 for i in range(n_space+2)]        
        z_ini_3 = [yN_ini, yO_ini, P_ini, T_ini, u_ini]
        z_ini_3 = [item for item in z_ini_3 for i in range(n_space+1)]
        z_ini_4 = [yA_ini, qN_Star_ini, qO_Star_ini, qA_Star_ini, EB_ini]
        z_ini_4 = [item for item in z_ini_4 for i in range(n_space+2)]
        z_ini = z_ini_1 + z_ini_2 + z_ini_3 + z_ini_4
        
        # %% iterations
        for s in range(0, 300):
    
            x_ini = x_ini[0:7*n_space] + [0]
            try:
                [x_1, ode_1, z_1, alg_1] = pressurization(n_space, yN_feed, yO_feed, T_feed, P_high)
                dae_1 = {'x':x_1, 'ode':ode_1, 'z':z_1, 'alg':alg_1}
                F = integrator('F', 'idas', dae_1, {'t0':0, 'tf':delta_t1})
                F_sol_1 = F(x0 = x_ini, z0 = z_ini)
            except:
                try:
                    [x_1, ode_1, z_1, alg_1] = pressurization(n_space, yN_feed, yO_feed, T_feed, P_high + 100)
                    dae_1 = {'x':x_1, 'ode':ode_1, 'z':z_1, 'alg':alg_1}
                    F = integrator('F', 'idas', dae_1, {'t0':0, 'tf':delta_t1})
                    F_sol_1 = F(x0 = x_ini, z0 = z_ini)
                except:
                    [x_1, ode_1, z_1, alg_1] = pressurization(n_space, yN_feed, yO_feed, T_feed, P_high - 100)
                    dae_1 = {'x':x_1, 'ode':ode_1, 'z':z_1, 'alg':alg_1}
                    F = integrator('F', 'idas', dae_1, {'t0':0, 'tf':delta_t1})
                    F_sol_1 = F(x0 = x_ini, z0 = z_ini)
            xf_1 = F_sol_1['xf'].elements()
            zf_1 = F_sol_1['zf'].elements()
            PRE_Input = np.float(xf_1[7*n_space])
    
            xf_1 = xf_1[0:7*n_space] + [0, 0, 0, 0]
            try:               
                [x_2, ode_2, z_2, alg_2] = adsorption_one(n_space, yN_feed, yO_feed, T_feed, P_high, P_ads_end)
                dae_2 = {'x':x_2, 'ode':ode_2, 'z':z_2, 'alg':alg_2}
                F = integrator('F', 'idas', dae_2, {'t0':0, 'tf': delta_t2})
                F_sol_2 = F(x0 = xf_1, z0 = zf_1)
            except:
                try:
                    [x_2, ode_2, z_2, alg_2] = adsorption_one(n_space, yN_feed, yO_feed, T_feed, P_high, P_ads_end + 100)
                    dae_2 = {'x':x_2, 'ode':ode_2, 'z':z_2, 'alg':alg_2}
                    F = integrator('F', 'idas', dae_2, {'t0':0, 'tf': delta_t2})
                    F_sol_2 = F(x0 = xf_1, z0 = zf_1)                
                except:
                    [x_2, ode_2, z_2, alg_2] = adsorption_one(n_space, yN_feed, yO_feed, T_feed, P_high, P_ads_end - 100)
                    dae_2 = {'x':x_2, 'ode':ode_2, 'z':z_2, 'alg':alg_2}
                    F = integrator('F', 'idas', dae_2, {'t0':0, 'tf': delta_t2})
                    F_sol_2 = F(x0 = xf_1, z0 = zf_1)
            xf_2 = F_sol_2['xf'].elements()
            zf_2 = F_sol_2['zf'].elements()
            ADS_Input = np.float(xf_2[7*n_space])
            ADS_Output = np.float(xf_2[7*n_space + 1])
            ADS_PN_out = np.float(xf_2[7*n_space + 2])
            ADS_PO_out = np.float(xf_2[7*n_space + 3])
                
            xf_2 = xf_2[0:7*n_space] + [0, 0, 0, 0]
            try:
                [x_3, ode_3, z_3, alg_3] = adsorption_two(n_space, yN_feed, yO_feed, T_feed, P_high, P_adp_end)
                dae_3 = {'x':x_3, 'ode':ode_3, 'z':z_3, 'alg':alg_3}
                F = integrator('F', 'idas', dae_3, {'t0':0, 'tf': delta_t3})
                F_sol_3 = F(x0 = xf_2, z0 = zf_2)
            except:
                try:
                    [x_3, ode_3, z_3, alg_3] = adsorption_two(n_space, yN_feed, yO_feed, T_feed, P_high, P_adp_end + 100)
                    dae_3 = {'x':x_3, 'ode':ode_3, 'z':z_3, 'alg':alg_3}
                    F = integrator('F', 'idas', dae_3, {'t0':0, 'tf': delta_t3 + 0.001})
                    F_sol_3 = F(x0 = xf_2, z0 = zf_2)
                except:
                    [x_3, ode_3, z_3, alg_3] = adsorption_two(n_space, yN_feed, yO_feed, T_feed, P_high, P_adp_end - 100)
                    dae_3 = {'x':x_3, 'ode':ode_3, 'z':z_3, 'alg':alg_3}
                    F = integrator('F', 'idas', dae_3, {'t0':0, 'tf': delta_t3 - 0.001})
                    F_sol_3 = F(x0 = xf_2, z0 = zf_2)
            xf_3 = F_sol_3['xf'].elements()
            zf_3 = F_sol_3['zf'].elements()
            ADP_Input = np.float(xf_3[7*n_space])
            ADP_Output = np.float(xf_3[7*n_space + 1])           
            ADP_PN_out = np.float(xf_3[7*n_space + 2])
            ADP_PO_out = np.float(xf_3[7*n_space + 3])
            ADP_yN_out = ADP_PN_out / ADP_Output
            ADP_yO_out = ADP_PO_out / ADP_Output
    
            xf_3 = xf_3[0:7*n_space] + [0, 0, 0]
            try:
                [x_4, ode_4, z_4, alg_4] = equalization_one(n_space, P_equ_end)
                dae_4 = {'x':x_4, 'ode':ode_4, 'z':z_4, 'alg':alg_4}
                F = integrator('F', 'idas', dae_4, {'t0':0, 'tf': delta_t4})
                F_sol_4 = F(x0 = xf_3, z0 = zf_3)
            except:
                try:
                    [x_4, ode_4, z_4, alg_4] = equalization_one(n_space, P_equ_end + 100)
                    dae_4 = {'x':x_4, 'ode':ode_4, 'z':z_4, 'alg':alg_4}
                    F = integrator('F', 'idas', dae_4, {'t0':0, 'tf': delta_t4})
                    F_sol_4 = F(x0 = xf_3, z0 = zf_3)
                except:
                    [x_4, ode_4, z_4, alg_4] = equalization_one(n_space, P_equ_end - 100)
                    dae_4 = {'x':x_4, 'ode':ode_4, 'z':z_4, 'alg':alg_4}
                    F = integrator('F', 'idas', dae_4, {'t0':0, 'tf': delta_t4})
                    F_sol_4 = F(x0 = xf_3, z0 = zf_3)
            xf_4 = F_sol_4['xf'].elements()
            zf_4 = F_sol_4['zf'].elements()
            EQU_Output = np.abs(np.float(xf_4[7*n_space] + xf_4[7*n_space + 1] + xf_4[7*n_space + 2]))
            EQU_yN_out = np.float(xf_4[7*n_space]) / EQU_Output
            EQU_yO_out = np.float(xf_4[7*n_space + 1]) / EQU_Output

            xf_4 = xf_4[0:7*n_space] + [0, 0, 0]
            try:
                [x_5, ode_5, z_5, alg_5] = desorption_one(n_space, P_low)
                dae_5 = {'x':x_5, 'ode':ode_5, 'z':z_5, 'alg':alg_5}
                F = integrator('F', 'idas', dae_5, {'t0':0, 'tf': delta_t1 + delta_t2})
                F_sol_5 = F(x0 = xf_4, z0 = zf_4)
            except:
                try:
                    [x_5, ode_5, z_5, alg_5] = desorption_one(n_space, P_low + 100)
                    dae_5 = {'x':x_5, 'ode':ode_5, 'z':z_5, 'alg':alg_5}
                    F = integrator('F', 'idas', dae_5, {'t0':0, 'tf': delta_t1 + delta_t2})
                    F_sol_5 = F(x0 = xf_4, z0 = zf_4)
                except:
                    [x_5, ode_5, z_5, alg_5] = desorption_one(n_space, P_low - 100)
                    dae_5 = {'x':x_5, 'ode':ode_5, 'z':z_5, 'alg':alg_5}
                    F = integrator('F', 'idas', dae_5, {'t0':0, 'tf': delta_t1 + delta_t2})
                    F_sol_5 = F(x0 = xf_4, z0 = zf_4)   
            xf_5 = F_sol_5['xf'].elements()
            zf_5 = F_sol_5['zf'].elements()
            DES_Output = np.float(xf_5[7*n_space])
            DES_PN_out = np.float(xf_5[7*n_space + 1])
            DES_PO_out = np.float(xf_5[7*n_space + 2])
            
            xf_6 = xf_5
            zf_6 = zf_5
            
            xf_6 = xf_6[0:7*n_space] + [0, 0, 0, 0, 0, 0]
            try:
                [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_adp_end)
                dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                F_sol_7 = F(x0 = xf_6, z0 = zf_6)
            except:
                try:
                    [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_adp_end + 100)
                    dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                    F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                    F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                except:
                    [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_adp_end - 100)
                    dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                    F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                    F_sol_7 = F(x0 = xf_6, z0 = zf_6)                
            xf_7 = F_sol_7['xf'].elements()
            PUR_Input = np.abs(np.float(xf_7[7*n_space] + xf_7[7*n_space + 1] + xf_7[7*n_space + 2]))
            if PUR_Input < ADP_Output:
                raise
            else:
                for i in range(1, 10):
                    P_pur_feed = int(P_adp_end - 0.1 * P_bar * i)
                    try:
                        [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed)
                        dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                        F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                        F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                    except:
                        try:
                            [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed + 100)
                            dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                            F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                            F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                        except:
                            [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed - 100)
                            dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                            F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                            F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                    xf_7 = F_sol_7['xf'].elements()
                    PUR_Input = np.abs(np.float(xf_7[7*n_space] + xf_7[7*n_space + 1] + xf_7[7*n_space + 2]))
                    if PUR_Input <= ADP_Output:
                        for j in range(1, 21):
                            P_pur_feed_final = int(P_pur_feed + (0.1 - 0.005 * j) * P_bar)
                            try:
                                [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed_final)
                                dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                                F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                                F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                            except:
                                try:
                                    [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed_final + 100)
                                    dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                                    F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                                    F_sol_7 = F(x0 = xf_6, z0 = zf_6)                            
                                except:
                                    [x_7, ode_7, z_7, alg_7] = purge(n_space, ADP_yN_out, ADP_yO_out, T_feed, P_low, P_pur_feed_final - 100)
                                    dae_7 = {'x':x_7, 'ode':ode_7, 'z':z_7, 'alg':alg_7}
                                    F = integrator('F', 'idas', dae_7, {'t0':0, 'tf': delta_t3})
                                    F_sol_7 = F(x0 = xf_6, z0 = zf_6)
                            xf_7 = F_sol_7['xf'].elements()
                            PUR_Input = np.abs(np.float(xf_7[7*n_space] + xf_7[7*n_space + 1] + xf_7[7*n_space + 2]))
                            if PUR_Input <= ADP_Output:
                                break
                        break
            zf_7 = F_sol_7['zf'].elements()
            PUR_PO_feed = np.float(xf_7[7*n_space + 1])
            PUR_Output = np.float(xf_7[7*n_space + 3])
            PUR_PN_out = np.float(xf_7[7*n_space + 4])
            PUR_PO_out = np.float(xf_7[7*n_space + 5])

            xf_7 = xf_7[0:7*n_space] + [0, 0, 0]
            try:
                [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_equ_end)    
                dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                F_sol_8 = F(x0 = xf_7, z0 = zf_7)
            except:
                try:
                    [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_equ_end + 100)    
                    dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                    F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                    F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                except:
                    [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_equ_end - 100)    
                    dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                    F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                    F_sol_8 = F(x0 = xf_7, z0 = zf_7)
            xf_8 = F_sol_8['xf'].elements()
            EQN_Input = np.float(xf_8[7*n_space] + xf_8[7*n_space + 1] + xf_8[7*n_space + 2])
            if EQN_Input < EQU_Output:
                raise
            else:
                for i in range(1, 10):
                    P_eqn_feed = int(P_equ_end - 0.1 * P_bar * i)
                    try:
                        [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed)
                        dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                        F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                        F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                    except:
                        try:
                            [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed + 100)
                            dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                            F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                            F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                        except:
                            [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed - 100)
                            dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                            F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                            F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                    xf_8 = F_sol_8['xf'].elements()
                    EQN_Input = np.float(xf_8[7*n_space] + xf_8[7*n_space + 1] + xf_8[7*n_space + 2])
                    if EQN_Input <= EQU_Output:
                        for j in range(1, 21):
                            P_eqn_feed_final = int(P_eqn_feed + (0.1 - 0.005 * j) * P_bar)
                            try:
                                [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed_final)
                                dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}
                                F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                                F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                            except:
                                try:
                                    [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed_final + 100)
                                    dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}                                    
                                    F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                                    F_sol_8 = F(x0 = xf_7, z0 = zf_7)
                                except:
                                    [x_8, ode_8, z_8, alg_8] = equalization_two(n_space, EQU_yN_out, EQU_yO_out, T_feed, P_eqn_feed_final - 100)
                                    dae_8 = {'x':x_8, 'ode':ode_8, 'z':z_8, 'alg':alg_8}                                    
                                    F = integrator('F', 'idas', dae_8, {'t0':0, 'tf': delta_t4})
                                    F_sol_8 = F(x0 = xf_7, z0 = zf_7)    
                            xf_8 = F_sol_8['xf'].elements()
                            EQN_Input = np.float(xf_8[7*n_space] + xf_8[7*n_space + 1] + xf_8[7*n_space + 2])
                            if EQN_Input <= EQU_Output:
                                break
                        break
            zf_8 = F_sol_8['zf'].elements()

            xstart = np.asarray(x_ini[0:6*n_space])
            xend = np.asarray(xf_8[0:6*n_space])
            error = np.max(np.abs((xstart - xend) / (xstart + 1e-10)))
            abs_error = np.max(np.abs(xstart - xend))
            
            print("%d th iter, max error: %f, %f" %(s, error, abs_error))
            if error <= 0.01:
                break
            else:
                x_ini = xf_8
                z_ini = zf_8
                z_ini[0] = yN_feed
                z_ini[1] = yO_feed
                z_ini[2] = P_high
                z_ini[3] = T_feed

# %% performance

        # PurityO = ADS_PO_out / ADS_Output
        # ProductivityO = 0.5 * (ADS_PO_out + ADP_PO_out - PUR_PO_feed) / (delta_t1 + delta_t2 + delta_t3 + delta_t4)         # mol/s
        # Volume = 3.1416 * R_inner * R_inner * L * (1 - Porosity_bed) * rho_s                                                # 7758.5 kg adsorbents
        # ProductivityO = ProductivityO * 80640 / Volume                                                                      # 1 mol/s/kg = 80640 Nm3/h/ton
        
        # P_pur_feed_final_result.append(P_pur_feed_final)
        # P_eqn_feed_final_result.append(P_eqn_feed_final)
        # PO_purity_result.append(PurityO)        
        # PO_productivity_result.append(ProductivityO)

        # PRE_Input_result.append(PRE_Input)
        # ADS_Input_result.append(ADS_Input)
        # ADS_Output_result.append(ADS_Output)
        # ADS_PO_out_result.append(ADS_PO_out)
        # ADP_Input_result.append(ADP_Input)
        # ADP_Output_result.append(ADP_Output)
        # EQU_Output_result.append(EQU_Output)
        # DES_Output_result.append(DES_Output)
        # PUR_Input_result.append(PUR_Input)
        # PUR_Output_result.append(PUR_Output)
        # EQN_Input_result.append(EQN_Input)
    # except:
    #     P_pur_feed_final_result.append(np.nan)
    #     P_eqn_feed_final_result.append(np.nan)
    #     PO_purity_result.append(np.nan)
    #     PO_productivity_result.append(np.nan)
        
    #     PRE_Input_result.append(np.nan)
    #     ADS_Input_result.append(np.nan)
    #     ADS_Output_result.append(np.nan)
    #     ADS_PO_out_result.append(np.nan)
    #     ADP_Input_result.append(np.nan)
    #     ADP_Output_result.append(np.nan)
    #     EQU_Output_result.append(np.nan)
    #     DES_Output_result.append(np.nan)
    #     PUR_Input_result.append(np.nan)
    #     PUR_Output_result.append(np.nan)
    #     EQN_Input_result.append(np.nan)

# %% save

# P_pur_feed_final_result = np.asarray(P_pur_feed_final_result) / P_bar
# P_eqn_feed_final_result = np.asarray(P_eqn_feed_final_result) / P_bar

# P_pur_feed_final_result = np.array(P_pur_feed_final_result).reshape(-1, 1)
# P_eqn_feed_final_result = np.array(P_eqn_feed_final_result).reshape(-1, 1)
# PO_purity_result = np.array(PO_purity_result).reshape(-1, 1)
# PO_productivity_result = np.array(PO_productivity_result).reshape(-1, 1)

# PRE_Input_result = np.array(PRE_Input_result).reshape(-1, 1)
# ADS_Input_result = np.array(ADS_Input_result).reshape(-1, 1)
# ADS_Output_result = np.array(ADS_Output_result).reshape(-1, 1)
# ADP_Input_result = np.array(ADP_Input_result).reshape(-1, 1)
# ADP_Output_result = np.array(ADP_Output_result).reshape(-1, 1)
# EQU_Output_result = np.array(EQU_Output_result).reshape(-1, 1)
# DES_Output_result = np.array(DES_Output_result).reshape(-1, 1)
# PUR_Input_result = np.array(PUR_Input_result).reshape(-1, 1)
# PUR_Output_result = np.array(PUR_Output_result).reshape(-1, 1)
# EQN_Input_result = np.array(EQN_Input_result).reshape(-1, 1)
# Balance = PRE_Input_result + ADS_Input_result - ADS_Output_result + ADP_Input_result - DES_Output_result - PUR_Output_result
# O_balance = (PRE_Input_result + ADS_Input_result + ADP_Input_result) * yO_feed - ADS_PO_out - DES_PO_out - PUR_PO_out

# df = df[z_min: z_max, :]
# Result = np.hstack([df, 
#                     PO_purity_result, PO_productivity_result, P_pur_feed_final_result, P_eqn_feed_final_result,
#                     PRE_Input_result, ADS_Input_result, ADS_Output_result, ADP_Input_result, 
#                     EQU_Output_result, DES_Output_result, PUR_Input_result, PUR_Output_result, EQN_Input_result, Balance, O_balance])

# Header = np.array(['P_high', 'P_ads', 'P_adp', 'P_equ', 'P_low', 'delta_t1', 'delta_t2', 'delta_t3', 'delta_t4',
#                     'PO_purity_result', 'PO_productivity_result', 'P_pur_feed', 'P_eqn_feed',
#                     'PRE_Input', 'ADS_Input', 'ADS_Output', 'ADP_Input', 
#                     'EQU_Output', 'DES_Output', 'PUR_Input', 'PUR_Output', 'EQN_Input', 'Balance', 'O_balance'])
# Result = np.vstack([Header, Result])

# # %% output to excel

# Result = pd.DataFrame(Result)

# save_name = int(rank)
# save_name = 'Result_' + str(save_name) + '.xlsx'
# Result.to_excel(save_name, header = False, index = False)  