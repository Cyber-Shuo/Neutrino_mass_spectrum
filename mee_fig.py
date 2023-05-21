import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import *
from tkinter import filedialog
import time
from pathlib import Path
import scipy.special as sps
from itertools import chain

M1 = np.arange(0.0001, 0.001, 0.0001)
M2 = np.arange(0.001, 0.01, 0.0001)
M3 = np.arange(0.01, 0.1, 0.001)
M4 = np.arange(0.1, 1, 0.01)
M = np.concatenate((M1,M2,M3,M4),axis=0)
MIN_m_ee_NO, MAX_m_ee_NO, MAX_m_ee_IO, MIN_m_ee_IO = [],[],[],[]
# 定义1sigma内质量平方差ev^2
delta_m_21 = 7.41*(10**(-5))
delta_m_31_NO = 2.507*(10**(-3))
delta_m_32_IO = -2.486*(10**(-3))

# 定义sin^2(\theta_12)
sin_theta_12_2_NO = 0.303
sin_theta_12_2_IO = 0.303

# 定义sin^2(\theta_13)
sin_theta_13_2_NO = 0.02225
sin_theta_13_2_IO = 0.02223

# 定义sin^2(\theta_23)
sin_theta_23_2_NO = 0.451
sin_theta_23_2_IO = 0.569

# 定义马拉约拿相位
def majorana_phases(t):
    m_phases = np.exp(t*2j)
    return m_phases

# 获得 正序的有效质量最小值 关于 最小中微子质量m_1 的函数
def get_min_m_ee_NO(m):
    m_ee_1_NO = m*(1-sin_theta_12_2_NO)*(1-sin_theta_13_2_NO)
    m_ee_2_NO = np.sqrt(m**2 + delta_m_21)*sin_theta_12_2_NO*(1-sin_theta_13_2_NO)
    m_ee_3_NO = np.sqrt(m**2 + delta_m_31_NO)*sin_theta_13_2_NO
    min_m_ee_NO = np.abs(- m_ee_1_NO + m_ee_2_NO - m_ee_3_NO)
    return min_m_ee_NO

# 获得 正序的有效质量最大值 关于 最小中微子质量m_1 的函数
def get_max_m_ee_NO(m):
    m_ee_1_NO = m*(1-sin_theta_12_2_NO)*(1-sin_theta_13_2_NO)
    m_ee_2_NO = np.sqrt(m**2 + delta_m_21)*sin_theta_12_2_NO*(1-sin_theta_13_2_NO)
    m_ee_3_NO = np.sqrt(m**2 + delta_m_31_NO)*sin_theta_13_2_NO
    max_m_ee_NO = m_ee_1_NO + m_ee_2_NO + m_ee_3_NO
    return max_m_ee_NO

# 获得 反序的有效质量最大值 关于 最小中微子质量m_3 的函数
def get_max_m_ee_IO(m):
    m_ee_1_IO = np.sqrt(m**2 - delta_m_32_IO - delta_m_21)*(1-sin_theta_12_2_IO)*(1-sin_theta_13_2_IO)
    m_ee_2_IO = np.sqrt(m**2 - delta_m_32_IO)*sin_theta_12_2_IO*(1-sin_theta_13_2_IO)
    m_ee_3_IO = m*sin_theta_13_2_IO
    max_m_ee_IO = m_ee_1_IO + m_ee_2_IO + m_ee_3_IO
    return max_m_ee_IO

# 获得 反序的有效质量最小值 关于 最小中微子质量m_3 的函数
def get_min_m_ee_IO(m):
    m_ee_1_IO = np.sqrt(m**2 - delta_m_32_IO - delta_m_21)*(1-sin_theta_12_2_IO)*(1-sin_theta_13_2_IO)
    m_ee_2_IO = np.sqrt(m**2 - delta_m_32_IO)*sin_theta_12_2_IO*(1-sin_theta_13_2_IO)
    m_ee_3_IO = m*sin_theta_13_2_IO
    max_m_ee_IO = m_ee_1_IO - m_ee_2_IO - m_ee_3_IO
    return max_m_ee_IO

for M_v in M:
    min_m_ee_NO_v = get_min_m_ee_NO(M_v)
    MIN_m_ee_NO.append(min_m_ee_NO_v)

    max_m_ee_NO_v = get_max_m_ee_NO(M_v)
    MAX_m_ee_NO.append(max_m_ee_NO_v)

    max_m_ee_IO_v = get_max_m_ee_IO(M_v)
    MAX_m_ee_IO.append(max_m_ee_IO_v)

    min_m_ee_IO_v = get_min_m_ee_IO(M_v)
    MIN_m_ee_IO.append(min_m_ee_IO_v)

print(MIN_m_ee_NO)
fig, ax = plt.subplots()

ax.plot(M, MIN_m_ee_NO, color = 'r', label = 'NO min')
ax.plot(M, MAX_m_ee_NO, color = 'g', label = 'NO max')
ax.plot(M, MAX_m_ee_IO, color = 'b', label = 'IO max')
ax.plot(M, MIN_m_ee_IO, color = 'y', label = 'IO min')

ax.set_xscale('log')
ax.set_xlabel('minimal neutrino mass/eV')
ax.set_yscale('log')
ax.set_ylabel(r'$|m_{ee}|$'+'/eV')
plt.fill_between(M, MIN_m_ee_IO, MAX_m_ee_IO, color = 'b', alpha=0.5)
plt.fill_between(M, MIN_m_ee_NO, MAX_m_ee_NO, color = 'orange', alpha=0.5)

plt.legend()
plt.show()