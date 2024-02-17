import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, AutoMinorLocator, FixedLocator, MultipleLocator,  MaxNLocator


def kth_cu_solid(T):
    Ascu = 0.00068
    Bscu = -0.30587
    Dscu = 452.1324
    return Ascu*T**3 + Bscu*T**2 + Dscu*T if T < 1300 else np.nan

def kth_cu_liquid(T):
    Alcu = 0.00156
    Blcu = -1.252
    Dlcu = 0.00156
    return Alcu*T**3 + Blcu*T**2 + Dlcu*T if T > 1390 else np.nan

#def kth_air_gas(T):
#    Agas = 1.2E-4 #0.00012
#    Bgas = -1.557E-4 #-0.001557
#    return Agas*(T-250)**2 + Bgas*(T-250) + 2.623E-02
    
# Reference
# https://www.cambridge.org/core/books/abs/gas-turbines/equations-of-air-thermophysical-properties/9572106E068EFF1B7C0896124C17A196
def kth_air_gas(T):
    Agas = 1.708186E-4 #0.00012
    Bgas = 0 #-2.3758E-5 #-0.001557
    Cgas = 0 #2.2E-10 
    Dgas = 0 #-9.45E-14 
    Egas = 0 # 1.58E-18 
    Fgas = -7.488E-3  #-7.488E-3
    return Agas*(T)  + Bgas*T**2 + Cgas*T**3 + Dgas*T**4 + Egas*T**5 + Fgas

def kth_tialv_solid(T):
    As = 1.46E-02
    Cs = -0.32
    return As*T + Cs

def kth_tialv_liquid(T):
    if 1950.0 < T < 2800.0:
        Al = 1.83E-02
        Cl = -6.66
        kth = Al*T + Cl
    else:
        kth = np.nan
    return kth
    
#def kth_ti6alv_liquid(T):
    #Al = 1.83E-02 #1.83E-02
    #Cl = -6.66
    #return Al*(T-0.0) + Cl if 298 < T < 2800 else np.nan

def main():
    st.title("Thermal Conductivity as Function of Temperature ")

    temperature_range = st.slider("Temperature Range (K)", min_value=250, max_value=3000, value=(250, 2500))

    T_values = np.linspace(temperature_range[0], temperature_range[1], 1000)
    kth_cu_solid_values = [kth_cu_solid(T) for T in T_values if T < 1300]
    kth_cu_liquid_values = [kth_cu_liquid(T) for T in T_values if T > 1400]
    kth_air_gas_values = [kth_air_gas(T) for T in T_values if 250 < T < 1358]
    kth_tialv_solid_values = [kth_tialv_solid(T) for T in T_values if 298 < T < 1800]
    kth_tialv_liquid_values = [kth_tialv_liquid(T) for T in T_values if 298 < T < 2700]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(T_values[:len(kth_cu_solid_values)], kth_cu_solid_values, label='Cu (s)', color='blue' , linewidth=5)
    ax.plot(T_values[len(kth_cu_solid_values):len(kth_cu_solid_values)+len(kth_cu_liquid_values)], kth_cu_liquid_values, label='Cu (l)', color='red', linewidth=5)
    ax.plot(T_values[:len(kth_air_gas_values)], kth_air_gas_values, label='Air', color='magenta', linewidth=5)
    ax.plot(T_values[:len(kth_tialv_solid_values)], kth_tialv_solid_values, label='Ti6Al4V (s)', color='green', linewidth=5)
    ax.plot(T_values[:len(kth_tialv_liquid_values)], kth_tialv_liquid_values, label='Ti6Al4V (l)', color='purple', linewidth=5)
    ax.set_xlabel('Temperature (K)', fontsize=22)
    ax.set_ylabel('k$_{th}$ (W/mK)', fontsize=22)
    #ax.set_title('Thermal Conductivity vs Temperature')
    #ax.legend()
    ax.legend(fontsize='large', loc='center left')  # Adjust fontsize and loc here
    ax.set_yscale('log')
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
    #ax.yaxis.set_minor_formatter(plt.NullFormatter())
    #minor_tick_positions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    #ax.yaxis.set_minor_locator(FixedLocator(minor_tick_positions))
    #ax.yaxis.set_minor_locator(MultipleLocator(1))
    #minor_ticks = np.logspace(np.log10(ax.get_ylim()[0]), np.log10(ax.get_ylim()[1]), num=9, base=10.0)[1:-1]
    #ax.yaxis.set_minor_locator(FixedLocator(minor_ticks))
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(0.1, 1, 0.1)))
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
    #ax.yaxis.set_minor_formatter(plt.NullFormatter())
    # Add subticks manually
    # Add subticks manually
    minor_ticks = np.logspace(np.log10(ax.get_ylim()[0]), np.log10(ax.get_ylim()[1]), num=9, base=10.0)
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=minor_ticks))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())

    
    ax.tick_params(axis='both', which='major', labelsize=19)  # Adjust tick font size
    ax.xaxis.set_tick_params(width=5)  # Adjust x-axis thickness
    ax.yaxis.set_tick_params(width=5)  # Adjust y-axis thickness
    ax.spines['bottom'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['left'].set_linewidth(2)  # Adjust left axis thickness
    ax.spines['top'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['right'].set_linewidth(2)  # Adjust left axis thickness

    st.pyplot(fig)

if __name__ == "__main__":
    main()

