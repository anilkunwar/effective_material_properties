import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,  MaxNLocator

def rho_cu_liquid(T):
    Alcu = -2.383E-11
    Blcu = 8.69E-8
    Clcu = -7.96E-04
    Dlcu = 7.89
    rho = Alcu*(T-1330)**3 + Blcu*(T-1330)**2 + Dlcu*(T-1330) + Clcu
    return rho if T > 1358 else np.nan

def rho_cu_solid(T):
    if 300 < T < 1358:
        Ascu = -2.77E-11
        Bscu = -6.12E-8
        Cscu = -4.28E-04
        Dscu = 881.0
        rho = Ascu*(T-300)**3 + Bscu*(T-300)**2 + Dscu*(T-300) + Cscu
    else:
        rho = np.nan
    return rho

def rho_air_gas(T):
    Agas = -7.65e-07
    Bgas = 3.16e-03
    Cgas = 1.14e+01
    Dgas = -1.94
    Egas = 1.0
    rho = Agas*T**2 + Bgas*T + Dgas*np.log(Egas*T) + Cgas
    #return rho if T > 298 else np.nan
    return rho if 298 < T < 1400 else np.nan

def rho_tialv_solid(T):
    if 298 < T < 1923:
        As = 1.152
        Bs = -4.1197e-04
        Cs = -4.254e+02
        Ds = 6.54263e+03
        rho = As*T + Bs*T**2 + Ds*np.log(T) + Cs
    else:
        rho = np.nan
    return rho

def rho_tialv_liquid(T):
    if 1923.0 < T < 2500.0:
        Al = -0.452
        Bl = 4955.0
        rho = Al*T + Bl
    else:
        rho = np.nan
    return rho

def main():
    st.title("Density vs Temperature")

    temperature_range = st.slider("Temperature Range (K)", min_value=100, max_value=2500, value=(300, 2500))

    T_values = np.linspace(temperature_range[0], temperature_range[1], 1000)
    rho_cu_liquid_values = [rho_cu_liquid(T) for T in T_values]
    rho_cu_solid_values = [rho_cu_solid(T) for T in T_values]
    rho_air_gas_values = [rho_air_gas(T) for T in T_values]
    rho_tialv_solid_values = [rho_tialv_solid(T) for T in T_values]
    rho_tialv_liquid_values = [rho_tialv_liquid(T) for T in T_values]

    #fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(T_values, rho_cu_liquid_values, label='Cu (l)', color='red', linewidth=5)
    ax.plot(T_values, rho_cu_solid_values, label='Cu (s)', color='blue', linewidth=5)
    ax.plot(T_values, rho_air_gas_values, label='Air', color='magenta', linewidth=5)
    ax.plot(T_values, rho_tialv_solid_values, label='Ti6Al4V (s)', color='green', linewidth=5)
    ax.plot(T_values, rho_tialv_liquid_values, label='Ti6Al4V (l)', color='purple', linewidth=5)
    ax.set_xlabel('Temperature (K)', fontsize=20)
    ax.set_ylabel('Density (kg/m$^3$)', fontsize=20)
    #ax.set_title('Density vs Temperature')
    ax.legend(fontsize='large', loc='lower right')  # Adjust fontsize and loc here
    ax.set_yscale('log')
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(0.1, 1, 0.1)))
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
    #ax.yaxis.set_minor_formatter(plt.NullFormatter())
    # Add subticks manually
    # Add subticks manually
    #minor_ticks = np.logspace(np.log10(ax.get_ylim()[0]), np.log10(ax.get_ylim()[1]), num=9, base=10.0)
    #ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=minor_ticks))
    #ax.yaxis.set_minor_formatter(plt.NullFormatter())
    

    
    ax.tick_params(axis='both', which='major', labelsize=15)  # Adjust tick font size
    ax.xaxis.set_tick_params(width=5)  # Adjust x-axis thickness
    ax.yaxis.set_tick_params(width=5)  # Adjust y-axis thickness
    ax.spines['bottom'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['left'].set_linewidth(2)  # Adjust left axis thickness
    ax.spines['top'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['right'].set_linewidth(2)  # Adjust left axis thickness
    
    # Custom formatter to show "10^0" label on y-axis ticks
    #ax.yaxis.set_major_formatter(ScalarFormatter())
    
    # Set the number of major ticks on the y-axis
    #ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
    # Allow user to set y-axis limits
    #y_min, y_max = st.slider("Y-axis Limits", min_value=0, max_value=100000, value=(0, 100000))
    # Set y-axis limits based on user input
    #ax.set_ylim(y_min, y_max)


    st.pyplot(fig)

if __name__ == "__main__":
    main()

