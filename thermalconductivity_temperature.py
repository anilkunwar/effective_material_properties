import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def thermal_conductivity(T, A, B, C, D, E, F, Tref, Tmin, Tmax):
    kth = A * (T - Tref) + B * (T - Tref)**2 + C * (T - Tref)**3 + D * (T - Tref)**4 + E * (T - Tref)**5 + F
    return kth if Tmin <= T <= Tmax else np.nan

def main():
    st.title("Thermal Conductivity vs Temperature")
    st.latex(r"k_{\phi} = A_{\phi}  \cdot (T - T_{\text{ref}}) + B_{\phi}  \cdot (T - T_{\text{ref}})^2 + C_{\phi}  \cdot (T - T_{\text{ref}})^3 + D_{\phi}  \cdot (T - T_{\text{ref}})^4 + E_{\phi}  \cdot (T - T_{\text{ref}})^5 + F_{\phi}")

    st.sidebar.title("Input Parameters")

    num_phases = st.sidebar.number_input("Number of Material Phases", min_value=1, step=1, value=1)

    phases = []
    for i in range(num_phases):
        phase_name = st.sidebar.text_input(f"Material Phase {i + 1} Name")
        key_prefix = f"{phase_name}_{i}"  # Unique key prefix based on phase name and index
        A = st.sidebar.number_input(f"A_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_A")
        B = st.sidebar.number_input(f"B_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_B")
        C = st.sidebar.number_input(f"C_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_C")
        D = st.sidebar.number_input(f"D_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_D")
        E = st.sidebar.number_input(f"E_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_E")
        F = st.sidebar.number_input(f"F_{phase_name}", value=0.0, format="%.8e", key=f"{key_prefix}_F")
        Tref = st.sidebar.number_input(f"Tref_{phase_name}", value=0.0, key=f"{key_prefix}_Tref")
        Tmin_phase, Tmax_phase = st.sidebar.slider(f"Temperature Range for {phase_name} (K)", min_value=100, max_value=2500, value=(300, 2500), key=f"{key_prefix}_slider")
        color = st.sidebar.color_picker(f"Color for {phase_name}", "#000000", key=f"{key_prefix}_color")
        linewidth = st.sidebar.slider(f"Linewidth for {phase_name}", min_value=1, max_value=10, value=2, key=f"{key_prefix}_linewidth")
        phases.append((phase_name, A, B, C, D, E, F, Tref, Tmin_phase, Tmax_phase, color, linewidth))

    T_values = np.linspace(100, 2500, 1000)

    fig, ax = plt.subplots(figsize=(8, 6))

    for phase in phases:
        phase_name, A, B, C, D, E, F, Tref,  Tmin_phase, Tmax_phase, color, linewidth = phase
        kth_values = [thermal_conductivity(T, A, B, C, D, E, F, Tref, Tmin_phase, Tmax_phase) for T in T_values]
        ax.plot(T_values, kth_values, label=phase_name, color=color, linewidth=linewidth)

        st.text("")
        st.text(f"Fitted Equation of Thermal Conductivity for {phase_name}:")
        st.latex(rf"k_{{{phase_name}}} = {A} \cdot (T - {Tref}) + {B} \cdot (T - {Tref})^2 + {C} \cdot (T - {Tref})^3 + {D} \cdot (T - {Tref})^4 + {E} \cdot (T - {Tref})^5 + {F}")

    ax.set_xlabel('Temperature (K)', fontsize=20) #14 
    ax.set_ylabel('Thermal Conductivity (W/m K)', fontsize=20)
    legend_loc = st.sidebar.radio(f"Legend Location", ("best", "upper right", "upper left", "lower left", "lower right", "right", "center left", "center right", "lower center", "upper center", "center"), index=0) 
    ax.legend(fontsize='large', loc=legend_loc)
    ax.set_yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=15)  # Adjust tick font size
    ax.xaxis.set_tick_params(width=5)  # Adjust x-axis thickness
    ax.yaxis.set_tick_params(width=5)  # Adjust y-axis thickness
    ax.spines['bottom'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['left'].set_linewidth(2)  # Adjust left axis thickness
    ax.spines['top'].set_linewidth(2)  # Adjust bottom axis thickness
    ax.spines['right'].set_linewidth(2)  # Adjust left axis thickness

    st.pyplot(fig)

if __name__ == "__main__":
    main()

