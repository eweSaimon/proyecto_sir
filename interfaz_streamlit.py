# -*- coding: utf-8 -*-
"""
interfaz_streamlit.py
Interfaz web sencilla para demostrar el modelo SIR de propagación de malware
SIN modificar el código. Reutiliza las funciones de programa_malware_sir.py.

Ejecución:
    pip install streamlit
    streamlit run interfaz_streamlit.py

Todos los valores son PARÁMETROS HIPOTÉTICOS de simulación (no datos reales).
"""

import math

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import programa_malware_sir as m

st.set_page_config(page_title="Propagación de malware: modelo SIR", layout="wide")
st.title("Propagación de malware en una red: modelo SIR con intervención")
st.caption("Simulación matemática con parámetros hipotéticos. No se ejecuta ni se manipula software malicioso real.")

# ------------------------- Barra lateral: parámetros -------------------------
with st.sidebar:
    st.header("Parámetros")
    # Límites fijos y claves estables: así cambiar N no reinicia los demás controles.
    # La coherencia entre I(0) y N (S(0) > 0) la comprueba la validación del programa.
    N = st.number_input("Cantidad de computadores N", min_value=2, max_value=1_000_000, value=1000, step=10, key="N")
    I0 = st.number_input("Infectados iniciales I(0)", min_value=1, max_value=1_000_000, value=10, step=1, key="I0")
    beta = st.number_input("Tasa de transmisión β [1/día]", min_value=0.0, max_value=5.0, value=0.30, step=0.01, format="%.3f")
    gamma = st.number_input("Tasa de recuperación γ [1/día]", min_value=0.001, max_value=5.0, value=0.10, step=0.01, format="%.3f")
    tmax = st.number_input("Tiempo de simulación [días]", min_value=1, max_value=1000, value=100, step=10)
    u_pct = st.slider("Intensidad de la intervención u [%]", min_value=0, max_value=100, value=50, step=1)
    comparar = st.checkbox("Incluir los escenarios del proyecto (u = 0 %, 50 %, 80 %)", value=True)

u = u_pct / 100.0
p = m.Parametros(N=float(N), I0=float(I0), beta=float(beta), gamma=float(gamma), tmax=float(tmax))

# Escenarios a simular
escenarios = [("Escenario del usuario (u = %d %%)" % u_pct, u)]
if comparar:
    escenarios = list(m.ESCENARIOS_POR_DEFECTO)
    if not any(math.isclose(u, v) for _, v in escenarios):
        escenarios.append((f"Escenario del usuario (u = {u_pct} %)", u))
else:
    escenarios.insert(0, ("Sin intervención (referencia)", 0.0))
    if math.isclose(u, 0.0):
        escenarios = escenarios[:1]

# ------------------------- Simulación con manejo de errores -------------------------
try:
    resultados = m.ejecutar_escenarios(p, escenarios)
except ValueError as exc:
    st.error(f"Parámetros no válidos:\n\n{exc}")
    st.stop()
except RuntimeError as exc:
    st.error(f"Error numérico: {exc}")
    st.stop()

# ------------------------- Indicadores principales -------------------------
u_c = m.umbral_intervencion(p.beta, p.gamma)
R0_base = p.beta / p.gamma
R0_u = m.numero_reproductivo(p.beta, p.gamma, u)

c1, c2, c3, c4 = st.columns(4)
c1.metric("R0 sin intervención (β/γ)", f"{R0_base:.3f}")
c2.metric("Umbral u* = 1 − γ/β", f"{u_c * 100:.2f} %")
c3.metric("R0(u) del usuario", f"{R0_u:.3f}")
c4.metric("Condición", "u > u*: controlado" if R0_u < 1 else ("u = u*: umbral" if math.isclose(R0_u, 1.0) else "u < u*: hay brote"))

if R0_u < 1:
    st.success(f"Con u = {u_pct} % se cumple R0(u) = {R0_u:.3f} < 1: la infección decrece y no hay brote epidémico.")
else:
    st.warning(f"Con u = {u_pct} % se tiene R0(u) = {R0_u:.3f} ≥ 1: hay crecimiento inicial de la infección.")

# ------------------------- Pestañas -------------------------
tab_I, tab_srI, tab_fase, tab_tabla, tab_verif, tab_sens = st.tabs(
    ["Infectados", "S, I y R por escenario", "Plano de fase", "Tabla y descarga", "Verificaciones", "Sensibilidad"])

with tab_I:
    st.pyplot(m.fig_compartimento(resultados, "I"))
    plt.close("all")
with tab_srI:
    st.pyplot(m.fig_comparacion(resultados, p))
    plt.close("all")
with tab_fase:
    st.pyplot(m.fig_plano_fase(resultados, p))
    plt.close("all")
with tab_tabla:
    resumen = m.construir_resumen(resultados, p)
    cols = ["escenario", "u", "R0_efectivo", "pico_I_numerico", "dia_pico", "S_en_tmax", "I_en_tmax",
            "R_en_tmax", "S_inf_teorico", "tasa_ataque_asintotica_pct"]
    st.dataframe(resumen[cols].round(3), width="stretch")
    serie = m.construir_serie_diaria(resultados, p)
    st.download_button("Descargar serie diaria (CSV)", serie.to_csv(index=False).encode("utf-8-sig"),
                       file_name="resultados_simulacion.csv", mime="text/csv")
    st.caption("S(t_max) puede diferir de S(∞) teórico si el brote aún no terminó dentro del horizonte simulado.")
with tab_verif:
    verif = m.verificar_resultados(resultados, p)
    st.dataframe(verif, width="stretch")
    fallas = int((verif["resultado"] == "FALLA").sum())
    (st.success if fallas == 0 else st.error)(f"{int((verif['resultado'] == 'OK').sum())} verificaciones OK, {fallas} fallas.")
    st.pyplot(m.fig_conservacion(resultados, p))
    plt.close("all")
with tab_sens:
    st.pyplot(m.fig_sensibilidad(p))
    plt.close("all")
    st.dataframe(m.tabla_sensibilidad(p).round(3), width="stretch")
