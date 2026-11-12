# Propagación de malware en redes: modelo SIR con intervención

**Proyecto de aula — Ecuaciones Diferenciales, Ingeniería Informática**

> **Aviso.** Todos los valores numéricos (N = 1000, I(0) = 10, β = 0.30, γ = 0.10, 0–100 días)
> son **parámetros hipotéticos de simulación**. No provienen de estadísticas de ninguna empresa
> u organización real. El programa es una simulación matemática: **no ejecuta, genera ni
> manipula software malicioso**.

## Modelo

```
dS/dt = -(1-u) β S I / N
dI/dt =  (1-u) β S I / N - γ I
dR/dt =  γ I
```

S: susceptibles · I: infectados · R: recuperados/aislados · N = S + I + R (constante)
β: transmisión [1/día] · γ: recuperación o aislamiento [1/día] · u ∈ [0,1]: reducción relativa de la
transmisión por medidas preventivas y de contención (efecto *agregado*, no la eficacia de una herramienta).

R₀(u) = (1-u)β/γ · umbral u\* = 1 − γ/β (con los parámetros base, u\* = 2/3 ≈ 66.67 %).

## Contenido de la carpeta

| Archivo | Descripción |
|---|---|
| `programa_malware_sir.py` | Programa principal (validación, `solve_ivp` RK45, escenarios, gráficas, CSV, verificaciones). |
| `interfaz_streamlit.py` | Interfaz web sencilla para la demostración (sin editar código). |
| `verificacion_simbolica.py` | Comprobación simbólica (SymPy): conservación, jacobiano, valores propios, pico, integral primera, umbral. |
| `resultados_simulacion.csv` | Serie diaria S, I, R de cada escenario (formato largo, 303 filas). |
| `resumen_escenarios.csv` | Tabla comparativa: R₀, pico, día del pico, valores finales, teoría vs. simulación. |
| `sensibilidad_beta_gamma.csv` | Sensibilidad a β y γ (un parámetro a la vez). |
| `barrido_intervencion.csv` | Pico y tamaño del brote en función de u. |
| `verificacion_automatica.csv` | Resultado de las comprobaciones internas. |
| `graficas/` | Todas las figuras en PNG. |

## Instalación

Requiere Python 3 (probado con 3.12; no se verificó en otras versiones).

```bash
python -m venv .venv
# Windows:      .venv\Scripts\activate
# Linux / macOS: source .venv/bin/activate
pip install -r requirements.txt
```

Versiones con las que se probó: numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2, pandas 3.0.6, sympy 1.14.0, streamlit 1.67.0.

## Ejecución

**1. Programa completo con los escenarios del proyecto (A: u=0, B: u=0.5, C: u=0.8)**

```bash
python programa_malware_sir.py
```

Genera los CSV y la carpeta `graficas/` en la carpeta actual e imprime la tabla de resultados y las verificaciones.

**2. Con otros parámetros**

```bash
python programa_malware_sir.py --beta 0.4 --gamma 0.2 --u 0 0.3 0.6 --N 5000 --I0 20 --tmax 150 --salida mis_resultados
python programa_malware_sir.py --help        # todas las opciones
python programa_malware_sir.py --mostrar     # abre ventanas de gráficas
```

Códigos de salida: 0 correcto · 1 alguna verificación falló o error inesperado · 2 entradas inválidas · 3 fallo del integrador · 4 error de archivo.

**3. Interfaz gráfica (para la sustentación)**

```bash
streamlit run interfaz_streamlit.py
```

Permite cambiar N, I(0), β, γ, tiempo de simulación y el porcentaje de intervención u, y muestra
gráficas, tabla, verificaciones y sensibilidad; permite descargar el CSV.

**4. Verificación simbólica**

```bash
python verificacion_simbolica.py
```

## Qué verifica el programa automáticamente

Integración exitosa · conservación S+I+R=N (error < 1e-12·N) · no negatividad · monotonía de S y R · R₀ de cada escenario
(3, 1.5 y 0.6 con los parámetros base) · tasa de crecimiento inicial numérica vs. teórica · coherencia
Reff(0) vs. 1 · pico numérico vs. teórico · condición S = γN/β_ef en el pico · constancia de la integral
primera I + S − S\* ln S · validez de S(∞) (función W de Lambert) · umbrales de intervención ·
coincidencia del CSV releído con los datos en memoria · existencia de los PNG.

## Notas de interpretación

* **Resultados teóricos** (R₀, S\*, pico, S(∞)) vs. **numéricos** (RK45): la tabla `resumen_escenarios.csv`
  los muestra lado a lado. S(t_max) y S(∞) no coinciden cuando el brote no ha terminado dentro del
  horizonte (p. ej. el escenario B a los 100 días).
* En el escenario C (R₀(u) < 1) **no hay pico epidémico**: el máximo de I ocurre en t = 0.
* El CSV usa punto decimal y coma como separador. En Excel con configuración regional en español,
  impórtelo con *Datos → Desde texto/CSV* para evitar errores de separador.
* El modelo supone mezcla homogénea y no incluye topología de red; es una aproximación, no una
  descripción de una red real.
