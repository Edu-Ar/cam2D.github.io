# params.py (Versión 1.5 - PowerNorm stable)
# Lee parámetros desde el DOM y arma el diccionario 'params' que usan la UI y el trazador.
# Salida coherente con el modelo 2D por metro de profundidad.

from dom_bindings import get_val_str, get_val_float, get_val_int
from materials import get_rho
import math

# Datos de tubos Philips UV-B (potencia UVB efectiva Φ a 5h IEC y longitud L)
PHILIPS_TUBE_DATA = {
    "TL20W": {
        "label": "Philips TL 20W/01 RS",
        "P_UVB": 3.23,     # [W]
        "L": 0.5898,       # [m]
    },
    "TL40W": {
        "label": "Philips TL 40W/01 RS",
        "P_UVB": 7.70,     # [W]
        "L": 1.1994,       # [m]
    },
}

def read_common_params():
    # --- Entradas principales desde la UI ---
    material     = get_val_str(  "material",   default="manual", required=False)
    rho_manual   = get_val_float("rho",        default=0.8,      required=False)

    W_base       = get_val_float("W_base",     required=True)
    W_ceil       = get_val_float("W_ceil",     required=True)
    H            = get_val_float("H",          required=True)
    h_tube       = get_val_float("h",          required=True)
    max_bounces  = get_val_int(  "max_bounces",required=True)

    # Selección de tubo (Φ, L)
    tube_key = get_val_str("tube_type", default="TL20W", required=True)
    tube = PHILIPS_TUBE_DATA.get(tube_key)
    if not tube:
        raise ValueError(f"Tipo de tubo '{tube_key}' no encontrado en PHILIPS_TUBE_DATA.")

    P_UVB  = float(tube["P_UVB"])   # [W]
    L_tube = float(tube["L"])       # [m]

    # Propiedades ópticas
    rho_val = float(get_rho(material, rho_manual))

    # Geometría del tubo (radio)
    radio_m = 0.02025  # Ø 40.5 mm ⇒ R = 0.02025 m

    # --- Armar diccionario de parámetros ---
    params = {
        "material":     material,
        "rho_val":      rho_val,
        "W_base":       float(W_base),
        "W_ceil":       float(W_ceil),
        "H":            float(H),
        "h_tube":       float(h_tube),
        "radio_m":      float(radio_m),
        "max_bounces":  int(max_bounces),
        "P_UVB":        P_UVB,
        "L_tube":       L_tube,
        "tube_label":   tube["label"],
    }

    # Reflector en V (opcional)
    reflector_enabled_str = get_val_str("reflector_enabled", default="no", required=False)
    params["reflector_enabled"] = (reflector_enabled_str == "si")
    if params["reflector_enabled"]:
        params["w_reflector_attach"] = get_val_float("w_reflector_attach", required=True)
        params["h_reflector_vertex"] = get_val_float("h_reflector_vertex", required=True)

    # --- Flags físicos ---
    params["PHYSICS_WEIGHTING"] = True
    params["CALIBRATE_TO_E0"]   = False

    # Kernel recomendado (rebotes radiance-conserving; la directa usa 1/(4π r^2) en el trazador)
    params["KERNEL_MODE"] = "radiance_conserving"
    # (Opcional) Exponente si usás el modo paramétrico de rebotes:
    # params["BOUNCE_KERNEL_EXPONENT"] = 1.2

    return params

# Analítico de línea finita (2D por metro): E(0) = q / [4 π h sqrt(h^2 + (L/2)^2)], con q = Φ/L
def E0_line_source(phi_e_W, L_m, h_m):
    if h_m <= 0:
        return float("inf")
    q = phi_e_W / L_m  # [W/m]
    return q / (4.0 * math.pi * h_m * math.sqrt(h_m*h_m + (0.5*L_m)**2))

def calculate_central_direct_irradiance(phi_e_W, L_m, h_m):
    # Alias usado por la UI
    return E0_line_source(phi_e_W, L_m, h_m)
