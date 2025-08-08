# dom_bindings.py (Versión 1.1)
# Centraliza los IDs del HTML y helpers mínimos para leer valores del DOM.

from js import document, console

# === EDITABLE: poné aquí los IDs reales de tu index.html ===
BIND = {
    "material":   "material",
    "rho":        "rho_manual",
    "W_base":     "ancho_base",
    "W_ceil":     "ancho_techo",
    "H":          "H",
    "h":          "h",
    "max_bounces":"max_bounces",
    "nrays_viz":  "nrays_viz",
    "nrays_calc": "nrays_calc",
    "canvas":     "ray-canvas",
    "tube_type":  "tube_type", # <--- AÑADIDO: ID para el selector de tipo de tubo
    "reflector_enabled": "reflector_enabled",
    "w_reflector_attach": "w_reflector_attach",
    "h_reflector_vertex": "h_reflector_vertex",    
}

def get_el(key: str):
    """Devuelve el elemento DOM para la clave dada según BIND, o None si no existe."""
    cid = BIND.get(key)
    el = document.getElementById(cid) if cid else None
    return el

def get_val_str(key: str, default=None, required=True) -> str:
    el = get_el(key)
    if el is None:
        if required:
            raise RuntimeError(f"Falta el elemento requerido: #{BIND.get(key, key)} (clave '{key}').")
        console.warn(f"[dom_bindings] No hallé #{BIND.get(key,key)}. Uso default={default!r}")
        return default
    v = getattr(el, "value", None)
    if (v is None or v == "") and required and default is None:
        raise RuntimeError(f"#{BIND[key]} está vacío y no hay default (clave '{key}').")
    return v if v not in (None, "") else default

def get_val_float(key: str, default=None, required=True) -> float:
    v = get_val_str(key, default=None, required=required)
    if v is None:
        return default
    try:
        return float(v)
    except Exception:
        raise RuntimeError(f"Valor inválido en #{BIND.get(key,key)}: {v!r} (esperaba float).")

def get_val_int(key: str, default=None, required=True) -> int:
    v = get_val_str(key, default=None, required=required)
    if v is None:
        return default
    try:
        return int(v)
    except Exception:
        raise RuntimeError(f"Valor inválido en #{BIND.get(key,key)}: {v!r} (esperaba int).")

def diagnose_required(keys=("material","rho","W_base","W_ceil","H","h","max_bounces","nrays_viz","nrays_calc","canvas", "tube_type")): # <--- Añadir aquí también para diagnóstico
    """Devuelve un dict {clave: True/False} indicando si existe cada id en el DOM."""
    status = {}
    for k in keys:
        status[k] = get_el(k) is not None
    return status