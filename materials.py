"""
materials.py
============
- Módulo mínimo y robusto para reflectividad.
- Si material == "manual", usa rho_manual.
- Si material está en MATERIALS_NK, calcula ρ con Fresnel (incidencia normal).
"""

import math

# Si más adelante querés usar (n,k) reales a ~311 nm, agregalos aquí:
# Ejemplo:
# MATERIALS_NK = {
#     "Al_311nm": (n_al, k_al),
#     "Ag_311nm": (n_ag, k_ag),
# }
MATERIALS_NK = {}

def rho_from_nk_normal(n: float, k: float) -> float:
    """
    Fresnel (conductor) en incidencia normal:
    R = ((n - 1)^2 + k^2) / ((n + 1)^2 + k^2)
    """
    num = (n - 1.0)**2 + k**2
    den = (n + 1.0)**2 + k**2
    return num / den if den != 0 else 0.0

def get_rho(material_key: str, rho_manual: float) -> float:
    """
    Devuelve la reflectividad efectiva:
    - 'manual' => rho_manual
    - material en MATERIALS_NK => Fresnel a incidencia normal
    - si no encuentra => rho_manual (fallback)
    """
    if material_key == "manual":
        return float(rho_manual)

    if material_key in MATERIALS_NK:
        n, k = MATERIALS_NK[material_key]
        return rho_from_nk_normal(n, k)

    return float(rho_manual)

