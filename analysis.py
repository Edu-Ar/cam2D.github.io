# analysis.py (Versión 1.2 - Devuelve perfiles sin normalizar)
import math

def find_inflection_point(profile, W_base):
    """Devuelve |x| del primer cruce de curvatura (d2) cercano al centro.
    Se suaviza con una media móvil corta para robustez.
    """
    N = len(profile)
    if N < 5:
        return None
    # Suavizado (ventana impar)
    k = 5
    half = k // 2
    smoothed = []
    for i in range(N):
        a = max(0, i - half)
        b = min(N, i + half + 1)
        smoothed.append(sum(profile[a:b]) / (b - a))
    # Segunda derivada discreta
    d2 = [smoothed[i+1] - 2*smoothed[i] + smoothed[i-1] for i in range(1, N-1)]
    # Buscar cruce de signo alrededor del centro
    center = (N - 1) // 2  # índice aproximado en d2
    def sign(x): return 1 if x > 0 else (-1 if x < 0 else 0)
    best_idx = None
    for delta in range(0, len(d2)//2):
        for idx in [center + delta, center - delta]:
            if 1 <= idx < len(d2):
                s1 = sign(d2[idx-1]); s2 = sign(d2[idx])
                if s1 != 0 and s2 != 0 and s1 != s2:
                    best_idx = idx
                    break
        if best_idx is not None:
            break
    if best_idx is None:
        return None
    sensor_width = W_base / N
    # Coordenada x (aprox. entre idx-1 y idx) respecto al centro
    x = (best_idx - (N/2 - 1)) * sensor_width
    return abs(x)

def process_hits(base_hits, W_base, num_sensors):
    if not base_hits or num_sensors == 0:
        zeros = [0] * (num_sensors or 1)
        return {
            "profiles": {"total": zeros, "direct": zeros, "first": zeros, "multi": zeros},
            "raw_profiles": {"total": zeros, "direct": zeros, "first": zeros, "multi": zeros}, # Añadir raw_profiles
            "energy_range": (0, 0), "inflection_x": None
        }

    sensors = { "total": [0]*num_sensors, "direct": [0]*num_sensors, "first": [0]*num_sensors, "multi": [0]*num_sensors }
    bin_width = W_base / num_sensors
    half_W_base = W_base / 2.0

    for hit in base_hits:
        bin_index = int((hit["x"] + half_W_base) / bin_width)
        if 0 <= bin_index < num_sensors:
            energy = hit["weight"]
            sensors["total"][bin_index] += energy
            if hit["bounce"] == 0: sensors["direct"][bin_index] += energy
            elif hit["bounce"] == 1: sensors["first"][bin_index] += energy
            else: sensors["multi"][bin_index] += energy

    # Almacenar los perfiles brutos (sin normalizar)
    raw_sensors = {
        "total":  list(sensors["total"]),
        "direct": list(sensors["direct"]),
        "first":  list(sensors["first"]),
        "multi":  list(sensors["multi"]),
    }

    # Normalización UNIFICADA para los perfiles que se graficarán en el eje derecho (0-1)
    max_total_energy = max(sensors["total"]) if sensors["total"] else 1.0
    scale_factor = max_total_energy if max_total_energy > 0 else 1.0

    profiles_normalized = {
        "total":  [e/scale_factor for e in sensors["total"]],
        "direct": [e/scale_factor for e in sensors["direct"]],
        "first":  [e/scale_factor for e in sensors["first"]],
        "multi":  [e/scale_factor for e in sensors["multi"]],
    }
    
    min_e_scaled = min(profiles_normalized["total"]) if profiles_normalized["total"] else 0.0
    max_e_scaled = max(profiles_normalized["total"]) if profiles_normalized["total"] else 0.0
    
    inflection_x = find_inflection_point(profiles_normalized["total"], W_base)

    return {
        "profiles": profiles_normalized,
        "raw_profiles": raw_sensors,
        "energy_range": (min_e_scaled, max_e_scaled),
        "inflection_x": inflection_x
    }
    