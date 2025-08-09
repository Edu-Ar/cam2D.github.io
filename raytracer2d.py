# raytracer2d.py (Versión 1.2 - Tubo opaco + Fresnel constante)
# ======================================
# Trazador de rayos 2D para visualización.
# - Paredes y techo reflectores, piso absorbente.
# - Reflector interno en V opcional.
# - NUEVO: Oclusión y reflexión especular (Fresnel constante) en el tubo emisor (visual).
# ======================================

import math

# --- Funciones de ayuda para vectores ---
def normalize(v):
    mag = math.hypot(v[0], v[1])
    return (v[0]/mag, v[1]/mag) if mag > 0 else v

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1]

def reflect(dir_vec, normal):
    # dir' = dir - 2 (dir·n) n
    d_dot_n = dot(dir_vec, normal)
    return (dir_vec[0] - 2.0*d_dot_n*normal[0],
            dir_vec[1] - 2.0*d_dot_n*normal[1])

# ------------------------------------------------------------
# Función principal (visualización de rayos)
# ------------------------------------------------------------
def trace_rays_simple(W_base, W_ceil, H, h_tube, radio_m, rho, nrays, max_bounces,
                      reflector_enabled=False, w_reflector_attach=0, h_reflector_vertex=0, **kwargs):
    """
    Devuelve una lista de 'rays', donde cada rayo es una lista de segmentos (x1, z1, x2, z2, bounce_index).
    Esta función es para visualización (no pondera energía).
    """

    rays = []

    # --- Switches del tubo (opcionales; por defecto interactúa con R constante) ---
    TUBE_INTERACTS = bool(kwargs.get("TUBE_INTERACTS", True))
    TUBE_OCCLUDES  = bool(kwargs.get("TUBE_OCCLUDES",  True))
    TUBE_FRESNEL_R = float(kwargs.get("TUBE_FRESNEL_R", 0.04))

    # Centro del tubo (en 2D)
    x_center, z_center = 0.0, h_tube

    # Geometría de paredes (trapecio)
    half_W_base = 0.5 * W_base
    half_W_ceil = 0.5 * W_ceil
    p_base_r = ( half_W_base, 0.0)
    p_base_l = (-half_W_base, 0.0)
    p_ceil_r = ( half_W_ceil, H)
    p_ceil_l = (-half_W_ceil, H)

    # Vectores y normales de paredes
    v_wall_r = (p_ceil_r[0]-p_base_r[0], p_ceil_r[1]-p_base_r[1])
    v_wall_l = (p_ceil_l[0]-p_base_l[0], p_ceil_l[1]-p_base_l[1])
    normal_r = normalize(( v_wall_r[1], -v_wall_r[0]))  # normal saliente
    normal_l = normalize((-v_wall_l[1],  v_wall_l[0]))  # normal saliente (simétrica)

    # Reflector en V (opcional)
    if reflector_enabled:
        x_attach = 0.5 * w_reflector_attach
        p_refl_l = (-x_attach, H)
        p_refl_r = ( x_attach, H)
        p_refl_v = (0.0, h_reflector_vertex)
        v_refl_l = (p_refl_v[0]-p_refl_l[0], p_refl_v[1]-p_refl_l[1])
        v_refl_r = (p_refl_v[0]-p_refl_r[0], p_refl_v[1]-p_refl_r[1])
        normal_refl_l = normalize(( v_refl_l[1], -v_refl_l[0]))
        normal_refl_r = normalize((-v_refl_r[1],  v_refl_r[0]))

    # Emisión desde el contorno del tubo (ángulo uniforme en 2D)
    for i in range(int(nrays)):
        angle = 2.0 * math.pi * i / max(1, int(nrays))
        dx, dz = math.cos(angle), math.sin(angle)
        x1, z1 = x_center + radio_m * dx, z_center + radio_m * dz

        segments = []

        for b in range(int(max_bounces) + 1):
            candidates = []

            # --- Base (z=0) absorbente (visual) ---
            if dz < -1e-6:
                t = -z1 / dz
                if t > 1e-6:
                    xb = x1 + dx * t
                    if abs(xb) <= half_W_base:
                        candidates.append({"surf": "base", "t": t})

            # --- Techo (z=H) reflectante ---
            if dz > 1e-6:
                t = (H - z1) / dz
                if t > 1e-6:
                    xc = x1 + dx * t
                    if abs(xc) <= half_W_ceil:
                        candidates.append({"surf": "ceiling", "t": t})

            # --- Pared derecha (segmento) ---
            det_r = dx * v_wall_r[1] - dz * v_wall_r[0]
            if abs(det_r) > 1e-6:
                t_r = ((p_base_r[0]-x1) * v_wall_r[1] - (p_base_r[1]-z1) * v_wall_r[0]) / det_r
                s_r = ((p_base_r[0]-x1) * dz          - (p_base_r[1]-z1) * dx         ) / det_r
                if t_r > 1e-6 and 0.0 <= s_r <= 1.0:
                    candidates.append({"surf": "wall_r", "t": t_r})

            # --- Pared izquierda (segmento) ---
            det_l = dx * v_wall_l[1] - dz * v_wall_l[0]
            if abs(det_l) > 1e-6:
                t_l = ((p_base_l[0]-x1) * v_wall_l[1] - (p_base_l[1]-z1) * v_wall_l[0]) / det_l
                s_l = ((p_base_l[0]-x1) * dz          - (p_base_l[1]-z1) * dx         ) / det_l
                if t_l > 1e-6 and 0.0 <= s_l <= 1.0:
                    candidates.append({"surf": "wall_l", "t": t_l})

            # --- Reflector en V (2 segmentos) ---
            if reflector_enabled:
                det_rl = dx * v_refl_l[1] - dz * v_refl_l[0]
                if abs(det_rl) > 1e-6:
                    t_rl = ((p_refl_l[0]-x1) * v_refl_l[1] - (p_refl_l[1]-z1) * v_refl_l[0]) / det_rl
                    s_rl = ((p_refl_l[0]-x1) * dz          - (p_refl_l[1]-z1) * dx         ) / det_rl
                    if t_rl > 1e-6 and 0.0 <= s_rl <= 1.0:
                        candidates.append({"surf": "refl_l", "t": t_rl})
                det_rr = dx * v_refl_r[1] - dz * v_refl_r[0]
                if abs(det_rr) > 1e-6:
                    t_rr = ((p_refl_r[0]-x1) * v_refl_r[1] - (p_refl_r[1]-z1) * v_refl_r[0]) / det_rr
                    s_rr = ((p_refl_r[0]-x1) * dz          - (p_refl_r[1]-z1) * dx         ) / det_rr
                    if t_rr > 1e-6 and 0.0 <= s_rr <= 1.0:
                        candidates.append({"surf": "refl_r", "t": t_rr})

            # --- Tubo (círculo) como obstáculo con Fresnel constante (visual) ---
            if TUBE_INTERACTS and TUBE_OCCLUDES:
                dx0, dz0 = (x1 - x_center), (z1 - z_center)
                if dx0*dx0 + dz0*dz0 > (radio_m + 1e-6)**2:
                    ox, oz = (x1 - x_center), (z1 - z_center)
                    a = dx*dx + dz*dz
                    bq = 2.0 * (ox*dx + oz*dz)
                    c = ox*ox + oz*oz - radio_m*radio_m
                    disc = bq*bq - 4.0*a*c
                    if disc > 0.0:
                        sqrt_disc = math.sqrt(disc)
                        t1 = (-bq - sqrt_disc) / (2.0*a)
                        t2 = (-bq + sqrt_disc) / (2.0*a)
                        t_hit = min(t for t in (t1, t2) if t > 1e-6) if (t1 > 1e-6 or t2 > 1e-6) else None
                        if t_hit is not None:
                            candidates.append({"surf": "tube", "t": t_hit})

            if not candidates:
                break  # no intersecta nada más

            # Elegir el primer impacto
            best_hit = min(candidates, key=lambda v: v["t"])
            t_min, surf = best_hit["t"], best_hit["surf"]
            x_hit, z_hit = (x1 + dx * t_min, z1 + dz * t_min)

            # Registrar el segmento hasta el punto de impacto
            segments.append((x1, z1, x_hit, z_hit, b))

            # --- Superficies ---
            if surf == "base":
                break  # absorbente (visual)

            if surf == "tube" and TUBE_INTERACTS and TUBE_OCCLUDES:
                nx, nz = (x_hit - x_center), (z_hit - z_center)
                nm = math.hypot(nx, nz)
                if nm > 0.0:
                    nx, nz = nx/nm, nz/nm
                    if TUBE_FRESNEL_R <= 0.0:
                        break  # opaco puro
                    dx, dz = reflect((dx, dz), (nx, nz))
                    eps = 1e-5
                    x1, z1 = (x_hit + eps*dx, z_hit + eps*dz)
                    continue  # no aplica rho del recinto aquí

            if   surf == "ceiling":
                n = (0.0, -1.0)
            elif surf == "wall_r":
                n = normal_r
            elif surf == "wall_l":
                n = normal_l
            elif surf == "refl_l":
                n = normal_refl_l
            elif surf == "refl_r":
                n = normal_refl_r
            else:
                break

            dx, dz = reflect((dx, dz), n)
            x1, z1 = x_hit, z_hit

        rays.append(segments)

    return rays
