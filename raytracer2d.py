# raytracer2d.py (Versión 1.1 - Reflector en V)
# ======================================
# Trazador de rayos 2D para visualización.
# - Paredes y techo reflectores, piso absorbente.
# - NUEVO: Añadida la lógica para un reflector interno en V opcional.
# ======================================
import math

# --- Funciones de ayuda para vectores ---
def normalize(v):
    mag = math.sqrt(v[0]**2 + v[1]**2)
    return (v[0]/mag, v[1]/mag) if mag > 0 else v

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1]

def reflect(v_in, normal):
    d = dot(v_in, normal)
    v_out = (v_in[0] - 2 * d * normal[0], v_in[1] - 2 * d * normal[1])
    return v_out

def trace_rays_simple(W_base, W_ceil, H, h_tube, radio_m, rho, nrays, max_bounces,
                      reflector_enabled=False, w_reflector_attach=0, h_reflector_vertex=0, **kwargs):
    rays = []
    x_center, z_center = 0.0, h_tube

    half_W_base = W_base / 2.0
    half_W_ceil = W_ceil / 2.0
    
    # Geometría de la caja
    p1_r, p2_r = (half_W_base, 0), (half_W_ceil, H)
    v_wall_r = (p2_r[0] - p1_r[0], p2_r[1] - p1_r[1])
    normal_r = normalize((v_wall_r[1], -v_wall_r[0]))

    p1_l, p2_l = (-half_W_base, 0), (-half_W_ceil, H)
    v_wall_l = (p2_l[0] - p1_l[0], p2_l[1] - p1_l[1])
    normal_l = normalize((-v_wall_l[1], v_wall_l[0]))

    # --- INICIO: Geometría del Reflector ---
    if reflector_enabled:
        p_vertex_refl = (0.0, h_reflector_vertex)
        p_left_refl = (-w_reflector_attach / 2.0, H)
        p_right_refl = (w_reflector_attach / 2.0, H)

        v_refl_l = (p_vertex_refl[0] - p_left_refl[0], p_vertex_refl[1] - p_left_refl[1])
        normal_refl_l = normalize((-v_refl_l[1], v_refl_l[0]))

        v_refl_r = (p_right_refl[0] - p_vertex_refl[0], p_right_refl[1] - p_vertex_refl[1])
        normal_refl_r = normalize((v_refl_r[1], -v_refl_r[0]))
    # --- FIN: Geometría del Reflector ---

    for i in range(nrays):
        angle = 2 * math.pi * i / nrays
        dx, dz = math.cos(angle), math.sin(angle)
        
        x1 = x_center + radio_m * dx
        z1 = z_center + radio_m * dz
        
        segments = []
        weight = 1.0

        for b in range(max_bounces + 1):
            candidates = []
            
            # Intersecciones con la caja
            if dz < -1e-6:
                t = -z1 / dz
                if t > 1e-6: candidates.append({"surf": "base", "t": t})
            if dz > 1e-6:
                t = (H - z1) / dz
                if t > 1e-6: candidates.append({"surf": "ceiling", "t": t})

            matrix_det_r = dx * v_wall_r[1] - dz * v_wall_r[0]
            if abs(matrix_det_r) > 1e-6:
                t_r = ((p1_r[0] - x1) * v_wall_r[1] - (p1_r[1] - z1) * v_wall_r[0]) / matrix_det_r
                s_r = ((p1_r[0] - x1) * dz - (p1_r[1] - z1) * dx) / matrix_det_r
                if t_r > 1e-6 and 0 <= s_r <= 1:
                    candidates.append({"surf": "wall_r", "t": t_r})

            matrix_det_l = dx * v_wall_l[1] - dz * v_wall_l[0]
            if abs(matrix_det_l) > 1e-6:
                t_l = ((p1_l[0] - x1) * v_wall_l[1] - (p1_l[1] - z1) * v_wall_l[0]) / matrix_det_l
                s_l = ((p1_l[0] - x1) * dz - (p1_l[1] - z1) * dx) / matrix_det_l
                if t_l > 1e-6 and 0 <= s_l <= 1:
                    candidates.append({"surf": "wall_l", "t": t_l})

            # --- INICIO: Intersecciones con el Reflector ---
            if reflector_enabled:
                matrix_det_refl_l = dx * v_refl_l[1] - dz * v_refl_l[0]
                if abs(matrix_det_refl_l) > 1e-6:
                    t_refl_l = ((p_left_refl[0] - x1) * v_refl_l[1] - (p_left_refl[1] - z1) * v_refl_l[0]) / matrix_det_refl_l
                    s_refl_l = ((p_left_refl[0] - x1) * dz - (p_left_refl[1] - z1) * dx) / matrix_det_refl_l
                    if t_refl_l > 1e-6 and 0 <= s_refl_l <= 1:
                        candidates.append({"surf": "refl_l", "t": t_refl_l})

                matrix_det_refl_r = dx * v_refl_r[1] - dz * v_refl_r[0]
                if abs(matrix_det_refl_r) > 1e-6:
                    t_refl_r = ((p_vertex_refl[0] - x1) * v_refl_r[1] - (p_vertex_refl[1] - z1) * v_refl_r[0]) / matrix_det_refl_r
                    s_refl_r = ((p_vertex_refl[0] - x1) * dz - (p_vertex_refl[1] - z1) * dx) / matrix_det_refl_r
                    if t_refl_r > 1e-6 and 0 <= s_refl_r <= 1:
                        candidates.append({"surf": "refl_r", "t": t_refl_r})
            # --- FIN: Intersecciones con el Reflector ---

            if not candidates: break
            
            # Validación de candidatos (existente)
            valid_candidates = []
            for hit in candidates:
                x_hit_check = x1 + hit["t"] * dx
                z_hit_check = z1 + hit["t"] * dz
                if hit["surf"] == "base" and abs(x_hit_check) <= half_W_base:
                    valid_candidates.append(hit)
                elif hit["surf"] == "ceiling" and abs(x_hit_check) <= half_W_ceil:
                    valid_candidates.append(hit)
                elif "wall" in hit["surf"] and 0 <= z_hit_check <= H:
                    valid_candidates.append(hit)
                elif "refl" in hit["surf"]: # El reflector es siempre válido si se intersecta
                    valid_candidates.append(hit)

            if not valid_candidates: break
            
            best_hit = min(valid_candidates, key=lambda v: v["t"])
            t_min, surf = best_hit["t"], best_hit["surf"]

            x_hit, z_hit = (x1 + dx * t_min, z1 + dz * t_min)
            segments.append((x1, z1, x_hit, z_hit, b))
            
            weight *= rho

            if surf == "base" or weight < 0.01:
                break
            
            # Lógica de reflexión
            if surf == "ceiling": dx, dz = reflect((dx, dz), (0, -1))
            elif surf == "wall_r": dx, dz = reflect((dx, dz), normal_r)
            elif surf == "wall_l": dx, dz = reflect((dx, dz), normal_l)
            # --- INICIO: Lógica de Reflexión del Reflector ---
            elif surf == "refl_l": dx, dz = reflect((dx, dz), normal_refl_l)
            elif surf == "refl_r": dx, dz = reflect((dx, dz), normal_refl_r)
            # --- FIN: Lógica de Reflexión del Reflector ---
            
            x1, z1 = x_hit, z_hit
        
        rays.append(segments)
    return rays
    