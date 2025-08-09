# quantitative_tracer.py (Versión 2.9 – PowerNormRadiance)
# Salida: base_hits[*]["weight"] en W/m² (modelo 2D: por metro de profundidad).
# Cambios clave respecto a versiones previas:
#  - Emisión isotrópica 2π desde el eje del tubo (consistente con fuente lineal 3D).
#  - Directa (b==0): estimador con 1/(4π r²) · cosθ / ΔA  (cuadra con analítico de línea finita).
#  - Rebotes (b>=1): estimador conservador de radiancia por defecto (cosθ / ΔA).
#  - Normalización absoluta por potencia: cada rayo arranca con p0 = (Φ/L)/N [W/m].
#  - Se atenúa potencia por reflectividad en cada rebote: p ← p · ρ.
#  - r² mínimo físico para evitar explosiones numéricas cuando r → 0.

import math

# Utilidades geométricas 2D
def normalize(v):
    m = math.sqrt(v[0]*v[0] + v[1]*v[1])
    return (v[0]/m, v[1]/m) if m > 0 else (0.0, 0.0)

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1]

def reflect(v, n):
    d = dot(v, n)
    return (v[0] - 2.0*d*n[0], v[1] - 2.0*d*n[1])

def trace_for_analysis(
    W_base, W_ceil, H, h_tube, radio_m, rho, nrays, max_bounces,
    reflector_enabled=False, w_reflector_attach=0.0, h_reflector_vertex=0.0,
    **kwargs
):
    """
    Devuelve lista de impactos en base:
      [{"x": x_hit, "bounce": b, "weight": contrib_W_per_m2}, ...]
    """

    # ----- Flags y kernel -----
    PHYSICS_WEIGHTING     = bool(kwargs.get("PHYSICS_WEIGHTING", True))
    KERNEL_MODE           = str(kwargs.get("KERNEL_MODE", "radiance_conserving"))
    BOUNCE_KERNEL_EXPONENT= float(kwargs.get("BOUNCE_KERNEL_EXPONENT", 1.0))

    # ----- Interacción del tubo (opcional) -----
    TUBE_INTERACTS = bool(kwargs.get("TUBE_INTERACTS", True))
    TUBE_OCCLUDES  = bool(kwargs.get("TUBE_OCCLUDES",  True))
    TUBE_FRESNEL_R = float(kwargs.get("TUBE_FRESNEL_R", 0.04))

    # ----- Normalización por potencia (OBLIGATORIA) -----
    # POWER_PER_RAY_P0: potencia lineal por rayo [W/m] = (Φ/L)/N
    # BIN_WIDTH_M: ancho del bin en base [m]  → ΔA = BIN_WIDTH_M · 1 m (2D)
    if "POWER_PER_RAY_P0" not in kwargs or "BIN_WIDTH_M" not in kwargs:
        raise ValueError("Faltan kwargs obligatorios: POWER_PER_RAY_P0 y/o BIN_WIDTH_M")
    POWER_PER_RAY_P0 = float(kwargs["POWER_PER_RAY_P0"])
    BIN_WIDTH_M      = float(kwargs["BIN_WIDTH_M"])
    DELTA_A          = BIN_WIDTH_M * 1.0  # [m² por metro de profundidad]

    # ----- Geometría -----
    base_n = (0.0, 1.0)
    half_W_base, half_W_ceil = 0.5*W_base, 0.5*W_ceil

    # Pared derecha: segmento (p1_r → p2_r) y normal saliente
    p1_r, p2_r = ( half_W_base, 0.0), ( half_W_ceil, H)
    v_r = (p2_r[0]-p1_r[0], p2_r[1]-p1_r[1])
    n_r = normalize(( v_r[1], -v_r[0]))

    # Pared izquierda
    p1_l, p2_l = (-half_W_base, 0.0), (-half_W_ceil, H)
    v_l = (p2_l[0]-p1_l[0], p2_l[1]-p1_l[1])
    n_l = normalize((-v_l[1],  v_l[0]))

    # Reflector en V (opcional)
    if reflector_enabled:
        pv = (0.0, h_reflector_vertex)
        pl = (-0.5*w_reflector_attach, H)
        pr = ( 0.5*w_reflector_attach, H)
        v_rl = (pv[0]-pl[0], pv[1]-pl[1]); n_rl = normalize((-v_rl[1],  v_rl[0]))
        v_rr = (pr[0]-pv[0], pr[1]-pv[1]); n_rr = normalize(( v_rr[1], -v_rr[0]))

    # ----- Clamp físico de distancia -----
    # r²_floor ~ 10% de (h_tube - radio_m)² para evitar r → 0
    r2_floor = max(1e-6, 0.1 * (max(h_tube - radio_m, 0.0)**2))

    base_hits = []
    x0, z0 = 0.0, h_tube  # eje del tubo

    # ===== Emisión isotrópica 2π desde el eje (modelo de línea 3D) =====
    for i in range(int(nrays)):
        ang = 2.0 * math.pi * (i + 0.5) / float(nrays)  # estratificación simple
        dx, dz = math.cos(ang), math.sin(ang)
        x1, z1 = x0, z0

        weight = 1.0                # peso relativo (compatibilidad)
        power  = POWER_PER_RAY_P0   # [W/m] potencia por rayo en 2D

        # Secuencia de rebotes hasta impactar base o terminar
        for b in range(int(max_bounces) + 1):
            candidates = []

            # Base (z=0)
            if dz < -1e-6:
                t = -z1 / dz
                if t > 1e-6:
                    xb = x1 + dx*t
                    if abs(xb) <= half_W_base:
                        candidates.append(("base", t, xb, 0.0))

            # Techo (z=H)
            if dz > 1e-6:
                t = (H - z1) / dz
                if t > 1e-6:
                    xc = x1 + dx*t
                    if abs(xc) <= half_W_ceil:
                        candidates.append(("ceiling", t, xc, H))

            # Pared derecha
            det_r = dx*v_r[1] - dz*v_r[0]
            if abs(det_r) > 1e-6:
                t_r = ((p1_r[0]-x1)*v_r[1] - (p1_r[1]-z1)*v_r[0]) / det_r
                s_r = ((p1_r[0]-x1)*dz     - (p1_r[1]-z1)*dx    ) / det_r
                if t_r > 1e-6 and 0.0 <= s_r <= 1.0:
                    candidates.append(("wall_r", t_r, None, None))

            # Pared izquierda
            det_l = dx*v_l[1] - dz*v_l[0]
            if abs(det_l) > 1e-6:
                t_l = ((p1_l[0]-x1)*v_l[1] - (p1_l[1]-z1)*v_l[0]) / det_l
                s_l = ((p1_l[0]-x1)*dz     - (p1_l[1]-z1)*dx    ) / det_l
                if t_l > 1e-6 and 0.0 <= s_l <= 1.0:
                    candidates.append(("wall_l", t_l, None, None))

            # Reflector en V (si habilitado)
            if reflector_enabled:
                det = dx * v_rl[1] - dz * v_rl[0]
                if abs(det) > 1e-6:
                    t = ((pl[0]-x1)*v_rl[1] - (pl[1]-z1)*v_rl[0]) / det
                    s = ((pl[0]-x1)*dz      - (pl[1]-z1)*dx     ) / det
                    if t > 1e-6 and 0.0 <= s <= 1.0:
                        candidates.append(("refl_l", t, None, None))
                det = dx * v_rr[1] - dz * v_rr[0]
                if abs(det) > 1e-6:
                    t = ((pv[0]-x1)*v_rr[1] - (pv[1]-z1)*v_rr[0]) / det
                    s = ((pv[0]-x1)*dz      - (pv[1]-z1)*dx     ) / det
                    if t > 1e-6 and 0.0 <= s <= 1.0:
                        candidates.append(("refl_r", t, None, None))

            # --- Tubo (círculo) como obstáculo con Fresnel constante ---
            if TUBE_INTERACTS and TUBE_OCCLUDES:
                # Sólo si el punto actual está fuera del tubo
                if (x1 - x0)*(x1 - x0) + (z1 - z0)*(z1 - z0) > (radio_m + 1e-6)**2:
                    ox, oz = (x1 - x0), (z1 - z0)
                    a = dx*dx + dz*dz
                    bq = 2.0*(ox*dx + oz*dz)
                    c = ox*ox + oz*oz - radio_m*radio_m
                    disc = bq*bq - 4.0*a*c
                    if disc > 0.0:
                        sdisc = math.sqrt(disc)
                        t1 = (-bq - sdisc) / (2.0*a)
                        t2 = (-bq + sdisc) / (2.0*a)
                        t_hit = min(t for t in (t1, t2) if t > 1e-6) if (t1 > 1e-6 or t2 > 1e-6) else None
                        if t_hit is not None:
                            xt, zt = (x1 + dx*t_hit, z1 + dz*t_hit)
                            candidates.append(("tube", t_hit, xt, zt))

            if not candidates:
                break  # no intersecta nada más

            # Elegir el primer impacto
            surf, tmin, x_hit, z_hit = min(candidates, key=lambda c: c[1])
            x2, z2 = (x1 + dx*tmin, z1 + dz*tmin)

            # ----- Contribución en base -----
            if surf == "base":
                cos_inc = max(0.0, -dot((dx, dz), base_n))  # (= -dz si dz<0)

                if not PHYSICS_WEIGHTING:
                    contrib = weight  # unidades arbitrarias (modo legado)
                else:
                    if b == 0:
                        # DIRECTA: 1/(4π r²) · cosθ / ΔA
                        r2 = max(tmin*tmin, r2_floor)
                        contrib = (power * cos_inc / (4.0 * math.pi * r2)) / DELTA_A
                    else:
                        # REBOTES: según KERNEL_MODE
                        if KERNEL_MODE == "inverse_square_all":
                            r2 = max(tmin*tmin, r2_floor)
                            contrib = (power * cos_inc / r2) / DELTA_A
                        elif KERNEL_MODE == "direct_inv_square_bounce_rp":
                            # Atenuación geométrica parametrizada (p≈1…2)
                            r_p = max(tmin, math.sqrt(r2_floor)) ** max(1.0, BOUNCE_KERNEL_EXPONENT)
                            contrib = (power * cos_inc / r_p) / DELTA_A
                        elif KERNEL_MODE == "radiance_conserving":
                            contrib = (power * cos_inc) / DELTA_A
                        else:
                            # Fallback seguro: 1/r²
                            r2 = max(tmin*tmin, r2_floor)
                            contrib = (power * cos_inc / r2) / DELTA_A

                base_hits.append({"x": x2, "bounce": b, "weight": contrib})
                break  # termina en la base

            # ----- Impacto en tubo: refleja R y absorbe 1-R; sin transmisión -----
            if surf == "tube" and TUBE_INTERACTS and TUBE_OCCLUDES:
                # Normal del círculo
                nx, nz = (x2 - x0), (z2 - z0)
                nm = math.hypot(nx, nz)
                if nm <= 0.0:
                    break
                nx, nz = nx/nm, nz/nm

                # Multiplicar energía por R constante; absorber 1-R
                weight *= TUBE_FRESNEL_R
                power  *= TUBE_FRESNEL_R
                if weight < 0.01:
                    break

                # Reflejar direccionalmente y avanzar un epsilon
                dx, dz = reflect((dx, dz), (nx, nz))
                eps = 1e-5
                x1, z1 = (x2 + eps*dx, z2 + eps*dz)
                continue  # no aplicar rho del recinto

            # ----- Atenuación por reflectividad y corte -----
            weight *= rho
            power  *= rho
            if weight < 0.01:
                break

            # ----- Reflexión especular -----
            if   surf == "ceiling": n = (0.0, -1.0)
            elif surf == "wall_r":  n = n_r
            elif surf == "wall_l":  n = n_l
            elif surf == "refl_l":  n = n_rl
            elif surf == "refl_r":  n = n_rr
            else: break

            dx, dz = reflect((dx, dz), n)
            x1, z1 = x2, z2

    return base_hits
