# ui.py (Versión 2.1 - PowerNorm pass-through)
# — un solo canvas: fondo + (rayos | curvas)
# - NUEVO: Se pasan todos los parámetros a las funciones de trazado usando **params.

import traceback
import math
from js import document, window, console
from pyscript import when
from pyodide.ffi import create_proxy

from materials import get_rho
from raytracer2d import trace_rays_simple
from quantitative_tracer import trace_for_analysis
from analysis import process_hits
from drawing import draw_main_scene, draw_zoom_rectangle, get_scene_transform
from params import read_common_params, calculate_central_direct_irradiance

# (El diccionario BIND y las funciones _get_el, _get_val se eliminan porque ahora están en params.py y dom_bindings.py)

def _ensure_initial_view(params):
    return {
        "x_min": -params["W_base"]/2 - 0.05, # Pequeño margen
        "x_max":  params["W_base"]/2 + 0.05,
        "z_min":  -0.02,
        "z_max":  params["H"] + 0.02,
    }

def set_status(msg_html: str):
    element = document.getElementById("results-summary-text")
    if element:
        element.innerHTML = msg_html

def set_summary_text(text: str):
    e1 = document.getElementById("irradiance-summary")
    if e1:
        e1.innerText = text
    e2 = document.getElementById("results-summary-text")
    if e2:
        e2.innerHTML = f"<p>{text}</p>"

view = {}
last_rays_cache = []
last_geom_cache = {}
is_dragging = False
drag_start, drag_end = {}, {}
curves_cache = None
MODE_RAYS = "RAYS"
MODE_CURVES = "CURVES"

@when("click", selector="#run_viz")
def run_visualization(event=None):
    try:
        params = read_common_params()
        nrays = int(document.getElementById("nrays_viz").value)
        
        global last_rays_cache, last_geom_cache, view, curves_cache
        
        # --- INICIO: LLAMADA CORREGIDA ---
        last_rays_cache = trace_rays_simple(
            rho=params["rho_val"], nrays=nrays, **params
        )
        # --- FIN: LLAMADA CORREGIDA ---

        view = _ensure_initial_view(params)
        curves_cache = None
        last_geom_cache = {**params, "render_mode": MODE_RAYS}
        set_summary_text("Modo RAYOS. Presioná «Calcular Irradiancia» para ver las curvas.")
        redraw_scene()
    except Exception:
        tb = traceback.format_exc()
        console.error("run_visualization exception:\n", tb)
        set_status(f"<pre style='color:#b00;white-space:pre-wrap'>Error en Visualizar Rayos:\n{tb}</pre>")

@when("click", selector="#run_calc")
def run_calculation(event=None):
    try:
        params = read_common_params()
        nrays_calc = int(document.getElementById("nrays_calc").value)
        
        P_UVB_selected = params["P_UVB"]
        L_tube_selected = params["L_tube"]
        h_tube_current = params["h_tube"]

        sensor_width_m = 0.002
        num_sensors = max(16, int(round(params["W_base"] / sensor_width_m)))
        # --- Normalización por potencia (2D) ---
        phi = float(params["P_UVB"])
        L   = float(params["L_tube"])
        q_W_per_m = phi / max(L, 1e-12)
        N   = int(nrays_calc)
        POWER_PER_RAY_P0 = q_W_per_m / max(N, 1)
        BIN_WIDTH_M = float(params["W_base"]) / float(num_sensors)

        E_central_direct_theoretical_wm2 = calculate_central_direct_irradiance(P_UVB_selected, L_tube_selected, h_tube_current) if params.get("CALIBRATE_TO_E0", True) else 1.0

        hits_for_calibration = trace_for_analysis(
            rho=params["rho_val"], nrays=nrays_calc,
            POWER_PER_RAY_P0=POWER_PER_RAY_P0, BIN_WIDTH_M=BIN_WIDTH_M,
            **params
        )
        
        results_for_calibration = process_hits(hits_for_calibration, params["W_base"], num_sensors)
 
        raw_direct_profile_ua_simulated = results_for_calibration["raw_profiles"]["direct"]
        central_bin_index = num_sensors // 2

        # --- PATCH CALIBRACIÓN ROBUSTA ---
        # Usamos una ventana central de ±2 bins y aplicamos un piso mínimo
        k = 2  # ancho de la ventana (puede ajustar según bins/rayos)
        i0 = max(0, central_bin_index - k)
        i1 = min(len(raw_direct_profile_ua_simulated), central_bin_index + k + 1)
        window_vals = raw_direct_profile_ua_simulated[i0:i1]
        simulated_central_direct_ua = sum(window_vals) / max(1, len(window_vals))
        piso = 1e-3
        if simulated_central_direct_ua < piso:
            print(f"[ADVERTENCIA] El bin central ({central_bin_index}) y su ventana están vacíos o bajos, usando piso de {piso}.")
            simulated_central_direct_ua = piso
        calibration_factor_wm2 = (E_central_direct_theoretical_wm2 / simulated_central_direct_ua) if params.get("CALIBRATE_TO_E0", True) else 1.0
        
        base_hits = trace_for_analysis(
            rho=params["rho_val"], nrays=nrays_calc,
            POWER_PER_RAY_P0=POWER_PER_RAY_P0, BIN_WIDTH_M=BIN_WIDTH_M,
            **params
        )
        
        results = process_hits(base_hits, params["W_base"], num_sensors)
        
        # El resto de la función sigue igual...
        raw_profiles_wm2 = {
            k: [e * calibration_factor_wm2 for e in v] 
            for k, v in results["raw_profiles"].items()
        }

        W_base, H, h_tube, r_tube = params["W_base"], params["H"], params["h_tube"], params["radio_m"]
        n = len(results["profiles"]["total"])
        bin_w = W_base / max(n, 1)
        x_centers = [-W_base/2 + (i + 0.5) * bin_w for i in range(n)]

        Z_bottom = 0.01
        Z_top = min(h_tube - r_tube - 0.01*H, 0.95*H)
        if not (Z_top > Z_bottom + 1e-6): Z_top = max(0.2*H, 0.1)

        y_min_real, y_max_real = 0.0, max(raw_profiles_wm2["total"]) if raw_profiles_wm2["total"] else 1.0
        if (y_max_real - y_min_real) < 1e-6: y_max_real = y_min_real + 1.0
        denom_real = max(1e-6, (y_max_real - y_min_real))
        
        def map_y_to_Z_real(y):
            t = (y - y_min_real) / denom_real
            return Z_bottom + t*(Z_top - Z_bottom)
        
        curves_real = {k: [(x_centers[i], map_y_to_Z_real(v[i])) for i in range(n)] for k, v in raw_profiles_wm2.items()}
        
        global curves_cache, last_geom_cache, view
        curves_cache = curves_real
        if not view: view = _ensure_initial_view(params)
        
        last_geom_cache = {
            **params, "render_mode": MODE_CURVES,
            "y_axis_real": {"Z_bottom": Z_bottom, "Z_top": Z_top, "y_min_real": y_min_real, "y_max_real": y_max_real, "unit_label": "Irradiancia (W/m²)"},
            "calibrated_profiles_wm2": raw_profiles_wm2,
        }

        min_e_real, max_e_real = min(raw_profiles_wm2["total"], default=0.0), max(raw_profiles_wm2["total"], default=0.0)
        summary_text = f"Irradiancia Absoluta — Mín: {min_e_real:.3f} W/m², Máx: {max_e_real:.3f} W/m² (Tubo: {params['tube_label']})"
        
        inflection_x = results.get("inflection_x")
        if inflection_x is not None:
            summary_text += f" | Pto. de inflexión (±x): ±{inflection_x:.4f} m"
            last_geom_cache["inflection_x"] = inflection_x

        set_summary_text(summary_text)
        redraw_scene()
    except Exception:
        tb = traceback.format_exc()
        console.error("run_calculation exception:\n", tb)
        set_status(f"<pre style='color:#b00;white-space:pre-wrap'>Error en Calcular Irradiancia:\n{tb}</pre>")

def redraw_scene():
    if not last_geom_cache: return
    draw_main_scene(last_rays_cache, curves_cache, view, last_geom_cache)
    draw_zoom_rectangle(is_dragging, drag_start, drag_end)

def on_mouse_down(event):
    global is_dragging, drag_start, drag_end
    is_dragging = True
    drag_start = {"x": event.offsetX, "y": event.offsetY}
    drag_end = {"x": event.offsetX, "y": event.offsetY}
    redraw_scene()

def on_mouse_move(event):
    global drag_end
    if is_dragging:
        drag_end = {"x": event.offsetX, "y": event.offsetY}
        redraw_scene()

def on_mouse_up(event):
    global is_dragging, view
    if not is_dragging: return
    is_dragging = False
    if abs(drag_end["x"]-drag_start["x"]) < 10 or abs(drag_end["y"]-drag_start["y"]) < 10:
        redraw_scene()
        return

    canvas = document.getElementById("ray-canvas")
    w, h = int(canvas.clientWidth), int(canvas.clientHeight)
    if w == 0 or h == 0: return

    scale, offset_x, offset_y = get_scene_transform(view, w, h)
    if scale == 0: return

    px_min_x, px_max_x = min(drag_start["x"], drag_end["x"]), max(drag_start["x"], drag_end["x"])
    px_min_y, px_max_y = min(drag_start["y"], drag_end["y"]), max(drag_start["y"], drag_end["y"])

    view = {
        "x_min": view["x_min"] + (px_min_x - offset_x) / scale,
        "x_max": view["x_min"] + (px_max_x - offset_x) / scale,
        "z_max": view["z_max"] - (px_min_y - offset_y) / scale,
        "z_min": view["z_max"] - (px_max_y - offset_y) / scale,
    }
    redraw_scene()

@when("click", selector="#reset_zoom")
def reset_zoom_view(event=None):
    global view
    if last_geom_cache:
        view = _ensure_initial_view(last_geom_cache)
        redraw_scene()

@when('change', selector='#reflector_enabled')
def toggle_reflector_params(event):
    reflector_params_div = document.getElementById('reflector_params')
    reflector_params_div.style.display = 'block' if event.target.value == 'si' else 'none'

def setup_zoom_events():
    canvas = document.getElementById("ray-canvas")
    if canvas:
        canvas.addEventListener("mousedown", create_proxy(on_mouse_down))
        canvas.addEventListener("mousemove", create_proxy(on_mouse_move))
        window.addEventListener("mouseup", create_proxy(on_mouse_up))

setup_zoom_events()
# Inicializa el estado del div de parámetros del reflector al cargar la página
toggle_reflector_params(type('obj', (object,), {'target': document.getElementById('reflector_enabled')})())