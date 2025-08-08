# drawing.py (Versión 1.17)
# - Centra la escena en el canvas y sincroniza tamaño CSS ↔ bitmap.
# - _draw_rays acepta segmentos, pares y polilíneas; convierte a float y valida finitud.
# - Eje Y  real a la izquierda.

import math
from js import document, window
try:
    from js import console
except Exception:
    console = None

# --------------------------
# Transformación de escena
# --------------------------
def get_scene_transform(view, w, h):
    xspan = max(1e-9, view["x_max"] - view["x_min"])
    zspan = max(1e-9, view["z_max"] - view["z_min"])
    scale_x = w / xspan
    scale_z = h / zspan
    scale = min(scale_x, scale_z)
    draw_w = xspan * scale
    draw_h = zspan * scale
    offset_x = (w - draw_w) * 0.5
    offset_y = (h - draw_h) * 0.5
    return scale, offset_x, offset_y

def _to_px(x, z, view, w, h):
    scale, offx, offy = get_scene_transform(view, w, h)
    px = offx + (x - view["x_min"]) * scale
    py = offy + (view["z_max"] - z) * scale
    return px, py

# --------------------------
# Dibujo de zoom overlay
# --------------------------
def draw_zoom_rectangle(is_dragging, drag_start, drag_end):
    canvas = document.getElementById("ray-canvas")
    if not canvas or not is_dragging:
        return
    ctx = canvas.getContext("2d")
    x0, y0 = drag_start.get("x", 0), drag_start.get("y", 0)
    x1, y1 = drag_end.get("x", 0), drag_end.get("y", 0)
    left, top = min(x0, x1), min(y0, y1)
    w, h = abs(x1 - x0), abs(y1 - y0)
    ctx.save()
    ctx.fillStyle = "rgba(0, 120, 215, 0.15)"
    ctx.strokeStyle = "rgba(0, 120, 215, 0.9)"
    ctx.lineWidth = 1.0
    ctx.fillRect(left, top, w, h)
    ctx.strokeRect(left, top, w, h)
    ctx.restore()

# --------------------------
# Primitivas internas
# --------------------------
def _draw_segment(ctx, x1, z1, x2, z2, view, w, h):
    x1p, y1p = _to_px(x1, z1, view, w, h)
    x2p, y2p = _to_px(x2, z2, view, w, h)
    ctx.beginPath()
    ctx.moveTo(x1p, y1p)
    ctx.lineTo(x2p, y2p)
    ctx.stroke()

def _draw_circle(ctx, x, z, r, view, w, h):
    cx, cy = _to_px(x, z, view, w, h)
    scale, _, _ = get_scene_transform(view, w, h)
    rp = max(1.0, r * scale)
    ctx.beginPath()
    ctx.arc(cx, cy, rp, 0, 2 * math.pi)
    ctx.stroke()

def _clear_canvas(ctx, w, h):
    ctx.clearRect(0, 0, w, h)

# --------------------------
# Fondo: caja + tubo
# --------------------------
def _draw_background(ctx, geom, view, w, h):
    Wb = geom["W_base"]; Wc = geom["W_ceil"]; H = geom["H"]
    h_tube = geom["h_tube"]; r = geom["radio_m"]
    _clear_canvas(ctx, w, h)
    ctx.lineWidth = 1.0
    ctx.strokeStyle = "#888"
    _draw_segment(ctx, -Wb/2, 0,  Wb/2, 0, view, w, h)
    _draw_segment(ctx, -Wc/2, H,  Wc/2, H, view, w, h)
    _draw_segment(ctx, -Wb/2, 0, -Wc/2, H, view, w, h)
    _draw_segment(ctx,  Wb/2, 0,  Wc/2, H, view, w, h)
    ctx.strokeStyle = "#444"
    _draw_circle(ctx, 0.0, h_tube, r, view, w, h)

    # --- INICIO: NUEVO CÓDIGO PARA DIBUJAR REFLECTOR ---
    if geom.get("reflector_enabled"):
        h_refl_v = geom["h_reflector_vertex"]
        w_refl_a = geom["w_reflector_attach"]
        H = geom["H"]
        
        p_vertex = (0.0, h_refl_v)
        p_left = (-w_refl_a / 2.0, H)
        p_right = (w_refl_a / 2.0, H)

        ctx.lineWidth = 1.0
        ctx.strokeStyle = "#888" # Mismo estilo que las paredes
        _draw_segment(ctx, p_left[0], p_left[1], p_vertex[0], p_vertex[1], view, w, h)
        _draw_segment(ctx, p_vertex[0], p_vertex[1], p_right[0], p_right[1], view, w, h)
    # --- FIN: NUEVO CÓDIGO ---    
    
# --------------------------
# Capa de rayos (robusta)
# --------------------------
def _draw_rays(ctx, rays, geom, view, w, h):
    if not rays:
        return
    colors = {
        0: "#2ca02c",  # Verde oscuro para rayos directos (0 rebotes)
        1: "#1f77b4",  # Azul para 1 rebote
        "multi": "#888888" # Gris para 2 o más rebotes
    }
    for ray_path in rays:
        if not ray_path:
            continue
        last_segment = ray_path[-1]
        hit_base = abs(last_segment[3]) < 1e-6
        if hit_base:
            ctx.lineWidth = 2.0
        else:
            ctx.lineWidth = 1.0
        for item in ray_path:
            if isinstance(item, (list, tuple)) and len(item) == 5:
                try:
                    x1, z1, x2, z2, bounce_count = map(float, item)
                    if not (math.isfinite(x1) and math.isfinite(z1) and math.isfinite(x2) and math.isfinite(z2)):
                        continue
                except Exception:
                    continue
                if bounce_count == 0:
                    ctx.strokeStyle = colors.get(0, "#000")
                elif bounce_count == 1:
                    ctx.strokeStyle = colors.get(1, "#000")
                else:
                    ctx.strokeStyle = colors.get("multi", "#000")
                _draw_segment(ctx, x1, z1, x2, z2, view, w, h)

# --------------------------
# Eje Y NORMALIZADO (derecha)
# --------------------------
def _draw_normalized_y_axis(ctx, y_axis_data, geom, view, w, h):
    if not y_axis_data:
        return
    x_axis = geom["W_base"]/2
    zb = y_axis_data["Z_bottom"]; zt = y_axis_data["Z_top"]
    y_min_display = y_axis_data["y_min_axis"]; y_max_display = y_axis_data["y_max_axis"]
    unit_label = y_axis_data.get("unit_label", "Fracción")
    ctx.lineWidth = 1.0
    ctx.strokeStyle = "#000"
    _draw_segment(ctx, x_axis, zb, x_axis, zt, view, w, h)
    def Z_of(y_val):
        t = (y_val - y_min_display) / max(1e-6, (y_max_display - y_min_display))
        return zb + t * (zt - zb)
    major_tick_step = 0.2
    minor_tick_step = 0.1
    all_tick_values = set()
    current_y = math.floor(y_min_display / minor_tick_step) * minor_tick_step
    while current_y <= y_max_display + 1e-9:
        if current_y >= y_min_display - 1e-9:
             all_tick_values.add(round(current_y, 2))
        current_y += minor_tick_step
        current_y = round(current_y, 2)
    sorted_tick_values = sorted(list(all_tick_values))
    for yv in sorted_tick_values:
        z_val = Z_of(yv)
        is_major_tick = False
        if abs(round(yv / major_tick_step) * major_tick_step - yv) < 1e-9 or \
           abs(yv - 0.0) < 1e-9 or abs(yv - 1.0) < 1e-9:
           is_major_tick = True
        if is_major_tick:
            tick_length_m = 0.02 * geom["W_base"]
            label_font_size = "14px"
            label_offset_m = 0.025 * geom["W_base"]
            label_text = f"{yv:.1f}"
            ctx.lineWidth = 1.0
            ctx.strokeStyle = "#000"
            _draw_segment(ctx, x_axis, z_val, x_axis + tick_length_m, z_val, view, w, h)
            px, py = _to_px(x_axis + label_offset_m, z_val, view, w, h)
            ctx.fillStyle = "#000"
            ctx.font = label_font_size + " sans-serif"
            ctx.textAlign = "left"
            ctx.textBaseline = "middle"
            ctx.fillText(label_text, px, py)
        elif abs(round(yv / minor_tick_step) * minor_tick_step - yv) < 1e-9:
            tick_length_m = 0.01 * geom["W_base"]
            ctx.lineWidth = 0.5
            ctx.strokeStyle = "#444"
            _draw_segment(ctx, x_axis, z_val, x_axis + tick_length_m, z_val, view, w, h)
    label_text = unit_label
    px_label, py_label = _to_px(x_axis + 0.06 * geom["W_base"], zb + (zt - zb) / 2, view, w, h)
    ctx.save()
    ctx.translate(px_label, py_label)
    ctx.rotate(math.pi / 2)
    ctx.font = "14px sans-serif"
    ctx.fillStyle = "#000"
    ctx.textAlign = "center"
    ctx.textBaseline = "bottom"
    ctx.fillText(label_text, 0, 0)
    ctx.restore()

# --------------------------
# Eje Y REAL (izquierda)
# --------------------------
def _draw_real_y_axis(ctx, y_axis_data, geom, view, w, h):
    if not y_axis_data:
        return
    x_axis = -geom["W_base"]/2
    zb = y_axis_data["Z_bottom"]; zt = y_axis_data["Z_top"]
    y_min_display_real = y_axis_data["y_min_real"]; y_max_display_real = y_axis_data["y_max_real"]
    unit_label = y_axis_data.get("unit_label", "Irradiancia (u.a.)")
    
    ctx.lineWidth = 1.0
    ctx.strokeStyle = "#000"
    _draw_segment(ctx, x_axis, zb, x_axis, zt, view, w, h)
    
    # --- INICIO DE LA LÓGICA MEJORADA ---
    # Calcular dinámicamente el número de ticks en función del espacio disponible
    _, py_bottom = _to_px(x_axis, zb, view, w, h)
    _, py_top = _to_px(x_axis, zt, view, w, h)
    axis_height_px = abs(py_top - py_bottom)
    
    # Intentar tener un tick cada 40-50 píxeles, pero no menos de 2.
    desired_tick_spacing_px = 50 
    num_ticks_real = max(2, int(round(axis_height_px / desired_tick_spacing_px)))
    
    range_real = y_max_display_real - y_min_display_real
    if range_real <= 0: range_real = 1.0
    
    # Algoritmo para encontrar un "paso" legible (múltiplo de 1, 2, 5)
    rough_step = range_real / num_ticks_real
    exponent = math.floor(math.log10(rough_step))
    fractional = rough_step / (10**exponent)
    
    if fractional <= 1.0: nice_fractional = 1.0
    elif fractional <= 2.0: nice_fractional = 2.0
    elif fractional <= 5.0: nice_fractional = 5.0
    else: nice_fractional = 10.0
        
    major_step_real = nice_fractional * (10**exponent)
    # --- FIN DE LA LÓGICA MEJORADA ---

    def Z_of_real(y_val_real):
        t_real = (y_val_real - y_min_display_real) / max(1e-6, (y_max_display_real - y_min_display_real))
        return zb + t_real * (zt - zb)

    first_tick_real = math.floor(y_min_display_real / major_step_real) * major_step_real
    
    real_tick_values = []
    current_real_tick = first_tick_real
    while current_real_tick <= y_max_display_real + 1e-9:
        if current_real_tick >= y_min_display_real - 1e-9:
            real_tick_values.append(current_real_tick)
        current_real_tick += major_step_real

    if not real_tick_values: # Asegurar al menos un tick si el bucle falla
        if y_max_display_real > y_min_display_real:
             real_tick_values.append(y_min_display_real)
             real_tick_values.append(y_max_display_real)
        else:
             real_tick_values.append(y_min_display_real)

    for yv_real in real_tick_values:
        z_val_real = Z_of_real(yv_real)
        tick_length_m = 0.02 * geom["W_base"]
        label_font_size = "14px"
        label_offset_m = 0.025 * geom["W_base"]
        
        # Formato de la etiqueta
        if major_step_real >= 1:
             label_text = f"{yv_real:.1f}"
        elif major_step_real >= 0.1:
             label_text = f"{yv_real:.2f}"
        else:
             label_text = f"{yv_real:.3f}"
        if yv_real == 0.0 : label_text = "0.0"

        ctx.lineWidth = 1.0
        ctx.strokeStyle = "#000"
        _draw_segment(ctx, x_axis, z_val_real, x_axis - tick_length_m, z_val_real, view, w, h)
        px, py = _to_px(x_axis - label_offset_m, z_val_real, view, w, h)
        ctx.fillStyle = "#000"
        ctx.font = label_font_size + " sans-serif"
        ctx.textAlign = "right"
        ctx.textBaseline = "middle"
        ctx.fillText(label_text, px, py)
        
    # Etiqueta del eje
    label_text = unit_label
    px_label, py_label = _to_px(x_axis - 0.08 * geom["W_base"], zb + (zt - zb) / 2, view, w, h) # Ajuste de posición
    ctx.save()
    ctx.translate(px_label, py_label)
    ctx.rotate(-math.pi / 2)
    ctx.font = "14px sans-serif"
    ctx.fillStyle = "#000"
    ctx.textAlign = "center"
    ctx.textBaseline = "bottom"
    ctx.fillText(label_text, 0, 0)
    ctx.restore()

# --------------------------
# Eje X en centímetros (en la base)
# --------------------------
def _draw_x_axis_cm(ctx, geom, view, w, h):
    Wb = geom["W_base"]
    axis_z = 0.0
    cm_step = 0.01
    start_x_m = -Wb / 2
    end_x_m = Wb / 2
    ctx.lineWidth = 1.0
    ctx.strokeStyle = "#000"
    _draw_segment(ctx, max(start_x_m, view["x_min"]), axis_z, min(end_x_m, view["x_max"]), axis_z, view, w, h)
    ctx.fillStyle = "#000"
    ctx.font = "10px sans-serif"
    ctx.textAlign = "center"
    ctx.textBaseline = "top"
    current_x = math.ceil(start_x_m / cm_step) * cm_step
    while current_x <= end_x_m:
        if current_x >= view["x_min"] and current_x <= view["x_max"]:
            tick_half_height_m_short = 0.0025
            tick_half_height_m_long = 0.005
            label_offset_m = 0.008
            is_main_tick = abs(round(current_x * 100) % 5) < 0.01
            if is_main_tick:
                tick_h = tick_half_height_m_long * 2
                z_start_tick = axis_z - tick_half_height_m_long
                z_end_tick = axis_z + tick_half_height_m_long
                label_text = f"{int(round(current_x * 100))}"
                px, py = _to_px(current_x, axis_z - label_offset_m, view, w, h)
                ctx.fillText(label_text, px, py)
            else:
                tick_h = tick_half_height_m_short * 2
                z_start_tick = axis_z - tick_half_height_m_short
                z_end_tick = axis_z + tick_half_height_m_short
            ctx.lineWidth = 1.0
            ctx.strokeStyle = "#000"
            _draw_segment(ctx, current_x, z_start_tick, current_x, z_end_tick, view, w, h)
        current_x += cm_step
    x_label_pos_m = end_x_m + 0.02 * geom["W_base"]
    y_label_pos_m = axis_z - 0.015 * geom["H"]
    px_label, py_label = _to_px(x_label_pos_m, y_label_pos_m, view, w, h)
    ctx.fillStyle = "#000"
    ctx.font = "10px sans-serif"
    ctx.textAlign = "center"
    ctx.textBaseline = "top"
    ctx.fillText("(cm)", px_label, py_label)

# --------------------------
# Curvas (total, directa, first, multi)
# --------------------------
def _draw_curves(ctx, curves, geom, view, w, h):
    if not curves:
        return
    styles = {
        "total":  {"color": "#d62728", "label": "Total"},
        "direct": {"color": "#2ca02c", "label": "Directa"},
        "first":  {"color": "#ff7f0e", "label": "1 Rebote"},
        "multi":  {"color": "#1f77b4", "label": "Multi Rebote"},
    }
    ctx.lineWidth = 1.5
    for key, pts in curves.items():
        if key in styles:
            col = styles.get(key)["color"]
            ctx.strokeStyle = col
            started = False
            ctx.beginPath()
            for (x, z) in pts:
                try:
                    x = float(x); z = float(z)
                    if not (math.isfinite(x) and math.isfinite(z)):
                        started = False
                        break
                except Exception:
                    started = False
                    break
                px, py = _to_px(x, z, view, w, h)
                if not started:
                    ctx.moveTo(px, py); started = True
                else:
                    ctx.lineTo(px, py)
            if started:
                ctx.stroke()
# --------------------------
# Marcador de inflexión (curva total)
# --------------------------
def _mark_inflection(ctx, curves, inflection_x, geom, view, w, h):
    if inflection_x is None:
        return
    pts = curves.get("total", [])
    if not pts:
        return
    def nearest_idx(x_target):
        return min(range(len(pts)), key=lambda i: abs(float(pts[i][0]) - x_target))
    for sx in (+1, -1):
        x_target = sx * abs(float(inflection_x))
        i = nearest_idx(x_target)
        x2, z2 = float(pts[i][0]), float(pts[i][1])
        scale, _, _ = get_scene_transform(view, w, h)
        rpx = max(3, int(0.004 * geom["W_base"] * scale))
        cx, cy = _to_px(x2, z2, view, w, h)
        ctx.save()
        ctx.strokeStyle = "#d62728"
        ctx.lineWidth = 2.0
        ctx.beginPath()
        ctx.arc(cx, cy, rpx, 0, 2 * math.pi)
        ctx.stroke()
        ctx.restore()

# --------------------------
# Leyenda de curvas
# --------------------------
def _draw_curve_legend(ctx, curves, geom, w_logical, h_logical):
    if not curves:
        return
    styles = {
        "total":  {"color": "#d62728", "label": "Total"},
        "direct": {"color": "#2ca02c", "label": "Directa"},
        "first":  {"color": "#ff7f0e", "label": "1 Rebote"},
        "multi":  {"color": "#1f77b4", "label": "Multi Rebote"},
    }
    margin_x_px = 15
    margin_y_px = 15
    line_height_px = 20
    box_size_px = 14
    text_offset_x_px = box_size_px + 8
    max_text_width = 0
    ctx.font = "14px sans-serif"
    for key in ["total", "direct", "first", "multi"]:
        if key in curves and curves[key]:
            max_text_width = max(max_text_width, ctx.measureText(styles[key]["label"]).width)
    legend_width = max_text_width + text_offset_x_px + box_size_px + 20
    legend_height = len(styles) * line_height_px + 10
    legend_x_px = w_logical - margin_x_px - legend_width
    legend_y_px = margin_y_px
    ctx.save()
    ctx.fillStyle = "rgba(255, 255, 255, 0.85)"
    ctx.fillRect(legend_x_px, legend_y_px, legend_width, legend_height)
    ctx.strokeStyle = "#ccc"
    ctx.lineWidth = 1.0
    ctx.strokeRect(legend_x_px, legend_y_px, legend_width, legend_height)
    ctx.font = "14px sans-serif"
    ctx.textBaseline = "middle"
    ctx.textAlign = "left"
    current_y_px = legend_y_px + line_height_px / 2 + 5
    for key in ["total", "direct", "first", "multi"]:
        if key in curves and curves[key]:
            style = styles.get(key)
            if style:
                color = style["color"]
                label = style["label"]
                color_box_x_px = legend_x_px + margin_x_px
                ctx.fillStyle = color
                ctx.fillRect(color_box_x_px, current_y_px - box_size_px / 2, box_size_px, box_size_px)
                ctx.fillStyle = "#000"
                ctx.fillText(label, color_box_x_px + text_offset_x_px, current_y_px)
                current_y_px += line_height_px
    ctx.restore()

# --------------------------
# Dibujo principal
# --------------------------
def draw_main_scene(rays, curves, view, geom):
    canvas = document.getElementById("ray-canvas")
    if not canvas:
        return
    device_pixel_ratio = window.devicePixelRatio or 1
    css_width = int(canvas.clientWidth)
    css_height = int(canvas.clientHeight)
    if css_width == 0 or css_height == 0:
        return
    if canvas.width != css_width * device_pixel_ratio or canvas.height != css_height * device_pixel_ratio:
        canvas.width = css_width * device_pixel_ratio
        canvas.height = css_height * device_pixel_ratio
    ctx = canvas.getContext("2d")
    ctx.setTransform(device_pixel_ratio, 0, 0, device_pixel_ratio, 0, 0)
    w_logical = css_width
    h_logical = css_height
    _draw_background(ctx, geom, view, w_logical, h_logical)
    _draw_x_axis_cm(ctx, geom, view, w_logical, h_logical)
    mode = geom.get("render_mode", "RAYS")
    if mode == "RAYS":
        _draw_rays(ctx, rays, geom, view, w_logical, h_logical)
    else: # mode == "CURVES"
        _draw_curves(ctx, curves, geom, view, w_logical, h_logical)
        if "y_axis_real" in geom and geom["y_axis_real"]:
            _draw_real_y_axis(ctx, geom["y_axis_real"], geom, view, w_logical, h_logical)
        if "y_axis_normalized" in geom and geom["y_axis_normalized"]:
            _draw_normalized_y_axis(ctx, geom["y_axis_normalized"], geom, view, w_logical, h_logical)
        _draw_curve_legend(ctx, curves, geom, w_logical, h_logical)
        if "inflection_x" in geom and geom["inflection_x"] is not None:
            _mark_inflection(ctx, curves, geom["inflection_x"], geom, view, w_logical, h_logical)
            
            