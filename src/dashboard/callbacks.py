"""Dash callbacks for the Nested RF Stimulus Ephys Dashboard."""

from __future__ import annotations

from collections import defaultdict

from dash import Input, Output, State, callback_context, html, no_update
import dash_bootstrap_components as dbc
from dash import dash_table

from dashboard.data_loader import DataStore

# Placeholder shown when an image is not available
_NOT_AVAILABLE = html.Div(
    "Not available",
    className="text-muted text-center py-3",
    style={"fontSize": "0.95rem"},
)

# Image type ↔ gallery tab mapping (updated to new image keys)
_GALLERY_TABS = {
    "tab-flash": "flash_rf_6px",
    "tab-bar-sweep": "bar_polar",
    "tab-bar-flash": "bar_flash_slow",
}


def _build_metrics_table(cell: dict) -> list:
    """Build a compact metrics card from the cell's 'metrics' dict.

    Returns a list of Dash components (dbc.Card with inner tables).
    If no metrics are available, returns a placeholder message.
    """
    metrics = cell.get("metrics", {})
    if not metrics:
        return html.Div(
            "No metrics available for this cell.",
            className="text-muted text-center py-2",
            style={"fontSize": "0.9rem"},
        )

    sections = []

    # --- Row 1: Metadata ---
    meta_items = []
    for key, label in [
        ("strain", "Strain"),
        ("frame", "Frame"),
        ("side", "Side"),
        ("age", "Age"),
        ("preferred_contrast", "Contrast"),
        ("resultant_angle_deg", "Pref. Dir (°)"),
    ]:
        val = metrics.get(key, cell.get(key, ""))
        if val != "" and val is not None:
            if isinstance(val, float):
                val = f"{val:.1f}"
            meta_items.append(html.Span(
                [html.Strong(f"{label}: "), str(val)],
                className="me-3",
                style={"fontSize": "0.85rem"},
            ))
    if meta_items:
        sections.append(html.Div(meta_items, className="mb-2"))

    # --- Row 2: Bar sweep metrics table (slow / fast / vfast) ---
    bar_speeds = []
    for speed_key, speed_label in [
        ("bar_slow", "Slow"), ("bar_fast", "Fast"), ("bar_vfast", "V.Fast"),
    ]:
        speed = metrics.get(speed_key)
        if speed and isinstance(speed, dict):
            bar_speeds.append({
                "Speed": speed_label,
                "Mag.": _fmt(speed.get("magnitude")),
                "DSI (vec)": _fmt(speed.get("DSI_vector")),
                "DSI (pd/nd)": _fmt(speed.get("DSI_pdnd")),
                "FWHM (°)": _fmt(speed.get("fwhm"), decimals=0),
                "CV": _fmt(speed.get("cv")),
            })
    if bar_speeds:
        sections.append(html.Div([
            html.Strong("Bar Sweep", style={"fontSize": "0.8rem"}),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in bar_speeds[0].keys()],
                data=bar_speeds,
                style_header={
                    "fontWeight": "bold",
                    "backgroundColor": "#f8f9fa",
                    "fontSize": "0.78rem",
                    "padding": "4px 6px",
                    "border": "1px solid #dee2e6",
                },
                style_cell={
                    "fontSize": "0.78rem",
                    "padding": "3px 6px",
                    "textAlign": "center",
                    "border": "1px solid #dee2e6",
                    "fontFamily": "inherit",
                },
                style_table={"marginBottom": "8px"},
            ),
        ], className="mb-2"))

    # --- Row 3: RF fit metrics table (4px / 6px) ---
    rf_rows = []
    for rf_key, rf_label in [("rf_4px", "4 px"), ("rf_6px", "6 px")]:
        rf = metrics.get(rf_key)
        if rf and isinstance(rf, dict):
            rf_rows.append({
                "Size": rf_label,
                "R² (exc)": _fmt(rf.get("R_squared"), 3),
                "R² (inh)": _fmt(rf.get("R_squaredi"), 3),
                "σx (exc)": _fmt(rf.get("sigma_x_exc")),
                "σy (exc)": _fmt(rf.get("sigma_y_exc")),
                "σx (inh)": _fmt(rf.get("sigma_x_inh")),
                "σy (inh)": _fmt(rf.get("sigma_y_inh")),
            })
    if rf_rows:
        sections.append(html.Div([
            html.Strong("RF Gaussian Fit", style={"fontSize": "0.8rem"}),
            dash_table.DataTable(
                columns=[{"name": c, "id": c} for c in rf_rows[0].keys()],
                data=rf_rows,
                style_header={
                    "fontWeight": "bold",
                    "backgroundColor": "#f8f9fa",
                    "fontSize": "0.78rem",
                    "padding": "4px 6px",
                    "border": "1px solid #dee2e6",
                },
                style_cell={
                    "fontSize": "0.78rem",
                    "padding": "3px 6px",
                    "textAlign": "center",
                    "border": "1px solid #dee2e6",
                    "fontFamily": "inherit",
                },
                style_table={"marginBottom": "4px"},
            ),
        ]))

    if not sections:
        return html.Div(
            "No metrics available for this cell.",
            className="text-muted text-center py-2",
            style={"fontSize": "0.9rem"},
        )

    return dbc.Card(
        dbc.CardBody(sections, style={"padding": "10px 14px"}),
        className="mb-3",
    )


def _fmt(value, decimals: int = 2) -> str:
    """Format a numeric value for table display, handling None/NaN."""
    if value is None:
        return "—"
    try:
        v = float(value)
        if v != v:  # NaN check
            return "—"
        if decimals == 0:
            return f"{v:.0f}"
        return f"{v:.{decimals}f}"
    except (TypeError, ValueError):
        return str(value)


def register_callbacks(app, data_store: DataStore) -> None:
    """Register all Dash callbacks."""

    # ==================================================================
    # 1. Data path validation + strain dropdown
    # ==================================================================
    @app.callback(
        [
            Output("data-status", "children"),
            Output("strain-dropdown", "options"),
            Output("strain-dropdown", "value"),
        ],
        Input("images-path-input", "value"),
    )
    def update_data_path(path):
        if not path:
            return "No path specified", [], None

        valid = data_store.reload(path)
        if not valid:
            return (
                "cell_index.json not found in this directory",
                [],
                None,
            )

        strains = data_store.get_strains()
        if not strains:
            return "No cells found in index", [], None

        n = len(data_store.load_cell_index())
        options = [{"label": s, "value": s} for s in strains]
        return f"{n} cells loaded", options, strains[0]

    # ==================================================================
    # 2. Cascading contrast dropdown
    # ==================================================================
    @app.callback(
        [
            Output("contrast-dropdown", "options"),
            Output("contrast-dropdown", "value"),
        ],
        Input("strain-dropdown", "value"),
    )
    def update_contrast_options(strain):
        if not strain:
            return [], None
        contrasts = data_store.get_contrasts_for_strain(strain)
        options = [{"label": c, "value": c} for c in contrasts]
        value = contrasts[0] if contrasts else None
        return options, value

    # ==================================================================
    # 3. Cascading cell dropdown
    # ==================================================================
    @app.callback(
        [
            Output("cell-dropdown", "options"),
            Output("cell-dropdown", "value"),
        ],
        [
            Input("strain-dropdown", "value"),
            Input("contrast-dropdown", "value"),
        ],
    )
    def update_cell_options(strain, contrast):
        if not strain or not contrast:
            return [], None
        cells = data_store.get_cells_for_filter(strain, contrast)
        options = [
            {"label": f"{c['date']}  {c['time']}", "value": c["cell_id"]}
            for c in cells
        ]
        value = cells[0]["cell_id"] if cells else None
        return options, value

    # ==================================================================
    # 4. Metadata table
    # ==================================================================
    @app.callback(
        [
            Output("metadata-table", "columns"),
            Output("metadata-table", "data"),
        ],
        Input("main-tabs", "active_tab"),
        Input("images-path-input", "value"),
    )
    def update_metadata_table(active_tab, _path):
        if active_tab != "tab-metadata":
            return no_update, no_update

        if not data_store.is_valid:
            return [], []

        summary = data_store.get_metadata_summary()
        if summary.empty:
            return [], []

        columns = [{"name": c, "id": c} for c in summary.columns]
        data = summary.to_dict("records")
        return columns, data

    # ==================================================================
    # 5. Per Cell — compact overview (16 outputs)
    # ==================================================================
    @app.callback(
        [
            Output("cell-title", "children"),                       # 0
            Output("cell-notes", "children"),                       # 1
            Output("metrics-table-container", "children"),          # 2
            Output("img-gridplot-1-container", "children"),         # 3
            Output("img-gridplot-2-container", "children"),         # 4
            Output("img-gridplot-3-container", "children"),         # 5
            Output("img-gridplot-4-container", "children"),         # 6
            Output("img-bar-polar-container", "children"),          # 7
            Output("img-polar-arrow-container", "children"),        # 8
            Output("img-flash-rf-4px-container", "children"),       # 9
            Output("img-flash-heatmap-4px-container", "children"),  # 10
            Output("img-flash-rf-6px-container", "children"),       # 11
            Output("img-flash-heatmap-6px-container", "children"),  # 12
            Output("img-gaussian-rf-container", "children"),        # 13
            Output("img-bar-flash-slow-container", "children"),     # 14
            Output("img-bar-flash-fast-container", "children"),     # 15
        ],
        Input("cell-dropdown", "value"),
        Input("main-tabs", "active_tab"),
    )
    def update_per_cell_images(cell_id, active_tab):
        n_outputs = 16
        if active_tab != "tab-per-cell" or not cell_id:
            return (no_update,) * n_outputs

        cell = data_store.get_cell_by_id(cell_id)
        if cell is None:
            return ("Cell not found", "", _NOT_AVAILABLE, *([_NOT_AVAILABLE] * 13))

        title = cell.get("display_label", cell_id)

        # Build notes string
        notes_parts = []
        ns = cell.get("notes_start", "")
        ne = cell.get("notes_end", "")
        if ns:
            notes_parts.append(ns)
        if ne:
            notes_parts.append(ne)
        notes = " | ".join(notes_parts) if notes_parts else ""

        # Metrics table
        metrics_card = _build_metrics_table(cell)

        def _img_or_na(image_type, max_width="100%"):
            """Return an html.Img component or 'Not available' placeholder."""
            url = data_store.get_image_url(cell, image_type)
            if url:
                return html.Img(
                    src=url,
                    style={"maxWidth": max_width, "height": "auto"},
                )
            return _NOT_AVAILABLE

        # GridPlot images (4 separate sub-figures)
        gp1 = _img_or_na("gridplot_1")
        gp2 = _img_or_na("gridplot_2")
        gp3 = _img_or_na("gridplot_3")
        gp4 = _img_or_na("gridplot_4")

        # Analysis plots (compact layout)
        bar_polar = _img_or_na("bar_polar")
        polar_arrow = _img_or_na("polar_arrow")
        gaussian_rf = _img_or_na("gaussian_rf")
        flash_rf_4px = _img_or_na("flash_rf_4px")
        flash_heatmap_4px = _img_or_na("flash_heatmap_4px")
        flash_rf_6px = _img_or_na("flash_rf_6px")
        flash_heatmap_6px = _img_or_na("flash_heatmap_6px")
        bar_flash_slow = _img_or_na("bar_flash_slow")
        bar_flash_fast = _img_or_na("bar_flash_fast")

        return (
            title, notes, metrics_card,
            gp1, gp2, gp3, gp4,
            bar_polar, polar_arrow,
            flash_rf_4px, flash_heatmap_4px,
            flash_rf_6px, flash_heatmap_6px,
            gaussian_rf,
            bar_flash_slow, bar_flash_fast,
        )

    # ==================================================================
    # 6. Gallery tab callbacks (one set per gallery tab)
    # ==================================================================
    for tab_id, image_type in _GALLERY_TABS.items():
        _register_gallery_callbacks(app, data_store, tab_id, image_type)


def _build_grouped_checklist(cells: list[dict], tab_id: str) -> list:
    """Build a list of HTML elements with cells grouped by strain and PC.

    Returns a list of Dash components with group headings and checkboxes.
    Each checkbox value is the cell_id; labels show date and time only.
    """
    # Group cells: strain → preferred_contrast → list of cells
    groups: dict[str, dict[str, list[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for c in cells:
        groups[c["strain"]][c["preferred_contrast"]].append(c)

    elements = []
    for strain in sorted(groups.keys()):
        # Strain heading
        # Shorten strain name for display (e.g. "42F06_T4T5_control" → "control")
        short_strain = strain.split("_")[-1] if "_" in strain else strain
        elements.append(
            html.Div(
                short_strain,
                className="fw-bold mt-2 mb-1",
                style={"fontSize": "0.9rem", "textTransform": "uppercase"},
            )
        )
        for pc in sorted(groups[strain].keys()):
            # PC sub-heading
            elements.append(
                html.Div(
                    pc,
                    className="fw-semibold ms-1 mb-1",
                    style={
                        "fontSize": "0.8rem",
                        "color": "#e74c3c" if pc == "ON" else "#3498db",
                    },
                )
            )
            # Checkboxes for this group
            group_cells = sorted(
                groups[strain][pc], key=lambda c: c["cell_id"]
            )
            options = [
                {
                    "label": f"{c['date']}  {c['time']}",
                    "value": c["cell_id"],
                }
                for c in group_cells
            ]
            elements.append(
                dbc.Checklist(
                    id={
                        "type": f"{tab_id}-group-check",
                        "index": f"{strain}-{pc}",
                    },
                    options=options,
                    value=[],
                    className="ms-2 mb-1",
                    style={"fontSize": "0.85rem"},
                )
            )

    return elements


def _register_gallery_callbacks(
    app, data_store: DataStore, tab_id: str, image_type: str
) -> None:
    """Register the cell-list and gallery callbacks for one gallery tab."""
    from dash import ALL

    # --- Populate grouped cell checklist when tab becomes active ---
    @app.callback(
        [
            Output(f"{tab_id}-grouped-cells", "children"),
            Output(f"{tab_id}-cell-list", "options"),
        ],
        Input("main-tabs", "active_tab"),
        Input("images-path-input", "value"),
    )
    def update_cell_list(active_tab, _path, _tab_id=tab_id, _img_type=image_type):
        if active_tab != _tab_id:
            return no_update, no_update

        if not data_store.is_valid:
            return [], []

        cells = data_store.get_cells_with_image(_img_type)
        grouped = _build_grouped_checklist(cells, _tab_id)
        # Also build flat options for the hidden checklist (used by show-all)
        flat_options = [
            {"label": c["display_label"], "value": c["cell_id"]}
            for c in cells
        ]
        return grouped, flat_options

    # --- "Show all" toggles all group checkboxes ---
    @app.callback(
        Output({"type": f"{tab_id}-group-check", "index": ALL}, "value"),
        Input(f"{tab_id}-show-all", "value"),
        State({"type": f"{tab_id}-group-check", "index": ALL}, "options"),
    )
    def toggle_show_all(show_all, all_options, _tab_id=tab_id):
        if not all_options:
            return no_update
        if "all" in (show_all or []):
            return [
                [opt["value"] for opt in opts]
                for opts in all_options
            ]
        return [[] for _ in all_options]

    # --- Sync group checkboxes → hidden flat cell-list ---
    @app.callback(
        Output(f"{tab_id}-cell-list", "value"),
        Input({"type": f"{tab_id}-group-check", "index": ALL}, "value"),
    )
    def sync_to_flat(group_values, _tab_id=tab_id):
        if not group_values:
            return []
        # Flatten all group selections into a single list
        return [v for vals in group_values for v in (vals or [])]

    # --- Render gallery images ---
    @app.callback(
        Output(f"{tab_id}-gallery", "children"),
        Input(f"{tab_id}-cell-list", "value"),
        Input("main-tabs", "active_tab"),
    )
    def update_gallery(selected, active_tab, _tab_id=tab_id, _img_type=image_type):
        if active_tab != _tab_id:
            return no_update

        if not selected:
            return html.Div(
                "Select cells from the list to view their plots.",
                className="text-muted text-center py-5",
            )

        cards = []
        for cell_id in selected:
            cell = data_store.get_cell_by_id(cell_id)
            if cell is None:
                continue

            # Use thumbnail for gallery, full-size as link target
            thumb_url = data_store.get_image_url(cell, _img_type, thumb=True)
            full_url = data_store.get_image_url(cell, _img_type, thumb=False)

            if thumb_url:
                img = html.A(
                    html.Img(
                        src=thumb_url,
                        style={"maxWidth": "100%", "height": "auto"},
                    ),
                    href=full_url,
                    target="_blank",
                )
            else:
                img = _NOT_AVAILABLE

            # For Square Flashes tab, also show heatmap below timeseries
            body_children = [img]
            if _img_type == "flash_rf_6px":
                hm_thumb = data_store.get_image_url(cell, "flash_heatmap_6px", thumb=True)
                hm_full = data_store.get_image_url(cell, "flash_heatmap_6px", thumb=False)
                if hm_thumb:
                    body_children.append(
                        html.Hr(style={"margin": "8px 0"}),
                    )
                    body_children.append(
                        html.A(
                            html.Img(
                                src=hm_thumb,
                                style={"maxWidth": "100%", "height": "auto"},
                            ),
                            href=hm_full,
                            target="_blank",
                        )
                    )

            # Short label: date time (strain-pc)
            short_strain = cell.get("strain", "").split("_")[-1]
            label = f"{cell['date']} {cell['time']} ({short_strain}-{cell.get('preferred_contrast', '')})"

            cards.append(
                dbc.Col(
                    dbc.Card(
                        [
                            dbc.CardHeader(
                                label,
                                style={"fontSize": "0.85rem"},
                            ),
                            dbc.CardBody(body_children),
                        ],
                        className="mb-3",
                    ),
                    xs=12,
                    md=6,
                    xl=4,
                )
            )

        if not cards:
            return html.Div(
                "No images available for the selected cells.",
                className="text-muted text-center py-5",
            )

        return dbc.Row(cards)
