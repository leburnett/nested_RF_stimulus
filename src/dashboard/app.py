"""Nested RF Stimulus — Electrophysiology Dashboard.

Usage:
    python -m dashboard
"""

import flask
import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from dashboard.callbacks import register_callbacks
from dashboard.constants import DEFAULT_IMAGES_DIR, IMAGE_CACHE_MAX_AGE, PORT
from dashboard.data_loader import DataStore

# ------------------------------------------------------------------
# Initialise app
# ------------------------------------------------------------------
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,
)
app.title = "Nested RF Stimulus — Ephys Dashboard"

data_store = DataStore(DEFAULT_IMAGES_DIR)

# ------------------------------------------------------------------
# Serve preprocessed images with caching
# ------------------------------------------------------------------

@app.server.route("/images/<path:filename>")
def serve_image(filename):
    resp = flask.send_from_directory(str(data_store.images_dir), filename)
    resp.cache_control.max_age = IMAGE_CACHE_MAX_AGE
    return resp


# ------------------------------------------------------------------
# Sidebar
# ------------------------------------------------------------------
sidebar = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Label("Images Directory", className="fw-bold"),
                dbc.Input(
                    id="images-path-input",
                    type="text",
                    value=DEFAULT_IMAGES_DIR,
                    debounce=True,
                    className="mb-1",
                    size="sm",
                ),
                html.Small(
                    id="data-status", className="text-muted d-block mb-3"
                ),
                html.Hr(),
                dbc.Label("Strain", className="fw-bold"),
                dcc.Dropdown(
                    id="strain-dropdown", className="mb-3", clearable=False
                ),
                dbc.Label("Preferred Contrast", className="fw-bold"),
                dcc.Dropdown(
                    id="contrast-dropdown", className="mb-3", clearable=False
                ),
                dbc.Label("Cell", className="fw-bold"),
                dcc.Dropdown(
                    id="cell-dropdown", className="mb-3", clearable=False
                ),
            ]
        ),
    ],
    className="sticky-top",
    style={"top": "10px"},
)

# ------------------------------------------------------------------
# Tab 1: Metadata
# ------------------------------------------------------------------
metadata_tab = dbc.Tab(
    label="Metadata",
    tab_id="tab-metadata",
    children=[
        html.Div(
            [
                html.H5("Dataset Overview", className="mb-3"),
                dash_table.DataTable(
                    id="metadata-table",
                    columns=[],
                    data=[],
                    sort_action="native",
                    style_header={
                        "fontWeight": "bold",
                        "backgroundColor": "#f8f9fa",
                        "border": "1px solid #dee2e6",
                    },
                    style_data_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "#f8f9fa",
                        }
                    ],
                    style_cell={
                        "padding": "8px 12px",
                        "fontSize": "14px",
                        "fontFamily": "inherit",
                        "border": "1px solid #dee2e6",
                        "textAlign": "left",
                    },
                    style_table={"marginBottom": "24px"},
                ),
            ],
            className="p-3",
        ),
    ],
)

# ------------------------------------------------------------------
# Tab 2: Per Cell — Compact overview layout
# ------------------------------------------------------------------
per_cell_tab = dbc.Tab(
    label="Per Cell",
    tab_id="tab-per-cell",
    children=[
        html.Div(
            [
                html.H5(id="cell-title", className="mb-1"),
                html.Div(id="cell-notes", className="text-muted mb-3",
                         style={"fontSize": "0.9rem"}),

                # --- Metrics Table (top of page) ---
                html.Div(id="metrics-table-container", className="mb-3"),

                # --- Bar Analysis (3 cols) ---
                html.H6("Bar Sweep Analysis", className="mt-3 mb-2"),
                html.Hr(style={"marginTop": "0", "marginBottom": "8px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-bar-polar-container",
                                     className="text-center"),
                            xs=12, md=4, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-polar-arrow-container",
                                     className="text-center"),
                            xs=12, md=4, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-gaussian-rf-container",
                                     className="text-center"),
                            xs=12, md=4, className="mb-2",
                        ),
                    ],
                ),

                # --- Flash RF 4px (2 cols) ---
                html.H6("Flash RF — 4 px", className="mt-3 mb-2"),
                html.Hr(style={"marginTop": "0", "marginBottom": "8px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-flash-rf-4px-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-flash-heatmap-4px-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                    ],
                ),

                # --- Flash RF 6px (2 cols) ---
                html.H6("Flash RF — 6 px", className="mt-3 mb-2"),
                html.Hr(style={"marginTop": "0", "marginBottom": "8px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-flash-rf-6px-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-flash-heatmap-6px-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                    ],
                ),

                # --- Bar Flash (2 cols) ---
                html.H6("Bar Flash", className="mt-3 mb-2"),
                html.Hr(style={"marginTop": "0", "marginBottom": "8px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-bar-flash-slow-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-bar-flash-fast-container",
                                     className="text-center"),
                            xs=12, md=6, className="mb-2",
                        ),
                    ],
                ),

                # --- Grid Plots (Protocol 1) at bottom ---
                html.H6("Grid Plots (Protocol 1)", className="mt-4 mb-2"),
                html.Hr(style={"marginTop": "0", "marginBottom": "8px"}),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-gridplot-1-container"),
                            width=6, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-gridplot-2-container"),
                            width=6, className="mb-2",
                        ),
                    ],
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(id="img-gridplot-3-container"),
                            width=6, className="mb-2",
                        ),
                        dbc.Col(
                            html.Div(id="img-gridplot-4-container"),
                            width=6, className="mb-2",
                        ),
                    ],
                ),
            ],
            className="p-3",
        ),
    ],
)


# ------------------------------------------------------------------
# Tabs 3-5: Gallery tabs
# ------------------------------------------------------------------
def _make_gallery_tab(tab_id: str, label: str) -> dbc.Tab:
    """Create a gallery tab with grouped cell checkboxes and image gallery."""
    return dbc.Tab(
        label=label,
        tab_id=tab_id,
        children=[
            html.Div(
                [
                    dbc.Row(
                        [
                            # Left: grouped cell checklist
                            dbc.Col(
                                [
                                    dbc.Checklist(
                                        id=f"{tab_id}-show-all",
                                        options=[
                                            {
                                                "label": "Show all",
                                                "value": "all",
                                            }
                                        ],
                                        value=[],
                                        className="fw-bold mb-2",
                                    ),
                                    # Hidden checklist stores all selected values
                                    dbc.Checklist(
                                        id=f"{tab_id}-cell-list",
                                        options=[],
                                        value=[],
                                        style={"display": "none"},
                                    ),
                                    # Grouped checkboxes rendered dynamically
                                    html.Div(
                                        id=f"{tab_id}-grouped-cells",
                                        style={
                                            "maxHeight": "75vh",
                                            "overflowY": "auto",
                                        },
                                    ),
                                ],
                                width=3,
                            ),
                            # Right: image gallery
                            dbc.Col(
                                dcc.Loading(
                                    html.Div(id=f"{tab_id}-gallery"),
                                    type="circle",
                                ),
                                width=9,
                            ),
                        ],
                    ),
                ],
                className="p-3",
            ),
        ],
    )


flash_tab = _make_gallery_tab("tab-flash", "Square Flashes")
bar_sweep_tab = _make_gallery_tab("tab-bar-sweep", "Bar Sweeps")
bar_flash_tab = _make_gallery_tab("tab-bar-flash", "Bar Flashes")


# ------------------------------------------------------------------
# Main layout
# ------------------------------------------------------------------
app.layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    html.H3(
                        "Nested RF Stimulus — Ephys Dashboard",
                        className="text-primary my-3",
                    ),
                ),
            ],
        ),
        dbc.Row(
            [
                dbc.Col(sidebar, width=3),
                dbc.Col(
                    dbc.Tabs(
                        [
                            metadata_tab,
                            per_cell_tab,
                            flash_tab,
                            bar_sweep_tab,
                            bar_flash_tab,
                        ],
                        id="main-tabs",
                        active_tab="tab-metadata",
                    ),
                    width=9,
                ),
            ],
        ),
    ],
    fluid=True,
)

# Register callbacks
register_callbacks(app, data_store)

if __name__ == "__main__":
    app.run(debug=True, port=PORT)
