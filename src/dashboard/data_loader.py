"""Load cell index and experiment metadata for the dashboard."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from dashboard.constants import (
    DEFAULT_IMAGES_DIR,
    EXCLUDED_STRAINS,
    EXPERIMENT_LOG_PATH,
    PROJECT_TYPES,
)


class DataStore:
    """Manages cell index data and image path resolution.

    Loads ``cell_index.json`` once and caches it in memory for fast
    repeated access during dashboard interactions.
    """

    def __init__(self, images_dir: str = DEFAULT_IMAGES_DIR) -> None:
        self.images_dir = Path(images_dir)
        self._cell_index: list[dict] | None = None
        self._metadata_df: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    @property
    def is_valid(self) -> bool:
        return (self.images_dir / "cell_index.json").exists()

    def reload(self, images_dir: str) -> bool:
        """Switch to a new images directory. Returns True if valid."""
        self.images_dir = Path(images_dir)
        self._cell_index = None
        self._metadata_df = None
        return self.is_valid

    # ------------------------------------------------------------------
    # Cell index
    # ------------------------------------------------------------------

    def load_cell_index(self) -> list[dict]:
        if self._cell_index is None:
            index_path = self.images_dir / "cell_index.json"
            with open(index_path) as f:
                data = json.load(f)
            self._cell_index = data["cells"]
        return self._cell_index

    # ------------------------------------------------------------------
    # Experiment log (Excel)
    # ------------------------------------------------------------------

    def load_experiment_log(self) -> pd.DataFrame:
        if self._metadata_df is None:
            log_path = Path(EXPERIMENT_LOG_PATH)
            if not log_path.exists():
                return pd.DataFrame()
            df = pd.read_excel(log_path, sheet_name="Form responses 1")
            # Filter to project types in scope
            mask = df["project type"].str.contains(
                "|".join(PROJECT_TYPES), na=False
            )
            # Exclude test strains
            for strain in EXCLUDED_STRAINS:
                mask = mask & ~df["strain"].str.contains(
                    strain, case=False, na=False
                )
            self._metadata_df = df[mask].reset_index(drop=True)
        return self._metadata_df

    # ------------------------------------------------------------------
    # Filtering helpers
    # ------------------------------------------------------------------

    def get_strains(self) -> list[str]:
        cells = self.load_cell_index()
        return sorted(set(c["strain"] for c in cells))

    def get_contrasts_for_strain(self, strain: str) -> list[str]:
        cells = self.load_cell_index()
        return sorted(
            set(
                c["preferred_contrast"]
                for c in cells
                if c["strain"] == strain
            )
        )

    def get_cells_for_filter(
        self, strain: str, contrast: str
    ) -> list[dict]:
        cells = self.load_cell_index()
        return [
            c
            for c in cells
            if c["strain"] == strain
            and c["preferred_contrast"] == contrast
        ]

    def get_cell_by_id(self, cell_id: str) -> dict | None:
        for c in self.load_cell_index():
            if c["cell_id"] == cell_id:
                return c
        return None

    def get_cells_with_image(self, image_type: str) -> list[dict]:
        """Return cells that have a specific image type available."""
        key = f"has_{image_type}"
        return [
            c for c in self.load_cell_index() if c.get(key, False)
        ]

    # ------------------------------------------------------------------
    # Image paths
    # ------------------------------------------------------------------

    def get_image_url(self, cell: dict, image_type: str, thumb: bool = False) -> str | None:
        """Return the URL path for an image, or None if unavailable."""
        if not cell.get(f"has_{image_type}", False):
            return None
        if thumb:
            key = f"{image_type}_thumb"
        else:
            key = f"{image_type}_path"
        rel_path = cell.get(key)
        if rel_path:
            return f"/images/{rel_path}"
        return None

    # ------------------------------------------------------------------
    # Summary table for metadata tab
    # ------------------------------------------------------------------

    def get_metadata_summary(self) -> pd.DataFrame:
        cells = self.load_cell_index()
        if not cells:
            return pd.DataFrame()

        df = pd.DataFrame(cells)

        # Ensure expected columns exist (handles old/new cell_index.json)
        _image_cols = [
            "has_gridplot_1", "has_bar_polar", "has_polar_arrow",
            "has_flash_rf_4px", "has_flash_heatmap_4px",
            "has_flash_rf_6px", "has_flash_heatmap_6px",
            "has_gaussian_rf",
            "has_bar_flash_slow", "has_bar_flash_fast",
            # Legacy names (old cell_index.json compatibility)
            "has_flash_rf", "has_flash_rf_heatmap", "has_bar_flash",
        ]
        for col in _image_cols:
            if col not in df.columns:
                df[col] = False

        # Use new column names if available, fall back to legacy
        flash_col = "has_flash_rf_6px" if "has_flash_rf_6px" in df.columns and df["has_flash_rf_6px"].any() else "has_flash_rf"
        bar_flash_col = "has_bar_flash_slow" if "has_bar_flash_slow" in df.columns and df["has_bar_flash_slow"].any() else "has_bar_flash"

        summary = (
            df.groupby("strain")
            .agg(
                ON_cells=(
                    "preferred_contrast",
                    lambda x: (x == "ON").sum(),
                ),
                OFF_cells=(
                    "preferred_contrast",
                    lambda x: (x == "OFF").sum(),
                ),
                total_cells=("cell_id", "count"),
                with_gridplot=(
                    "has_gridplot_1",
                    lambda x: x.sum(),
                ),
                with_sweep=(
                    "has_bar_polar",
                    lambda x: x.sum(),
                ),
                with_flash_rf=(
                    flash_col,
                    lambda x: x.sum(),
                ),
                with_gaussian=(
                    "has_gaussian_rf",
                    lambda x: x.sum(),
                ),
                with_bar_flash=(
                    bar_flash_col,
                    lambda x: x.sum(),
                ),
                date_range=(
                    "date",
                    lambda x: f"{x.min()} to {x.max()}",
                ),
            )
            .reset_index()
        )
        summary.columns = [
            "Strain",
            "ON Cells",
            "OFF Cells",
            "Total",
            "With GridPlot",
            "With Sweep",
            "With Flash RF",
            "With Gaussian",
            "With Bar Flash",
            "Date Range",
        ]
        return summary
