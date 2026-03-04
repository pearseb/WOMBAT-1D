from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field, model_validator


class ColumnConfig(BaseModel):
    depth_m: float = Field(200.0, gt=0)
    nz: int = Field(100, ge=10)
    dt_seconds: float = Field(1800.0, gt=0)


class DataConfig(BaseModel):
    jra55_uri: str | None = None
    init_uri: str | None = None
    cache_dir: Path = Path("./data_cache")


class ModelConfig(BaseModel):
    scheme: Literal["wombat-lite", "wombat-mid"] = "wombat-lite"
    use_gotm: bool = True
    gotm_executable: str = "gotm"
    gotm_setup_dir: Path | None = None
    gotm_run_dir: Path = Path("./runs/gotm")
    gotm_output_path: str = "output.nc"
    parameters: dict[str, float] = Field(default_factory=dict)


class RuntimeConfig(BaseModel):
    year: int = Field(..., ge=1958)
    latitude: float = Field(..., ge=-90, le=90)
    longitude: float = Field(..., ge=-180, le=180)
    days: int = Field(365, gt=0, le=366)
    column: ColumnConfig = Field(default_factory=ColumnConfig)
    data: DataConfig = Field(default_factory=DataConfig)
    model: ModelConfig = Field(default_factory=ModelConfig)

    @model_validator(mode="after")
    def _normalize_lon(self) -> "RuntimeConfig":
        if self.longitude > 180:
            self.longitude -= 360
        return self
