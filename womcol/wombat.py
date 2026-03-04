from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class WombatState:
    tracers: dict[str, np.ndarray]


DEFAULT_TRACERS = {
    "wombat-lite": ("o2", "no3", "phy", "zoo", "det"),
    "wombat-mid": ("o2", "no3", "po4", "dfe", "phy", "zoo", "det", "doc"),
}


class WombatModel:
    """Simplified Python implementation point for WOMBAT parameterization.

    The biological source/sink terms should be kept synchronized with the canonical
    Fortran in ACCESS-NRI/GFDL-generic-tracers. See scripts/sync_wombat_sources.py.
    """

    def __init__(self, scheme: str, parameters: dict[str, float] | None = None):
        if scheme not in DEFAULT_TRACERS:
            raise ValueError(f"Unsupported scheme: {scheme}")
        self.scheme = scheme
        self.params = parameters or {}

    def empty_state(self, nz: int) -> WombatState:
        return WombatState(tracers={name: np.zeros(nz, dtype=float) for name in DEFAULT_TRACERS[self.scheme]})

    def step(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        new = {k: v.copy() for k, v in state.tracers.items()}
        growth = self.params.get("phy_growth", 0.6) / 86400.0
        mortality = self.params.get("phy_mortality", 0.1) / 86400.0

        if "phy" in new:
            limiter = np.clip(light, 0.0, 1.0)
            dphy = dt_seconds * (growth * limiter * new["phy"] - mortality * new["phy"])
            new["phy"] = np.maximum(new["phy"] + dphy, 0.0)
        return WombatState(tracers=new)
