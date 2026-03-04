from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class WombatState:
    tracers: dict[str, np.ndarray]


DEFAULT_TRACERS = {
    # WOMBAT-lite is treated as O2-NO3-Fe NPZD in this Python implementation.
    "wombat-lite": ("o2", "no3", "dfe", "phy", "zoo", "det"),
    "wombat-mid": ("o2", "no3", "po4", "dfe", "phy", "zoo", "det", "doc"),
}


class WombatModel:
    """Simplified Python implementation point for WOMBAT parameterization.

    NOTE: this implementation now includes coupled nutrient/oxygen NPZD tendencies
    for WOMBAT-lite so tracer concentrations evolve coherently.
    """

    def __init__(self, scheme: str, parameters: dict[str, float] | None = None):
        if scheme not in DEFAULT_TRACERS:
            raise ValueError(f"Unsupported scheme: {scheme}")
        self.scheme = scheme
        self.params = parameters or {}

    def empty_state(self, nz: int) -> WombatState:
        return WombatState(tracers={name: np.zeros(nz, dtype=float) for name in DEFAULT_TRACERS[self.scheme]})

    def _step_lite(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        tr = {k: v.copy() for k, v in state.tracers.items()}

        phy = tr["phy"]
        zoo = tr["zoo"]
        det = tr["det"]
        no3 = tr["no3"]
        dfe = tr["dfe"]
        o2 = tr["o2"]

        day = 86400.0
        # Biological rates (day^-1 defaults)
        mu_max = self.params.get("mu_max", self.params.get("phy_growth", 0.8)) / day
        phy_mort = self.params.get("phy_mortality", 0.05) / day
        zoo_mort = self.params.get("zoo_mortality", 0.03) / day
        g_max = self.params.get("grazing_max", 0.6) / day
        k_g = self.params.get("grazing_half_sat", 0.2)
        assim = self.params.get("grazing_assimilation", 0.7)
        remin = self.params.get("det_remin", 0.03) / day

        # Half saturation constants
        k_no3 = self.params.get("k_no3", 0.5)
        k_dfe = self.params.get("k_dfe", 1e-4)

        # Stoichiometry (per phy C-like unit)
        r_no3 = self.params.get("redfield_no3", 16.0 / 106.0)
        r_o2 = self.params.get("redfield_o2", 172.0 / 106.0)
        r_dfe = self.params.get("redfield_dfe", 1e-5)

        light_lim = np.clip(light, 0.0, 1.0)
        no3_lim = no3 / (no3 + k_no3 + 1e-12)
        dfe_lim = dfe / (dfe + k_dfe + 1e-12)
        nut_lim = np.minimum(no3_lim, dfe_lim)

        phy_growth = mu_max * light_lim * nut_lim * phy

        grazing = g_max * (phy * phy / (phy * phy + k_g * k_g + 1e-12)) * zoo
        phy_loss_mort = phy_mort * phy
        zoo_loss_mort = zoo_mort * zoo

        det_remin = remin * det

        dphy = phy_growth - grazing - phy_loss_mort
        dzoo = assim * grazing - zoo_loss_mort
        ddet = (1.0 - assim) * grazing + phy_loss_mort + zoo_loss_mort - det_remin

        dno3 = -r_no3 * phy_growth + r_no3 * det_remin
        ddfe = -r_dfe * phy_growth + r_dfe * det_remin
        do2 = r_o2 * phy_growth - r_o2 * det_remin

        phy_new = np.maximum(phy + dt_seconds * dphy, 0.0)
        zoo_new = np.maximum(zoo + dt_seconds * dzoo, 0.0)
        det_new = np.maximum(det + dt_seconds * ddet, 0.0)
        no3_new = np.maximum(no3 + dt_seconds * dno3, 0.0)
        dfe_new = np.maximum(dfe + dt_seconds * ddfe, 0.0)
        o2_new = np.maximum(o2 + dt_seconds * do2, 0.0)

        tr["phy"] = phy_new
        tr["zoo"] = zoo_new
        tr["det"] = det_new
        tr["no3"] = no3_new
        tr["dfe"] = dfe_new
        tr["o2"] = o2_new

        return WombatState(tracers=tr)

    def step(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        if self.scheme == "wombat-lite":
            return self._step_lite(state, dt_seconds=dt_seconds, light=light)

        # Keep wombat-mid placeholder behavior for now.
        tr = {k: v.copy() for k, v in state.tracers.items()}
        growth = self.params.get("phy_growth", 0.6) / 86400.0
        mortality = self.params.get("phy_mortality", 0.1) / 86400.0
        limiter = np.clip(light, 0.0, 1.0)
        dphy = dt_seconds * (growth * limiter * tr["phy"] - mortality * tr["phy"])
        tr["phy"] = np.maximum(tr["phy"] + dphy, 0.0)
        return WombatState(tracers=tr)
