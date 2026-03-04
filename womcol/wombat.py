from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class WombatState:
    tracers: dict[str, np.ndarray]


DEFAULT_TRACERS = {
    # WOMBAT-lite in this Python scaffold tracks core O2-NO3-Fe NPZD pools.
    "wombat-lite": ("o2", "no3", "dfe", "phy", "zoo", "det"),
    # WOMBAT-mid placeholder tracer set used by this column framework.
    "wombat-mid": ("o2", "no3", "po4", "dfe", "phy", "zoo", "det", "doc"),
}


class WombatModel:
    """Python WOMBAT process model aligned to upstream WOMBAT-lite/mid structure.

    This implementation remains reduced-order versus the full Fortran modules but
    now follows the same process blocks used in both upstream files:
    * temperature-scaled maximum growth,
    * multiplicative light and nutrient limitation,
    * quadratic detrital remineralisation,
    * adaptive zooplankton grazing efficiency,
    * linear + quadratic mortality terms,
    * explicit nutrient/oxygen/doc source-sink coupling.
    """

    def __init__(self, scheme: str, parameters: dict[str, float] | None = None):
        if scheme not in DEFAULT_TRACERS:
            raise ValueError(f"Unsupported scheme: {scheme}")
        self.scheme = scheme
        self.params = parameters or {}

    def empty_state(self, nz: int) -> WombatState:
        return WombatState(tracers={name: np.zeros(nz, dtype=float) for name in DEFAULT_TRACERS[self.scheme]})

    def _temperature_growth_scale(self) -> float:
        """Arrhenius-like growth scaling using WOMBAT-style `abioa * bbioa**T`.

        A scalar fallback is used because this lightweight interface does not pass
        a temperature profile into `step()`.
        """

        temp_c = self.params.get("temp_c", 15.0)
        abioa = self.params.get("abioa", self.params.get("mu_0", 0.6 / 86400.0))
        bbioa = self.params.get("bbioa", 1.066)
        return abioa * bbioa**temp_c

    @staticmethod
    def _zooplankton_grazing_rate(
        phy: np.ndarray,
        det: np.ndarray,
        zoo: np.ndarray,
        *,
        zprefphy: float,
        zprefdet: float,
        zoogmax: float,
        zooepsmin: float,
        zooepsmax: float,
        zooepsrat: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        total_pref = max(zprefphy + zprefdet, 1e-30)
        zooprefphy = zprefphy / total_pref
        zooprefdet = zprefdet / total_pref

        biophy = np.maximum(phy, 0.0)
        biodet = np.maximum(det, 0.0)
        zooprey = zooprefphy * biophy + zooprefdet * biodet

        g_peffect = np.exp(-zooprey * zooepsrat)
        zooeps = zooepsmin + (zooepsmax - zooepsmin) * g_peffect
        g_npz = zoogmax * (zooeps * zooprey * zooprey) / (zoogmax + zooeps * zooprey * zooprey + 1e-30)

        mask = zooprey > 1e-12
        graz_phy = np.where(mask, g_npz * zoo * (zooprefphy * biophy) / (zooprey + 1e-30), 0.0)
        graz_det = np.where(mask, g_npz * zoo * (zooprefdet * biodet) / (zooprey + 1e-30), 0.0)
        return graz_phy, graz_det, zooeps

    def _step_lite(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        tr = {k: v.copy() for k, v in state.tracers.items()}

        phy = tr["phy"]
        zoo = tr["zoo"]
        det = tr["det"]
        no3 = tr["no3"]
        dfe = tr["dfe"]
        o2 = tr["o2"]

        # Fortran-style coefficients (day^-1 converted to s^-1 where applicable).
        mumax = self.params.get("phy_mumax", self._temperature_growth_scale())
        phykn = self.params.get("phykn", self.params.get("k_no3", 0.5))
        phykf = self.params.get("phykf", self.params.get("k_dfe", 1e-4))
        phylmor = self.params.get("phylmor", self.params.get("phy_mortality", 0.05) / 86400.0)
        phyqmor = self.params.get("phyqmor", 0.03 / 86400.0)

        zoogmax = self.params.get("zoogmax", self.params.get("grazing_max", 0.6) / 86400.0)
        zooepsmin = self.params.get("zooepsmin", 0.35)
        zooepsmax = self.params.get("zooepsmax", 0.85)
        zooepsrat = self.params.get("zooepsrat", 0.5)
        zprefphy = self.params.get("zprefphy", 1.0)
        zprefdet = self.params.get("zprefdet", 0.35)
        zooCingest = self.params.get("zooCingest", 0.7)
        zooCassim = self.params.get("zooCassim", 0.7)
        zoolmor = self.params.get("zoolmor", self.params.get("zoo_mortality", 0.03) / 86400.0)
        zooqmor = self.params.get("zooqmor", 0.02 / 86400.0)
        zookz = self.params.get("zookz", 0.1)

        reminr = self.params.get("detlrem", self.params.get("det_remin", 0.03) / 86400.0)

        # C:N:Fe:O2 stoichiometry for reduced tracer set.
        r_no3 = self.params.get("redfield_no3", 16.0 / 106.0)
        r_dfe = self.params.get("redfield_dfe", 1.0e-5)
        r_o2 = self.params.get("redfield_o2", 172.0 / 106.0)

        # Support internal sub-stepping, similar to dt_npzd behaviour in Fortran.
        max_internal_dt = self.params.get("dt_npzd", dt_seconds)
        nsub = max(1, int(np.ceil(dt_seconds / max_internal_dt)))
        dts = dt_seconds / nsub

        light_lim = np.clip(light, 0.0, 1.0)

        for _ in range(nsub):
            phy_lnit = no3 / (no3 + phykn + 1e-30)
            phy_lfer = dfe / (dfe + phykf + 1e-30)
            phy_mu = mumax * light_lim * np.minimum(phy_lnit, phy_lfer)
            phygrow = phy_mu * phy

            graz_phy, graz_det, _ = self._zooplankton_grazing_rate(
                phy,
                det,
                zoo,
                zprefphy=zprefphy,
                zprefdet=zprefdet,
                zoogmax=zoogmax,
                zooepsmin=zooepsmin,
                zooepsmax=zooepsmax,
                zooepsrat=zooepsrat,
            )

            phymorl = phylmor * phy
            phymorq = phyqmor * phy * phy
            zoo_slmor = zoo / (zoo + zookz + 1e-30)
            zoomorl = zoolmor * zoo * zoo_slmor
            zoomorq = zooqmor * zoo * zoo
            detremi = reminr * det * det

            zoo_ingest = graz_phy + graz_det
            zoo_assim = zoo_ingest * zooCingest * zooCassim
            zoo_excr = zoo_ingest * zooCingest * (1.0 - zooCassim)
            zoo_eges = zoo_ingest * (1.0 - zooCingest)

            dphy = phygrow - graz_phy - phymorl - phymorq
            dzoo = zoo_assim - zoomorl - zoomorq
            ddet = phymorl + phymorq + zoomorl + zoomorq + zoo_eges - graz_det - detremi

            # Dissolved nutrients and oxygen.
            remin_source = detremi + zoo_excr
            dno3 = -r_no3 * phygrow + r_no3 * remin_source
            ddfe = -r_dfe * phygrow + r_dfe * remin_source
            do2 = r_o2 * phygrow - r_o2 * remin_source

            phy = np.maximum(phy + dts * dphy, 0.0)
            zoo = np.maximum(zoo + dts * dzoo, 0.0)
            det = np.maximum(det + dts * ddet, 0.0)
            no3 = np.maximum(no3 + dts * dno3, 0.0)
            dfe = np.maximum(dfe + dts * ddfe, 0.0)
            o2 = np.maximum(o2 + dts * do2, 0.0)

        tr["phy"] = phy
        tr["zoo"] = zoo
        tr["det"] = det
        tr["no3"] = no3
        tr["dfe"] = dfe
        tr["o2"] = o2

        return WombatState(tracers=tr)

    def _step_mid(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        tr = {k: v.copy() for k, v in state.tracers.items()}

        phy = tr["phy"]
        zoo = tr["zoo"]
        det = tr["det"]
        doc = tr["doc"]
        no3 = tr["no3"]
        po4 = tr["po4"]
        dfe = tr["dfe"]
        o2 = tr["o2"]

        mumax = self.params.get("phy_mumax", self._temperature_growth_scale())
        k_no3 = self.params.get("phykn", 0.5)
        k_po4 = self.params.get("phykp", 0.03)
        k_dfe = self.params.get("phykf", 1e-4)
        overflow = self.params.get("overflow", 0.08)

        phylmor = self.params.get("phylmor", 0.04 / 86400.0)
        phyqmor = self.params.get("phyqmor", 0.02 / 86400.0)
        reminr = self.params.get("detlrem", 0.03 / 86400.0)
        doc_remin = self.params.get("doc_remin", 0.015 / 86400.0)

        zoogmax = self.params.get("zoogmax", 0.7 / 86400.0)
        zooepsmin = self.params.get("zooepsmin", 0.35)
        zooepsmax = self.params.get("zooepsmax", 0.85)
        zooepsrat = self.params.get("zooepsrat", 0.5)
        zprefphy = self.params.get("zprefphy", 1.0)
        zprefdet = self.params.get("zprefdet", 0.35)
        zooCingest = self.params.get("zooCingest", 0.7)
        zooCassim = self.params.get("zooCassim", 0.7)
        zoolmor = self.params.get("zoolmor", 0.03 / 86400.0)
        zooqmor = self.params.get("zooqmor", 0.02 / 86400.0)
        zookz = self.params.get("zookz", 0.1)

        r_no3 = self.params.get("redfield_no3", 16.0 / 106.0)
        r_po4 = self.params.get("redfield_po4", 1.0 / 106.0)
        r_dfe = self.params.get("redfield_dfe", 1.0e-5)
        r_o2 = self.params.get("redfield_o2", 172.0 / 106.0)

        light_lim = np.clip(light, 0.0, 1.0)

        l_no3 = no3 / (no3 + k_no3 + 1e-30)
        l_po4 = po4 / (po4 + k_po4 + 1e-30)
        l_dfe = dfe / (dfe + k_dfe + 1e-30)
        nut_lim = np.minimum(np.minimum(l_no3, l_po4), l_dfe)

        phy_mu = mumax * light_lim * nut_lim
        gross_growth = phy_mu * phy
        phydoc = overflow * gross_growth
        phygrow = np.maximum(gross_growth - phydoc, 0.0)

        graz_phy, graz_det, _ = self._zooplankton_grazing_rate(
            phy,
            det,
            zoo,
            zprefphy=zprefphy,
            zprefdet=zprefdet,
            zoogmax=zoogmax,
            zooepsmin=zooepsmin,
            zooepsmax=zooepsmax,
            zooepsrat=zooepsrat,
        )

        phymorl = phylmor * phy
        phymorq = phyqmor * phy * phy
        zoo_slmor = zoo / (zoo + zookz + 1e-30)
        zoomorl = zoolmor * zoo * zoo_slmor
        zoomorq = zooqmor * zoo * zoo

        detremi = reminr * det * det
        docremi = doc_remin * doc

        zoo_ingest = graz_phy + graz_det
        zoo_assim = zoo_ingest * zooCingest * zooCassim
        zoo_excr = zoo_ingest * zooCingest * (1.0 - zooCassim)
        zoo_eges = zoo_ingest * (1.0 - zooCingest)

        dphy = phygrow - graz_phy - phymorl - phymorq
        dzoo = zoo_assim - zoomorl - zoomorq
        ddet = phymorl + phymorq + zoomorl + zoomorq + zoo_eges - graz_det - detremi
        ddoc = phydoc + zoo_excr + detremi - docremi

        remin_source = docremi
        dno3 = -r_no3 * gross_growth + r_no3 * remin_source
        dpo4 = -r_po4 * gross_growth + r_po4 * remin_source
        ddfe = -r_dfe * gross_growth + r_dfe * remin_source
        do2 = r_o2 * gross_growth - r_o2 * remin_source

        tr["phy"] = np.maximum(phy + dt_seconds * dphy, 0.0)
        tr["zoo"] = np.maximum(zoo + dt_seconds * dzoo, 0.0)
        tr["det"] = np.maximum(det + dt_seconds * ddet, 0.0)
        tr["doc"] = np.maximum(doc + dt_seconds * ddoc, 0.0)
        tr["no3"] = np.maximum(no3 + dt_seconds * dno3, 0.0)
        tr["po4"] = np.maximum(po4 + dt_seconds * dpo4, 0.0)
        tr["dfe"] = np.maximum(dfe + dt_seconds * ddfe, 0.0)
        tr["o2"] = np.maximum(o2 + dt_seconds * do2, 0.0)

        return WombatState(tracers=tr)

    def step(self, state: WombatState, *, dt_seconds: float, light: np.ndarray) -> WombatState:
        if self.scheme == "wombat-lite":
            return self._step_lite(state, dt_seconds=dt_seconds, light=light)
        return self._step_mid(state, dt_seconds=dt_seconds, light=light)
