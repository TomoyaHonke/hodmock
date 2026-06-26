"""
halomodel.py - Halo model calculator for galaxy correlation functions.

Designed for repeated HOD evaluations in MCMC:
  - HOD-independent quantities (HMF, Pmm, NFW profiles, bias) are precomputed once.
  - Only HOD-dependent integrals are recomputed each call.
  - NFW Fourier transform vectorized over k for each M (eliminates inner loop).

Usage
-----
    hm = HaloModel(k_arr, M_arr, z)
    xi = hm.compute_xi(hod, r_arr)
    xi, Nc, Ns = hm.compute_xi(hod, r_arr, return_Nc_Ns=True)
    ng = hm.compute_ng(hod)

Recompute after changing halo model parameters:
    hm.update_halo_model(hmf_model='press74')
"""

import numpy as np
import camb
from colossus.cosmology import cosmology as col_cosmo
from colossus.lss import bias as col_bias, mass_function
from colossus.halo import concentration, mass_so
from scipy.special import erf
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.interpolate import interp1d
from mcfit import P2xi


# Default: Planck 2018
PLANCK18 = {
    'H0': 67.4,
    'ombh2': 0.0224,
    'omch2': 0.12,
    'As': 2.1e-9,
    'ns': 0.965,
}


# ---------------------------------------------------------------------------
# HOD occupation functions
# ---------------------------------------------------------------------------

def N_cen(M, hod):
    """
    Central galaxy occupation <Nc>(M).

    hod['model'] selects the functional form:

    'more'  : standard step-function HOD (More+ 2015)
              requires: logMmin, sigma_logM
              Nc = 0.5 * (1 + erf((logM - logMmin) / sigma_logM))

    'ghod'  : Gaussian
              requires: Ac, logMc, sigma_M

    'lnhod' : log-normal in logM - (logMc - 1)
              requires: Ac, logMc, sigma_M

    'sfhod' : Gaussian below logMc, power law above
              requires: Ac, logMc, sigma_M, gamma

    'mhmq'  : modified Gaussian with skewness via erf
              requires: Ac, logMc, sigma_M, gamma
    """
    logM = np.log10(M)
    model = hod['model']

    if model == 'more':
        return 0.5 * (1.0 + erf((logM - hod['logMmin']) / hod['sigma_logM']))

    Ac = hod['Ac']
    sigma_M = hod['sigma_M']
    logMc = hod['logMc']

    if model == 'ghod':
        return Ac / (np.sqrt(2 * np.pi) * sigma_M) * np.exp(-0.5 * ((logM - logMc) / sigma_M)**2)

    elif model == 'lnhod':
        x = logM - (logMc - 1.0)
        Nc = np.zeros_like(M, dtype=float)
        mask = x > 0
        Nc[mask] = (
            Ac / (np.sqrt(2 * np.pi) * sigma_M * x[mask])
            * np.exp(-0.5 * (np.log(x[mask]) / sigma_M)**2)
        )
        return Nc

    elif model == 'sfhod':
        Nc = np.zeros_like(M, dtype=float)
        lo = logM < logMc
        Nc[lo] = Ac / (np.sqrt(2 * np.pi) * sigma_M) * np.exp(-0.5 * ((logM[lo] - logMc) / sigma_M)**2)
        Nc[~lo] = Ac / (np.sqrt(2 * np.pi) * sigma_M) * (M[~lo] / 10**logMc)**hod['gamma']
        return Nc

    elif model == 'mhmq':
        Nc_g = Ac / (np.sqrt(2 * np.pi) * sigma_M) * np.exp(-0.5 * ((logM - logMc) / sigma_M)**2)
        return Nc_g * (1.0 + erf(hod['gamma'] * (logM - logMc) / (np.sqrt(2) * sigma_M)))

    else:
        raise ValueError(f"Unknown Nc model: {model!r}")


def N_sat(M, hod, Nc=None):
    """
    Satellite galaxy occupation <Ns>(M).

    hod.get('sat_model', 'elg') selects the prefactor:

    'elg'  : Ns = As * ((M - M0) / M1)^alpha   [ELG-style, As is free]
    'more' : Ns = Nc * ((M - M0) / M1)^alpha    [LRG-style, coupled to Nc]

    If hod.get('conformity', False) is True or 'strict' and sat_model=='elg':
        Ns = Nc * Ns_base  (strict central-satellite conformity)
    """
    M0 = 10**hod['logM0']
    M1 = 10**hod['logM1']
    alpha = hod['alpha']
    sat_model = hod.get('sat_model', 'elg')

    Ns = np.zeros_like(M, dtype=float)
    mask = M > M0

    if sat_model == 'more':
        if Nc is None:
            Nc = N_cen(M, hod)
        Ns[mask] = Nc[mask] * ((M[mask] - M0) / M1)**alpha
    else:
        As = hod['As']
        Ns[mask] = As * ((M[mask] - M0) / M1)**alpha
        if hod.get('conformity', False) is True or hod.get('conformity', False) == 'strict':
            if Nc is None:
                Nc = N_cen(M, hod)
            Ns = Nc * Ns

    return Ns


def galaxy_number_density(M, dndM, Nc, Ns):
    """Integrate dndM * (Nc + Ns) over M."""
    return float(simpson(dndM * (Nc + Ns), x=M))


# ---------------------------------------------------------------------------
# Internal precomputation helpers
# ---------------------------------------------------------------------------

def _compute_Pmm(k_arr, z, cosmo_params):
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=cosmo_params['H0'],
        ombh2=cosmo_params['ombh2'],
        omch2=cosmo_params['omch2'],
    )
    pars.InitPower.set_params(As=cosmo_params['As'], ns=cosmo_params['ns'])
    pars.set_matter_power(redshifts=[z], kmax=float(np.max(k_arr)))
    results = camb.get_results(pars)
    k_camb, _, pk = results.get_matter_power_spectrum(
        minkh=float(np.min(k_arr)),
        maxkh=float(np.max(k_arr)),
        npoints=500,
    )
    Pmm_interp = interp1d(k_camb, pk[0], kind='linear', bounds_error=False, fill_value='extrapolate')
    return Pmm_interp(k_arr)


def _compute_us(k_arr, M_arr, z, rs_arr, Nr=512):
    """
    NFW profile Fourier transform u(k, M), shape (Nk, NM).

    Vectorized over k for each M: reduces the Python loop from
    O(Nk * NM) to O(NM), with numpy handling the k dimension.
    """
    Nk = len(k_arr)
    NM = len(M_arr)
    us = np.zeros((Nk, NM))

    for iM in range(NM):
        R200_Mpc = mass_so.M_to_R(M_arr[iM], z, '200m') / 1000.0  # kpc/h -> Mpc/h
        rs_Mpc = rs_arr[iM] / 1000.0

        r = np.linspace(1e-6 * rs_Mpc, R200_Mpc, Nr)
        x = r / rs_Mpc
        rho = 1.0 / (x * (1.0 + x)**2)
        weight = 4.0 * np.pi * r**2 * rho  # (Nr,)

        # kr: (Nk, Nr) — vectorize k integration
        kr = k_arr[:, None] * r[None, :]
        sinkr = np.where(np.abs(kr) < 1e-10, 1.0, np.sin(kr) / kr)

        integral = simpson(weight[None, :] * sinkr, x=r, axis=1)  # (Nk,)
        norm = simpson(weight, x=r)
        us[:, iM] = integral / norm

    return us


# ---------------------------------------------------------------------------
# HaloModel class
# ---------------------------------------------------------------------------

class HaloModel:
    """
    Halo model calculator for galaxy correlation functions.

    Precomputes all HOD-independent quantities on construction.
    Call compute_xi() / compute_ng() repeatedly (e.g. in MCMC) at low cost.

    Parameters
    ----------
    k_arr : array_like
        Wavenumber array [h/Mpc]. Log-spaced, e.g. np.logspace(-3, 2, 200).
    M_arr : array_like
        Halo mass array [Msun/h]. Log-spaced, e.g. np.logspace(10, 15, 200).
    z : float
        Redshift.
    cosmo_params : dict, optional
        CAMB cosmological parameters. Keys: H0, ombh2, omch2, As, ns.
        Defaults to Planck 2018.
    colossus_cosmo : str, optional
        Colossus cosmology name (must be consistent with cosmo_params).
        Default: 'planck18'.
    hmf_model : str
        Colossus HMF model. Default: 'tinker08'.
    conc_model : str
        Colossus concentration model. Default: 'diemer19'.
    bias_model : str
        Colossus halo bias model. Default: 'tinker10'.
    Nr : int
        Radial grid points for NFW profile integration. Default: 512.

    Notes
    -----
    update_halo_model(**kwargs) resets any of the above parameters and
    triggers a full recompute. Use this when switching HMF / cosmology.
    """

    def __init__(
        self,
        k_arr,
        M_arr,
        z,
        cosmo_params=None,
        colossus_cosmo='planck18',
        hmf_model='tinker08',
        conc_model='diemer19',
        bias_model='tinker10',
        Nr=512,
    ):
        self.k_arr = np.asarray(k_arr, dtype=float)
        self.M_arr = np.asarray(M_arr, dtype=float)
        self.z = float(z)
        self.cosmo_params = PLANCK18.copy() if cosmo_params is None else dict(cosmo_params)
        self.colossus_cosmo = colossus_cosmo
        self.hmf_model = hmf_model
        self.conc_model = conc_model
        self.bias_model = bias_model
        self.Nr = Nr

        col_cosmo.setCosmology(self.colossus_cosmo)
        self._precompute()

    # ------------------------------------------------------------------
    # Precomputation
    # ------------------------------------------------------------------

    def _precompute(self):
        """Compute all HOD-independent quantities and cache them."""
        M, k, z = self.M_arr, self.k_arr, self.z

        print("  [HaloModel] computing HMF ...", end=' ', flush=True)
        M2dndM = mass_function.massFunction(M, z, mdef='200m', model=self.hmf_model, q_out='dndlnM')
        self.dndM = M2dndM / M
        print("done")

        print("  [HaloModel] computing Pmm ...", end=' ', flush=True)
        self.Pmm = _compute_Pmm(k, z, self.cosmo_params)
        print("done")

        print("  [HaloModel] computing concentration / scale radius ...", end=' ', flush=True)
        self.c_arr = concentration.concentration(M, '200m', z, model=self.conc_model)
        R200 = mass_so.M_to_R(M, z, '200m')   # kpc/h
        self.rs_arr = R200 / self.c_arr         # kpc/h
        print("done")

        print(f"  [HaloModel] computing NFW u(k,M) [{len(k)}k x {len(M)}M] ...", end=' ', flush=True)
        self.us = _compute_us(k, M, z, self.rs_arr, Nr=self.Nr)
        print("done")

        print("  [HaloModel] computing halo bias ...", end=' ', flush=True)
        self.bias_arr = col_bias.haloBias(M, z, mdef='200m', model=self.bias_model)
        print("done")

        # P(k) -> xi(r) transformer
        self.xi_func = P2xi(k, l=0)

        # Products that are constant w.r.t. HOD
        self._dndM_bias = self.dndM * self.bias_arr  # (NM,)
        self._rs_Mpc = self.rs_arr / 1000.0           # Mpc/h, for off-centering

        # Growth rate f = Omega_m(z)^0.55
        cosmo = col_cosmo.getCurrent()
        H0 = cosmo.Hz(0.0)
        Hz = cosmo.Hz(z)
        Ez2 = (Hz / H0) ** 2
        Om_z = cosmo.Om0 * (1.0 + z) ** 3 / Ez2
        self.f_growth = Om_z ** 0.55

    def update_halo_model(self, **kwargs):
        """
        Update halo model parameters and recompute precomputed quantities.

        Accepted keyword arguments
        --------------------------
        hmf_model, conc_model, bias_model, Nr : str / int
            Halo structure parameters.
        cosmo_params : dict
            CAMB cosmological parameters.
        colossus_cosmo : str
            Colossus cosmology name.
        z : float
            Redshift.
        k_arr, M_arr : array_like
            Grid arrays.
        """
        valid = {'hmf_model', 'conc_model', 'bias_model', 'Nr',
                 'cosmo_params', 'colossus_cosmo', 'z', 'k_arr', 'M_arr'}
        for key, val in kwargs.items():
            if key not in valid:
                raise ValueError(f"Unknown parameter: {key!r}")
            if key in ('k_arr', 'M_arr'):
                setattr(self, key, np.asarray(val, dtype=float))
            else:
                setattr(self, key, val)

        if 'colossus_cosmo' in kwargs:
            col_cosmo.setCosmology(self.colossus_cosmo)

        self._precompute()

    # ------------------------------------------------------------------
    # Internal HOD / power spectrum computation
    # ------------------------------------------------------------------

    def _compute_hod(self, hod):
        """Return Nc(M), Ns(M) arrays."""
        Nc = N_cen(self.M_arr, hod)
        Ns = N_sat(self.M_arr, hod, Nc=Nc)
        return Nc, Ns

    def _make_Pgg(self, Nc, Ns, hod, rsd=False, rsd_mode='2h'):
        """
        Compute P_gg(k) from precomputed arrays and HOD.

        Parameters
        ----------
        rsd : bool
            If True, apply Kaiser RSD monopole correction.
        rsd_mode : {'2h', 'full'}
            '2h'   — apply Kaiser factor only to the 2-halo term (physically motivated).
            'full' — apply Kaiser factor to the entire P_gg (simple approximation).
        """
        M, k = self.M_arr, self.k_arr

        strict_conformity = (
            (hod.get('conformity', False) is True or hod.get('conformity', False) == 'strict')
            and hod.get('sat_model', 'elg') == 'elg'
        )

        ng = galaxy_number_density(M, self.dndM, Nc, Ns)

        # Central off-centering factor: (Nk, NM) if poff > 0, else 1.0
        poff = hod.get('poff', 0.0)
        if poff > 0.0:
            Roff = hod['Roff']
            f_off = (
                1.0 - poff
                + poff * np.exp(-0.5 * k[:, None]**2 * (self._rs_Mpc[None, :] * Roff)**2)
            )
        else:
            f_off = 1.0  # scalar, broadcasts for free

        dndM = self.dndM  # (NM,)

        # Window functions (Nk, NM)
        Hc = Nc[None, :] * f_off / ng
        Hs = Ns[None, :] * self.us / ng

        # 1-halo term
        if strict_conformity:
            M0 = 10**hod['logM0']
            M1 = 10**hod['logM1']
            Ns_base = np.zeros_like(M, dtype=float)
            mask = M > M0
            Ns_base[mask] = hod['As'] * ((M[mask] - M0) / M1)**hod['alpha']

            P1h_cs = simpson(
                Nc[None, :] * Ns_base[None, :] * f_off * self.us * dndM[None, :] / ng**2,
                x=M,
                axis=1,
            )
            P1h_ss = simpson(
                Nc[None, :] * Ns_base[None, :]**2 * self.us**2 * dndM[None, :] / ng**2,
                x=M,
                axis=1,
            )
        else:
            P1h_cs = simpson(Hc * Hs * dndM[None, :], x=M, axis=1)
            P1h_ss = simpson(Hs * Hs * dndM[None, :], x=M, axis=1)
        P1h = 2.0 * P1h_cs + P1h_ss

        # 2-halo term (use precomputed dndM * bias)
        Ic = simpson(Hc * self._dndM_bias[None, :], x=M, axis=1)
        Is = simpson(Hs * self._dndM_bias[None, :], x=M, axis=1)
        P2h = self.Pmm * (Ic + Is)**2

        if not rsd:
            return P1h + P2h

        # Finger-of-God damping on the 1-halo term (isotropic, k-only approximation
        # of satellite/virial velocity dispersion smearing pairs along the line of
        # sight in redshift space). Real-space pairs (rsd=False, handled above) are
        # unaffected. Off by default: hod.get('sigma_fog', 0.0) == 0 reproduces the
        # previous behaviour exactly (damping factor = 1).
        sigma_fog = hod.get('sigma_fog', 0.0)
        if sigma_fog > 0.0:
            P1h = P1h * np.exp(-0.5 * (k * sigma_fog) ** 2)

        # Kaiser monopole factor: (1 + 2/3 beta + 1/5 beta^2), beta = f / b_eff
        b_eff = float(simpson((Nc + Ns) * self.dndM * self.bias_arr, x=M) / ng)
        beta = self.f_growth / b_eff
        kaiser = 1.0 + (2.0 / 3.0) * beta + (1.0 / 5.0) * beta ** 2

        if rsd_mode == '2h':
            return P1h + P2h * kaiser
        elif rsd_mode == 'full':
            return (P1h + P2h) * kaiser
        else:
            raise ValueError(f"Unknown rsd_mode: {rsd_mode!r}. Choose '2h' or 'full'.")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def compute_xi(self, hod, r_arr, return_Nc_Ns=False, rsd=False, rsd_mode='2h'):
        """
        Compute the galaxy two-point correlation function xi(r).

        Parameters
        ----------
        hod : dict
            HOD parameters. Must include 'model' key.
            See N_cen / N_sat docstrings for required keys per model.
        r_arr : array_like
            Separation array [Mpc/h].
        return_Nc_Ns : bool
            If True, also return Nc(M) and Ns(M).
        rsd : bool
            If True, apply Kaiser RSD monopole correction (large-scale approximation).
            Small-scale (1-halo) RSD is not modelled; restrict to r > ~10 Mpc/h.
        rsd_mode : {'2h', 'full'}
            '2h'   — Kaiser factor applied only to the 2-halo term (default).
            'full' — Kaiser factor applied to the entire P_gg.

        Returns
        -------
        xi : ndarray  shape (len(r_arr),)
        Nc  : ndarray  shape (len(M_arr),)   — only if return_Nc_Ns=True
        Ns  : ndarray  shape (len(M_arr),)   — only if return_Nc_Ns=True
        """
        r_arr = np.asarray(r_arr, dtype=float)
        Nc, Ns = self._compute_hod(hod)
        Pgg = self._make_Pgg(Nc, Ns, hod, rsd=rsd, rsd_mode=rsd_mode)

        r_fft, xi_fft = self.xi_func(Pgg)
        xi = interp1d(r_fft, xi_fft, bounds_error=False, fill_value='extrapolate')(r_arr)

        if return_Nc_Ns:
            return xi, Nc, Ns
        return xi

    def compute_ng(self, hod, return_Nc_Ns=False):
        """
        Compute galaxy number density ng [h^3/Mpc^3].

        Parameters
        ----------
        hod : dict
        return_Nc_Ns : bool

        Returns
        -------
        ng : float
        Nc, Ns : ndarray  — only if return_Nc_Ns=True
        """
        Nc, Ns = self._compute_hod(hod)
        ng = galaxy_number_density(self.M_arr, self.dndM, Nc, Ns)
        if return_Nc_Ns:
            return ng, Nc, Ns
        return ng

    def compute_Pgg(self, hod, return_Nc_Ns=False, rsd=False, rsd_mode='2h'):
        """
        Compute galaxy power spectrum P_gg(k) [Mpc^3/h^3].

        Parameters
        ----------
        hod : dict
        return_Nc_Ns : bool
        rsd : bool
            If True, apply Kaiser RSD monopole correction.
        rsd_mode : {'2h', 'full'}
            '2h'   — Kaiser factor applied only to the 2-halo term (default).
            'full' — Kaiser factor applied to the entire P_gg.

        Returns
        -------
        Pgg : ndarray  shape (len(k_arr),)
        Nc, Ns : ndarray  — only if return_Nc_Ns=True
        """
        Nc, Ns = self._compute_hod(hod)
        Pgg = self._make_Pgg(Nc, Ns, hod, rsd=rsd, rsd_mode=rsd_mode)
        if return_Nc_Ns:
            return Pgg, Nc, Ns
        return Pgg

    def compute_xi_multipoles(self, hod, r_arr):
        """
        Compute the real-space monopole xi0(r) and the linear-Kaiser quadrupole
        xi2(r) predicted by the model, following Eq. (53) of Cacciato et al. 2013
        (arXiv:1206.6890):

            xi2(r) = (4/3 beta + 4/7 beta^2) * xi0(r) - 3 * J3(r)
            J3(r)  = (1/r^3) * integral_0^r  xi0(y) y^2 dy
            beta   = f_growth / b_eff

        This is a *model-predicted* quadrupole (no observational quadrupole data
        is needed): it only uses quantities already computed for the monopole
        (real-space Pgg, halo bias, growth rate). It is the model-side ingredient
        needed to build the anisotropic xi(rp, pi) used by compute_wp(rsd_quad=True).

        Parameters
        ----------
        hod : dict
        r_arr : array_like
            Separation array [Mpc/h] at which to evaluate xi0, xi2.

        Returns
        -------
        xi0 : ndarray  shape (len(r_arr),)
        xi2 : ndarray  shape (len(r_arr),)
        """
        r_arr = np.asarray(r_arr, dtype=float)
        M = self.M_arr

        Nc, Ns = self._compute_hod(hod)
        Pgg_real = self._make_Pgg(Nc, Ns, hod, rsd=False)

        r_fft, xi0_fft = self.xi_func(Pgg_real)

        ng = galaxy_number_density(M, self.dndM, Nc, Ns)
        b_eff = float(simpson((Nc + Ns) * self.dndM * self.bias_arr, x=M) / ng)
        beta = self.f_growth / b_eff

        integrand = xi0_fft * r_fft ** 2
        cum = cumulative_trapezoid(integrand, x=r_fft, initial=0.0)
        J3 = cum / r_fft ** 3
        xi2_fft = (4.0 / 3.0 * beta + 4.0 / 7.0 * beta ** 2) * xi0_fft - 3.0 * J3

        xi0 = interp1d(r_fft, xi0_fft, bounds_error=False, fill_value='extrapolate')(r_arr)
        xi2 = interp1d(r_fft, xi2_fft, bounds_error=False, fill_value='extrapolate')(r_arr)
        return xi0, xi2

    def compute_wp(self, hod, rp_arr, pi_arr, rsd_quad=False):
        """
        Compute the projected correlation function wp(rp) by integrating
        xi(rp, pi) over the (finite) line-of-sight range given by pi_arr
        (Davis & Peebles 1983 Abel transform; Eq. 45-46 of Cacciato et al. 2013,
        arXiv:1206.6890).

        Parameters
        ----------
        hod : dict
        rp_arr : array_like
            Projected separation array [Mpc/h].
        pi_arr : array_like
            Line-of-sight separation array [Mpc/h] (finite range, matching the
            range used for the observational wp measurement).
        rsd_quad : bool
            False (default) — legacy behaviour, kept for backward compatibility:
            project the *real-space* monopole xi0(r) only, i.e. treat xi(rp, pi)
            as isotropic. This is exact only in the pi_max -> infinity limit
            (Eq. 46) and is what compute_xi(hod, r, rsd=False) + a manual
            line-of-sight integral has always computed in run_mcmc.py /
            analyze notebooks.
            True — include the model-predicted Kaiser quadrupole correction
            (Eq. 53), xi(rp, pi) = xi0(r) + xi2(r) * P2(mu), before integrating
            over the finite pi_arr range. This is a closer approximation to the
            finite-rmax redshift-space wp than the rsd_quad=False default,
            without requiring any observational quadrupole data.

        Returns
        -------
        wp : ndarray  shape (len(rp_arr),)
        """
        rp_arr = np.asarray(rp_arr, dtype=float)
        pi_arr = np.asarray(pi_arr, dtype=float)
        r_2d = np.sqrt(rp_arr[:, None] ** 2 + pi_arr[None, :] ** 2)

        if not rsd_quad:
            xi0_flat = self.compute_xi(hod, r_2d.ravel(), rsd=False)
            xi_2d = xi0_flat.reshape(r_2d.shape)
        else:
            r_flat = r_2d.ravel()
            xi0_flat, xi2_flat = self.compute_xi_multipoles(hod, r_flat)
            pi_2d = np.broadcast_to(pi_arr[None, :], r_2d.shape)
            mu_flat = np.divide(
                pi_2d.ravel(), r_flat,
                out=np.zeros_like(r_flat), where=(r_flat > 0),
            )
            P2 = 0.5 * (3.0 * mu_flat ** 2 - 1.0)
            xi_2d = (xi0_flat + xi2_flat * P2).reshape(r_2d.shape)

        wp = np.trapz(xi_2d, x=pi_arr, axis=1)
        return wp

    def __repr__(self):
        return (
            f"HaloModel(z={self.z}, Nk={len(self.k_arr)}, NM={len(self.M_arr)}, "
            f"hmf={self.hmf_model!r}, conc={self.conc_model!r}, bias={self.bias_model!r})"
        )
