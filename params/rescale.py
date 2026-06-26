"""
rescale.py — Ac/As の rescaling によって PFS ターゲット ng に合わせる。

使い方:
    from hodmock.params.rescale import compute_all_params, save_params
    from hodmock.config import HodMockConfig

    cfg = HodMockConfig()                  # デフォルト（Uchuu + PFS）
    params = compute_all_params(cfg)       # 全 z でベースライン ng 計算 + rescaling
    save_params(params, cfg)               # → hod_params_rescale.yaml

    # z_list だけ差し替えて別カタログにも使える
    cfg2 = HodMockConfig(z_list=[0.5, 1.0, 1.5])
    params2 = compute_all_params(cfg2)
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import numpy as np
import yaml

from hodmock.config import DEFAULT_CONFIG, HodMockConfig

# hodmock.halomodel（CAMB が必要）は compute_baseline_ng 内で遅延インポートする。
# これにより、rescale_params / save_params など CAMB 不要の関数は
# CAMB がない環境でも使える。


# ─────────────────────────────────────────────────────────────────────────────
# 内部ユーティリティ
# ─────────────────────────────────────────────────────────────────────────────

def _find_bin_index(
    z: float,
    bins: list[tuple[float, float]],
    clamp: bool = False,
) -> int:
    """z がどの bin に入るかを返す。

    Parameters
    ----------
    z : float
    bins : list of (zmin, zmax)
        最後の bin は zmax を含む（閉区間）。それ以外は zmax を含まない（半開区間）。
    clamp : bool
        True のとき、範囲外は端の bin のインデックスを返す（外挿）。
    """
    for i, (zmin, zmax) in enumerate(bins):
        if i == len(bins) - 1:
            if zmin <= z <= zmax:
                return i
        else:
            if zmin <= z < zmax:
                return i

    if clamp:
        return 0 if z < bins[0][0] else len(bins) - 1

    raise ValueError(f"z={z} は bins の範囲外: {bins}")


# ─────────────────────────────────────────────────────────────────────────────
# パブリック API
# ─────────────────────────────────────────────────────────────────────────────

def compute_baseline_ng(
    cfg: HodMockConfig = DEFAULT_CONFIG,
    verbose: bool = True,
) -> dict[float, float]:
    """cfg.z_list の各 z でベースライン HOD の ng を計算する。

    halomodel_module.HaloModel を使って積分し、その ng 値を返す。
    rescale_params() の分母（ng_hod）として使う。

    Parameters
    ----------
    cfg : HodMockConfig
    verbose : bool
        True のとき進捗を表示する。

    Returns
    -------
    dict[float, float]
        {z: ng [(Mpc/h)^-3]}
    """
    from hodmock.halomodel import HaloModel  # noqa: PLC0415  (CAMB が必要)

    ng_dict: dict[float, float] = {}

    for z in cfg.z_list:
        hod_i = _find_bin_index(z, cfg.hod_bins, clamp=True)
        hod = cfg.hod_baseline[hod_i]

        with contextlib.redirect_stdout(io.StringIO()):
            hm = HaloModel(
                cfg.k_arr, cfg.m_arr, z,
                cosmo_params=cfg.cosmo_params,
                colossus_cosmo=cfg.colossus_cosmo,
                hmf_model=cfg.hmf_model,
                conc_model=cfg.conc_model,
                bias_model=cfg.bias_model,
                Nr=cfg.halomodel_Nr,
            )

        ng = hm.compute_ng(hod)
        ng_dict[z] = ng

        if verbose:
            print(f"  z={z:.2f}  ng={ng:.4e}", flush=True)

    return ng_dict


def rescale_params(
    z: float,
    baseline_ng: dict[float, float],
    cfg: HodMockConfig = DEFAULT_CONFIG,
) -> dict:
    """指定した z の Ac/As を PFS ターゲット ng に合わせて rescale する。

    HODDIES では Ac と As が独立した絶対振幅なので、両者を同率でスケールする:
        scale    = ng_pfs / ng_hod
        Ac       = hod_ac_base[bin] * scale
        As       = hod_as_base[bin] * scale   （HODDIES YAML にそのまま書く値）

    z が hod_bins の範囲外の場合は端の bin に外挿する。

    Parameters
    ----------
    z : float
        対象スナップショットの赤方偏移。baseline_ng のキーに含まれている必要がある。
    baseline_ng : dict[float, float]
        compute_baseline_ng() が返す {z: ng} dict。
    cfg : HodMockConfig

    Returns
    -------
    dict
        キー: z, pfs_bin, hod_bin, ng_pfs, ng_hod, scale,
              Ac, As, density, extrapolated
    """
    if z not in baseline_ng:
        raise KeyError(f"z={z} は baseline_ng に含まれていない。compute_baseline_ng() を先に実行すること。")

    pfs_i = _find_bin_index(z, cfg.pfs_bins)
    hod_i = _find_bin_index(z, cfg.hod_bins, clamp=True)

    ng_pfs = float(cfg.pfs_ng[pfs_i])
    ng_hod = float(baseline_ng[z])
    scale  = ng_pfs / ng_hod

    Ac = float(cfg.hod_ac_base[hod_i]) * scale
    As = float(cfg.hod_as_base[hod_i]) * scale  # HODDIES 絶対振幅

    extrapolated = hod_i in (0, len(cfg.hod_bins) - 1) and not (
        cfg.hod_bins[hod_i][0] <= z <= cfg.hod_bins[hod_i][1]
    )

    return {
        "z":           z,
        "pfs_bin":     cfg.pfs_bins[pfs_i],
        "hod_bin":     cfg.hod_bins[hod_i],
        "ng_pfs":      ng_pfs,
        "ng_hod":      ng_hod,
        "scale":       scale,
        "Ac":          Ac,
        "As":          As,
        "density":     ng_pfs,
        "extrapolated": extrapolated,
    }


def compute_all_params(
    cfg: HodMockConfig = DEFAULT_CONFIG,
    verbose: bool = True,
) -> dict[float, dict]:
    """全スナップショットの rescaling を実行する。

    ① ベースライン ng を計算（compute_baseline_ng）
    ② 各 z で Ac/As を rescale（rescale_params）

    Parameters
    ----------
    cfg : HodMockConfig
    verbose : bool

    Returns
    -------
    dict[float, dict]
        {z: rescale_params の戻り値}
    """
    if verbose:
        print("[rescale] computing baseline ng ...", flush=True)
    baseline_ng = compute_baseline_ng(cfg, verbose=verbose)

    if verbose:
        print("[rescale] rescaling Ac/As ...", flush=True)
        print(f"  {'z':>5}  {'ng_pfs':>10}  {'ng_hod':>10}  {'scale':>8}  {'Ac':>8}  {'As':>8}", flush=True)

    results: dict[float, dict] = {}
    for z in cfg.z_list:
        r = rescale_params(z, baseline_ng, cfg)
        results[z] = r
        if verbose:
            tag = " [extrap]" if r["extrapolated"] else ""
            print(
                f"  {z:>5.2f}  {r['ng_pfs']:>10.3e}  {r['ng_hod']:>10.3e}"
                f"  {r['scale']:>8.4f}  {r['Ac']:>8.5f}  {r['As']:>8.5f}{tag}",
                flush=True,
            )

    return results


def save_params(
    params: dict[float, dict],
    cfg: HodMockConfig = DEFAULT_CONFIG,
) -> Path:
    """rescale 結果をコンパクト YAML に保存する（mock/maker.py が読む形式）。

    Parameters
    ----------
    params : dict[float, dict]
        compute_all_params() の戻り値。
    cfg : HodMockConfig

    Returns
    -------
    Path
        保存先のパス（cfg.hod_params_yaml）。
    """
    data = {
        "method": "rescale",
        "snapshots": {
            round(z, 2): {
                "Ac":      float(p["Ac"]),
                "As":      float(p["As"]),
                "density": float(p["density"]),
            }
            for z, p in params.items()
        },
    }

    cfg.hod_params_yaml.parent.mkdir(parents=True, exist_ok=True)
    with open(cfg.hod_params_yaml, "w") as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True)

    print(f"[rescale] saved → {cfg.hod_params_yaml}", flush=True)
    return cfg.hod_params_yaml
