"""
config.py — hodmock パイプライン全体の設定値を一元管理する。

使い方:
    # デフォルト（Uchuu + PFS 設定）
    from hodmock.config import HodMockConfig
    cfg = HodMockConfig()

    # z_list だけ別のカタログ用に変える
    cfg = HodMockConfig(z_list=[0.5, 1.0, 1.5])

    # 一部だけ上書き（dataclasses.replace を使う）
    from dataclasses import replace
    cfg2 = replace(cfg, nthreads=8, mock_outdir=Path("/tmp/mock"))
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# デフォルト値を返すヘルパー関数
# （dataclass の mutable default は field(default_factory=...) が必要なため）
# ─────────────────────────────────────────────────────────────────────────────

def _default_z_list() -> list[float]:
    """Uchuu スナップショット赤方偏移リスト（デフォルト）。"""
    return [
        0.63, 0.70, 0.78, 0.86, 0.94, 1.03,
        1.12, 1.22, 1.32, 1.43, 1.54, 1.65,
        1.77, 1.90, 2.03, 2.17, 2.31,
    ]


def _default_z_to_halodir() -> dict[float, str]:
    """z → ハローカタログサブディレクトリ名の対応表（Uchuu デフォルト）。"""
    return {
        0.63: "halodir_038",
        0.70: "halodir_037",
        0.78: "halodir_036",
        0.86: "halodir_035",
        0.94: "halodir_034",
        1.03: "halodir_033",
        1.12: "halodir_032",
        1.22: "halodir_031",
        1.32: "halodir_030",
        1.43: "halodir_029",
        1.54: "halodir_028",
        1.65: "halodir_027",
        1.77: "halodir_026",
        1.90: "halodir_025",
        2.03: "halodir_024",
        2.17: "halodir_023",
        2.31: "halodir_022",
    }


def _default_pfs_bins() -> list[tuple[float, float]]:
    """PFS ELG z ビン定義（デフォルト）。"""
    return [
        (0.6, 0.8),
        (0.8, 1.0),
        (1.0, 1.2),
        (1.2, 1.4),
        (1.4, 1.6),
        (1.6, 2.0),
        (2.0, 2.4),
    ]


def _default_pfs_ng() -> np.ndarray:
    """PFS ELG ターゲット銀河数密度 [(Mpc/h)^-3]（デフォルト）。
    pfs_bins と同じ順序で並べること。
    """
    return np.array([
        5.4e-4,
        11.2e-4,
        10.6e-4,
        14.4e-4,
        10.8e-4,
        4.3e-4,
        4.8e-4,
    ])


def _default_hod_bins() -> list[tuple[float, float]]:
    """rescaling 用 HOD パラメータ bin 定義（デフォルト）。"""
    return [
        (0.8, 1.1),  # z1 bin
        (1.1, 1.6),  # z2 bin
    ]


def _default_hod_baseline() -> list[dict]:
    """hodmock.halomodel での ng 計算に使うベースライン HOD パラメータ。
    hod_bins の各 bin に対応する dict のリスト。Ac=0.1 のときの ng が基準値 ng_hod になる。
    """
    base = {
        "model":      "mhmq",
        "logMc":      11.62,
        "sigma_M":    0.21,
        "Ac":         0.1,
        "gamma":      6.49,
        "logM0":      11.10,
        "logM1":      13.0,
        "alpha":      0.74,
        "As":         0.38,
        "poff":       0.0,
        "Roff":       0.0,
        "conformity": True,
        "sat_model":  "elg",
    }
    return [dict(base), dict(base)]  # z1 / z2 bin（必要なら個別に変更可）


def _default_hod_ac_base() -> np.ndarray:
    return np.array([0.1, 0.1])


def _default_hod_as_base() -> np.ndarray:
    """HODDIES ベースライン As 値（hod_bins に対応）。
    HODDIES では Ac と As は独立した絶対振幅。ng に合わせて両者を同率でスケールする:
        Ac_scaled = hod_ac_base * scale
        As_scaled = hod_as_base * scale  （HODDIES YAML にそのまま書く値）
    """
    return np.array([0.38, 0.47])


def _default_k_arr() -> np.ndarray:
    return np.logspace(-4, 3, 1000)  # [h/Mpc]


def _default_m_arr() -> np.ndarray:
    return np.logspace(10, 15, 200)  # [Msun/h]


def _default_cosmo_params() -> dict:
    """Planck 2018 コスモロジーパラメータ（CAMB 形式）。"""
    return {
        'H0':    67.4,
        'ombh2': 0.0224,
        'omch2': 0.12,
        'As':    2.1e-9,
        'ns':    0.965,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Config dataclass
# ─────────────────────────────────────────────────────────────────────────────

_HERE = Path(__file__).parent  # hodmock/ ディレクトリ


@dataclass
class HodMockConfig:
    """hodmock パイプライン全体の設定。

    全フィールドにデフォルト値があるので、変えたいフィールドだけ指定すればよい:
        cfg = HodMockConfig(z_list=[0.5, 1.0], nthreads=8)

    dataclasses.replace() で一部上書きしたコピーも作れる:
        from dataclasses import replace
        cfg2 = replace(cfg, mock_outdir=Path("/tmp/out"))
    """

    # ── スナップショット設定 ──────────────────────────────────────────────────
    z_list: list[float] = field(default_factory=_default_z_list)
    """対象とするスナップショットの赤方偏移リスト。
    別のハローカタログを使うときはここを変える。"""

    z_to_halodir: dict[float, str] = field(default_factory=_default_z_to_halodir)
    """z → ハローカタログサブディレクトリ名。halobase / z_to_halodir[z] がフルパス。"""

    # ── サーベイターゲット ng ─────────────────────────────────────────────────
    pfs_bins: list[tuple[float, float]] = field(default_factory=_default_pfs_bins)
    """サーベイ z ビン定義 [(zmin, zmax), ...]。"""

    pfs_ng: np.ndarray = field(default_factory=_default_pfs_ng)
    """ターゲット銀河数密度 [(Mpc/h)^-3]。pfs_bins と同順。"""

    # ── HOD ベースライン（rescaling 用）───────────────────────────────────────
    hod_bins: list[tuple[float, float]] = field(default_factory=_default_hod_bins)
    """rescaling 用 HOD パラメータ bin。範囲外の z は端の bin に外挿する。"""

    hod_baseline: list[dict] = field(default_factory=_default_hod_baseline)
    """hod_bins 各 bin のベースライン HOD パラメータ（hodmock.halomodel 形式）。"""

    hod_ac_base: np.ndarray = field(default_factory=_default_hod_ac_base)
    """hod_bins 各 bin の Ac ベースライン値（ng_hod 計算の基準値）。"""

    hod_as_base: np.ndarray = field(default_factory=_default_hod_as_base)
    """hod_bins 各 bin の HODDIES ベースライン As 値（絶対振幅）。
    halomodel での検証時は hod_ac_base * hod_as_base を渡す。"""

    # ── Halo model 計算グリッド ───────────────────────────────────────────────
    k_arr: np.ndarray = field(default_factory=_default_k_arr)
    """波数グリッド [h/Mpc]。"""

    m_arr: np.ndarray = field(default_factory=_default_m_arr)
    """ハロー質量グリッド [Msun/h]。"""

    # ── コスモロジー / Halo model モデル選択 ──────────────────────────────────
    cosmo_params: dict = field(default_factory=_default_cosmo_params)
    """CAMB に渡すコスモロジーパラメータ。キー: H0, ombh2, omch2, As, ns。"""

    colossus_cosmo: str = "planck18"
    """Colossus のコスモロジー名。cosmo_params と整合させること。"""

    hmf_model: str = "tinker08"
    """ハロー質量関数モデル（Colossus）。例: 'press74', 'sheth99', 'tinker08'。"""

    conc_model: str = "diemer19"
    """ハロー集中度モデル（Colossus）。例: 'duffy08', 'diemer19'。"""

    bias_model: str = "tinker10"
    """ハローバイアスモデル（Colossus）。例: 'tinker10', 'sheth01'。"""

    halomodel_Nr: int = 512
    """NFW プロファイルの Fourier 変換で使う動径グリッド点数。"""

    # ── シミュレーション / ジョブ設定 ────────────────────────────────────────
    mass_cut: float = 11.0
    """ハロー質量下限カット log10(M200c / [Msun/h])。"""

    nthreads: int = 27
    """ハローカタログ読み込みの並列スレッド数。"""

    seed: int = 42
    tracer: str = "ELG"

    # ── パス設定 ─────────────────────────────────────────────────────────────
    halobase: Path = field(default=Path("/data/PFS/Uchuu/RockstarExtendedM200c1e11"))
    """ハローカタログのベースディレクトリ。z_to_halodir と組み合わせてフルパスになる。"""

    mock_outdir: Path = field(default=Path("/home/honke/data/HOD_mock"))
    """生成した mock カタログの出力先ディレクトリ。"""

    template_yaml: Path = field(default=_HERE / "mock" / "template" / "Uchuu_ELG_template.yaml")
    """HODDIES に渡す YAML テンプレートのパス。"""

    hod_params_yaml: Path = field(default=_HERE / "hod_params_rescale.yaml")
    """rescaling が出力するコンパクト HOD パラメータ YAML のパス。"""

    # ── ユーティリティ ────────────────────────────────────────────────────────
    def halodir(self, z: float) -> Path:
        """指定した z のハローカタログディレクトリのフルパスを返す。"""
        name = self.z_to_halodir.get(z)
        if name is None:
            raise KeyError(f"z={z} は z_to_halodir に登録されていない")
        return self.halobase / name

    def mock_output(self, z: float) -> Path:
        """指定した z の mock カタログ出力パスを返す。"""
        return self.mock_outdir / f"HODmock_{self.tracer}_z{z:.2f}.fits"


# デフォルト設定インスタンス（そのまま import して使える）
DEFAULT_CONFIG = HodMockConfig()
