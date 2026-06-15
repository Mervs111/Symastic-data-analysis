#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analisis Data Suhu TCN4S — Pencarian Parameter PID
==================================================

Tujuan
------
1. Mencoba menentukan parameter PID dari kurva respons suhu TCN4S.
2. Mengaitkan setiap langkah dengan materi mata kuliah METODE NUMERIK:
     - Diferensiasi numerik (beda hingga) ............ cari laju reaksi / kurva reaksi
     - Regresi kuadrat terkecil (least squares) ...... fit model garis/FOPDT
     - Interpolasi linear ............................ cari waktu lintas SP & dead time
     - Integrasi numerik (trapesium) ................. hitung error integral (IAE/ISE)
3. Mendiagnosis KENDALA: kenapa "K" (gain proses) tidak bisa dicari dan
   kenapa data terlihat "tidak stabil".

Cara pakai
----------
    python analisis_tcn4s.py

Output
------
    - Ringkasan analisis di terminal.
    - Grafik di output/kurva_reaksi_tcn4s.png

Catatan: TCN4S (Autonics) memakai istilah Proportional Band (P, dalam %),
Integral time (I, detik), Derivative time (D, detik). Hubungan ke gain:
    Kp = 100 / P(%)
"""

import os
import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(__file__), "data", "data_suhu_tcn4s_14Juni2026.csv")
OUTDIR = os.path.join(os.path.dirname(__file__), "output")


def muat_data(path=DATA):
    """Muat data dan kembalikan t (detik), PV (C), SV (C)."""
    df = pd.read_csv(path)
    t = df["Waktu_detik"].astype(float).to_numpy()
    pv = df["PV_C"].astype(float).to_numpy()
    sv = df["SV_C"].astype(float).to_numpy()
    return df, t, pv, sv


# ---------------------------------------------------------------------------
# METODE NUMERIK 1 — Diferensiasi numerik (beda hingga pusat)
#   Turunan PV terhadap waktu = "laju reaksi" proses. Dipakai untuk menemukan
#   titik infleksi (slope maksimum) pada metode kurva reaksi Ziegler-Nichols.
# ---------------------------------------------------------------------------
def laju_reaksi(t, pv):
    dpvdt = np.gradient(pv, t)            # beda hingga pusat (np.gradient)
    i_inf = int(np.argmax(dpvdt))        # titik slope maksimum
    return dpvdt, i_inf


# ---------------------------------------------------------------------------
# METODE NUMERIK 2 — Regresi kuadrat terkecil (least squares)
#   Mencocokkan garis PV = m*t + c. Jika R^2 ~ 1, respons = RAMP linear
#   (ciri proses integrating / open-loop), bukan kurva-S orde-satu.
# ---------------------------------------------------------------------------
def regresi_linear(t, pv):
    A = np.vstack([t, np.ones_like(t)]).T
    (m, c), *_ = np.linalg.lstsq(A, pv, rcond=None)
    pred = m * t + c
    ss_res = float(np.sum((pv - pred) ** 2))
    ss_tot = float(np.sum((pv - pv.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot
    return m, c, r2


# ---------------------------------------------------------------------------
# METODE NUMERIK 3 — Interpolasi linear
#   Cari waktu saat PV memotong SV (lewat dua titik yang mengapit).
# ---------------------------------------------------------------------------
def waktu_lintas_sp(t, pv, sv):
    for i in range(len(t) - 1):
        if (pv[i] - sv[i]) * (pv[i + 1] - sv[i + 1]) <= 0 and pv[i] != sv[i]:
            return t[i] + (sv[i] - pv[i]) * (t[i + 1] - t[i]) / (pv[i + 1] - pv[i])
    return None


# ---------------------------------------------------------------------------
# METODE NUMERIK 4 — Integrasi numerik (aturan trapesium)
#   Error integral IAE = ∫|SV-PV| dt, ISE = ∫(SV-PV)^2 dt. Indikator performa.
# ---------------------------------------------------------------------------
def error_integral(t, pv, sv):
    e = sv - pv
    iae = float(np.trapezoid(np.abs(e), t))
    ise = float(np.trapezoid(e ** 2, t))
    return iae, ise


def cek_osilasi(dpvdt):
    """Jumlah pergantian arah turunan; 0 = monoton (tak ada osilasi)."""
    return int(np.sum(np.diff(np.sign(dpvdt)) != 0))


def analisis():
    df, t, pv, sv = muat_data()
    dpvdt, i_inf = laju_reaksi(t, pv)
    m, c, r2 = regresi_linear(t, pv)
    tc = waktu_lintas_sp(t, pv, sv)
    iae, ise = error_integral(t, pv, sv)
    osc = cek_osilasi(dpvdt)

    # Dead time L: ekstrapolasi garis fit ke baseline awal PV[0]
    L = (pv[0] - c) / m

    print("=" * 70)
    print("ANALISIS DATA SUHU TCN4S — PARAMETER PID")
    print("=" * 70)
    print(f"Jumlah titik         : {len(t)}  (t = {t[0]:.0f}..{t[-1]:.0f} s, dt = {t[1]-t[0]:.0f} s)")
    print(f"SV (setpoint)        : {np.unique(sv)} C  (konstan)")
    print(f"PV awal / akhir / maks: {pv[0]:.0f} / {pv[-1]:.0f} / {pv.max():.0f} C")
    print("-" * 70)
    print("[1] Diferensiasi numerik (laju reaksi dPV/dt):")
    print(f"    rata-rata = {dpvdt.mean():.4f} C/s, std = {dpvdt.std():.4f} (hampir konstan)")
    print(f"    slope maks (infleksi) di t={t[i_inf]:.0f}s, PV={pv[i_inf]:.0f}C, R={dpvdt[i_inf]:.4f} C/s")
    print("[2] Regresi least-squares:")
    print(f"    PV = {m:.5f}*t + {c:.3f}   R^2 = {r2:.6f}  -> RAMP linear (integrating)")
    print(f"[3] Interpolasi: PV memotong SV(100C) pada t = {tc:.1f} s")
    print(f"    Dead time L (ekstrapolasi ke baseline {pv[0]:.0f}C) = {L:.2f} s")
    print(f"[4] Integrasi numerik: IAE = {iae:.1f} C.s, ISE = {ise:.1f} C^2.s")
    print(f"    Pergantian arah (osilasi) = {osc}  -> {'TIDAK ADA osilasi (monoton)' if osc == 0 else 'ada osilasi'}")
    print("-" * 70)

    # --- DIAGNOSIS KENDALA -------------------------------------------------
    print("DIAGNOSIS: KENAPA 'K' TIDAK BISA DICARI")
    print("-" * 70)
    print(
        "K = gain proses = (ΔPV steady-state) / (ΔMV output kontroler).\n"
        "Untuk mendapatkannya WAJIB ada DUA hal yang TIDAK ADA pada data ini:\n"
        "  (a) Kolom OUTPUT/MV (% heater) — di file hanya ada PV & SV.\n"
        "  (b) Nilai AKHIR (steady-state) — PV tak pernah mendatar (R^2~1, ramp).\n"
        "Karena PV menembus SP dan terus naik (runaway, overshoot "
        f"{(pv.max()-100)/100*100:.0f}%), tak ada plateau untuk mengukur ΔPV.\n"
        "Tanpa osilasi pula, metode Ultimate Gain (Ku, Pu) juga tak bisa dipakai.\n"
        "=> 'K tidak bisa dicari' & 'data tidak stabil' adalah BENAR: ini data\n"
        "   open-loop/manual (kontrol tidak menutup), bukan respons PID tertutup."
    )
    print("-" * 70)

    # --- JIKA OUTPUT/MV DIKETAHUI ------------------------------------------
    print("SOLUSI (perlu data tambahan output %MV) — contoh Z-N integrating:")
    dMV = 100.0  # ASUMSI: step output 100% (full power). GANTI sesuai data nyata.
    Rn = m / dMV
    if L > 0:
        Kp = 1.2 / (Rn * L)
        print(f"    Asumsi step MV={dMV:.0f}% -> R_norm={Rn:.5f} C/(s·%), L={L:.2f}s")
        print(f"    PID Z-N: Kp={Kp:.2f}, Ti={2*L:.1f}s, Td={0.5*L:.1f}s")
        print(f"    (TCN4S: P-band = {100/Kp:.2f}%, I = {2*L:.0f}s, D = {0.5*L:.0f}s)")
    print("    *Angka di atas hanya valid bila ΔMV benar; ganti dMV dengan nilai nyata.")
    print("=" * 70)

    buat_grafik(t, pv, sv, dpvdt, i_inf, m, c, L)
    return dict(m=m, c=c, r2=r2, R=dpvdt[i_inf], L=L, tc=tc, iae=iae, ise=ise, osc=osc)


def buat_grafik(t, pv, sv, dpvdt, i_inf, m, c, L):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(OUTDIR, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    # Kiri: kurva reaksi + garis fit + tangen + dead time
    ax[0].plot(t, pv, "o-", ms=4, label="PV (suhu aktual)")
    ax[0].axhline(sv[0], ls="--", color="r", label=f"SV = {sv[0]:.0f} C")
    ax[0].plot(t, m * t + c, "g--", lw=1.5, label=f"Fit: {m:.3f}t+{c:.1f} (R²~1)")
    ax[0].plot(t[i_inf], pv[i_inf], "ks", ms=8, label="titik infleksi")
    ax[0].axvline(L, color="purple", ls=":", label=f"dead time L≈{L:.1f}s")
    ax[0].set_xlabel("waktu (s)"); ax[0].set_ylabel("suhu (°C)")
    ax[0].set_title("Kurva Reaksi TCN4S — RAMP (tak ada steady-state)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    # Kanan: laju reaksi dPV/dt
    ax[1].plot(t, dpvdt, "o-", ms=4, color="darkorange")
    ax[1].axhline(dpvdt.mean(), ls="--", color="gray", label=f"rata² {dpvdt.mean():.3f} C/s")
    ax[1].set_xlabel("waktu (s)"); ax[1].set_ylabel("dPV/dt (°C/s)")
    ax[1].set_title("Diferensiasi numerik: laju ~konstan → integrating")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    fig.suptitle("Analisis Numerik Data Suhu TCN4S (14 Juni 2026)", fontweight="bold")
    fig.tight_layout()
    out = os.path.join(OUTDIR, "kurva_reaksi_tcn4s.png")
    fig.savefig(out, dpi=120)
    print(f"[grafik] disimpan -> {out}")


if __name__ == "__main__":
    analisis()
