#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulasi & Koreksi PID untuk Proses Suhu TCN4S
==============================================

Kenapa PID Z-N/Cohen-Coon kemarin "salah" (stuck / naik terus)?
---------------------------------------------------------------
1. Pasangan (K=1.18, tau=250) TIDAK konsisten dengan kurva:
   model itu memprediksi laju 0.472 C/s, padahal data nyata hanya 0.248 C/s
   (meleset ~2x). Identifikasi modelnya salah -> PID ikut salah.
2. Kurva data berbentuk GARIS LURUS (ramp) => proses ini INTEGRATING,
   bukan FOPDT. Rumus Z-N/Cohen-Coon versi FOPDT tidak tepat di sini.
3. Pada model integrator murni (tanpa rugi panas, heater-only) begitu PV
   melewati setpoint, output=0 tapi heater tak bisa mendinginkan ->
   nyangkut di atas / divergen. Inilah gejala "stuck 77 / naik terus".

Pendekatan yang benar
---------------------
* Model proses dibuat REALISTIS (self-regulating, ada rugi panas) dan
  dikalibrasi agar laju awalnya cocok dengan data (~0.248 C/s).
* PID dirancang & DIVERIFIKASI lewat simulasi loop tertutup (overshoot,
  settling, steady-state error dihitung dengan integrasi/iterasi numerik).

ASUMSI (wajib ditulis di laporan; data tak cukup untuk menetapkannya):
  Tamb = 30 C, suhu maksimum pada output 100% ~ 200 C, dead time ~ 10 s.

Cara pakai:  python simulasi_pid.py
"""

import os
import numpy as np

OUTDIR = os.path.join(os.path.dirname(__file__), "output")

# ---- Parameter model proses (kalibrasi dari data) -------------------------
TAMB = 30.0          # suhu lingkungan (asumsi)
SLOPE0 = 0.248       # laju awal terukur dari data (C/s) saat heating penuh
TMAX_FULL = 200.0    # suhu maksimum pada output 100% (asumsi)
# dPV/dt = g*u - c*(PV - Tamb).  Kalibrasi:
#   g/c = (TMAX_FULL - Tamb)/100 ;  g*100 - c*(39-Tamb) = SLOPE0
_ratio = (TMAX_FULL - TAMB) / 100.0          # = g/c
C = SLOPE0 / (100.0 * _ratio - (39.0 - TAMB))
G = _ratio * C
L_PROC = 10.0        # dead time proses (s)


def simulasi(Kp, Ki, Kd, SP=100.0, Tstart=39.0, Tend=1500.0, dt=1.0,
             cooling=True):
    """PID diskrit (anti-windup + clamp output 0..100). Kembalikan t, PV, u."""
    n = int(Tend / dt)
    PV, I, e_prev = Tstart, 0.0, SP - Tstart
    dbuf = [0.0] * int(L_PROC / dt)
    ts, pvs, us = [], [], []
    for k in range(n):
        e = SP - PV
        I += e * dt
        d = (e - e_prev) / dt
        u = Kp * e + Ki * I + Kd * d
        if u > 100.0:            # clamp + anti-windup
            u = 100.0; I -= e * dt
        elif u < 0.0:
            u = 0.0; I -= e * dt
        e_prev = e
        u_del = dbuf.pop(0); dbuf.append(u)
        loss = C * (PV - TAMB) if cooling else 0.0
        PV += (G * u_del - loss) * dt
        ts.append(k * dt); pvs.append(PV); us.append(u)
    return np.array(ts), np.array(pvs), np.array(us)


def metrik(ts, pv, SP=100.0):
    final = pv[-1]
    sse = SP - final
    os = max(0.0, (pv.max() - SP) / SP * 100.0)
    band, settle = 0.02 * SP, None
    for i in range(len(pv)):
        if np.all(np.abs(pv[i:] - SP) <= band):
            settle = ts[i]; break
    return final, sse, os, settle


# Tuning yang dibandingkan ----------------------------------------------------
TUNINGS = {
    "Z-N (kemarin, SALAH)":      (25.41, 1.27, 127.05),
    "Cohen-Coon (kemarin)":      (28.46, 1.18, 102.76),
    "PID KOREKSI (disarankan)":  (12.0, 12.0 / 150.0, 12.0 * 6.0),  # Kp=12,Ti=150,Td=6
}


def main():
    print("=" * 70)
    print("SIMULASI & KOREKSI PID — PROSES SUHU TCN4S")
    print("=" * 70)
    print(f"Model proses (kalibrasi dari data): dPV/dt = {G:.5f}*u - {C:.6f}*(PV-{TAMB:.0f})")
    print(f"  tau_proc = {1/C:.0f}s, K_proc = {G/C:.2f} C/%, dead time = {L_PROC:.0f}s")
    print(f"  PID koreksi: Kp=12, Ti=150s (Ki=0.080), Td=6s (Kd=72)")
    print("-" * 70)
    print("%-26s %8s %8s %7s %9s" % ("Tuning", "akhir", "SS-err", "OS%", "settle"))
    print("-" * 70)
    hasil = {}
    for nm, (Kp, Ki, Kd) in TUNINGS.items():
        ts, pv, u = simulasi(Kp, Ki, Kd)
        f, sse, os, st = metrik(ts, pv)
        hasil[nm] = (ts, pv, u)
        print("%-26s %8.1f %8.2f %7.1f %9s"
              % (nm, f, sse, os, (f"{st:.0f}s" if st else "tak settle")))
    print("=" * 70)
    print("Catatan: model di atas pakai ASUMSI (Tamb, suhu maks, dead time).")
    print("Untuk PID yang benar-benar sahih: rekam kolom output %MV + biarkan")
    print("suhu mencapai steady-state, ATAU pakai Auto-Tuning (AT) bawaan TCN4S.")
    buat_grafik(hasil)


def buat_grafik(hasil):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(OUTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhline(100, ls="--", color="k", lw=1, label="Setpoint 100°C")
    for nm, (ts, pv, u) in hasil.items():
        ax.plot(ts, pv, label=nm)
    ax.set_xlabel("waktu (s)"); ax.set_ylabel("suhu PV (°C)")
    ax.set_title("Respons loop tertutup PID (model proses realistis)")
    ax.legend(); ax.grid(alpha=0.3); ax.set_ylim(30, 130)
    out = os.path.join(OUTDIR, "respons_pid.png")
    fig.tight_layout(); fig.savefig(out, dpi=120)
    print(f"[grafik] -> {out}")


if __name__ == "__main__":
    main()
