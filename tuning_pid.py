#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tuning PID dari Kurva Reaksi TCN4S — Ziegler-Nichols & Cohen-Coon
================================================================

Cara mendapatkan Kp, Ki, Kd dari data suhu yang berbentuk kurva reaksi,
memakai model FOPDT (First Order Plus Dead Time):

        G(s) = K * e^(-L*s) / (tau*s + 1)

PENTING — asumsi yang dipakai (wajib ditulis di laporan):
  * Data tidak punya kolom output/MV dan tidak punya steady-state asli.
  * Maka K TIDAK BISA DIUKUR, hanya bisa DIESTIMASI dengan mengasumsikan
    besar step input/output = 100 (mis. 100% daya heater atau rentang SV).
  * Dengan asumsi itu, 3 parameter diestimasi dari kurva:
        L   = dead time      ~ 1 interval sampling          = 10 s
        tau = time constant  ~ waktu PV mencapai setpoint    = 250 s
        K   = gain proses    = (PV_akhir - PV_awal)/step100  = 1.18

Hubungan gain paralel:  Ki = Kp/Ti ,  Kd = Kp*Td

Cara pakai:
    python tuning_pid.py
"""

import os
import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(__file__), "data", "data_suhu_tcn4s_14Juni2026.csv")


# ---------------------------------------------------------------------------
# Estimasi parameter FOPDT (K, tau, L) dari data + asumsi step input
# ---------------------------------------------------------------------------
def estimasi_fopdt(t, pv, step_input=100.0, sp=100.0):
    L = float(t[1] - t[0])                     # dead time ~ 1 sampling
    # tau ~ waktu PV memotong setpoint (interpolasi linear)
    tau = float(t[-1])
    for i in range(len(t) - 1):
        if pv[i] < sp <= pv[i + 1]:
            tau = t[i] + (sp - pv[i]) * (t[i + 1] - t[i]) / (pv[i + 1] - pv[i])
            break
    K = (pv[-1] - pv[0]) / step_input          # gain (ASUMSI step=step_input)
    return K, tau, L


# ---------------------------------------------------------------------------
# Rumus tuning (PID, bentuk paralel)
# ---------------------------------------------------------------------------
def ziegler_nichols(K, tau, L):
    Kp = 1.2 * tau / (K * L)
    Ti, Td = 2.0 * L, 0.5 * L
    return dict(Kp=Kp, Ki=Kp / Ti, Kd=Kp * Td, Ti=Ti, Td=Td)


def cohen_coon(K, tau, L):
    r = L / tau
    Kp = (1.0 / K) * (tau / L) * (4.0 / 3.0 + r / 4.0)
    Ti = L * (32.0 + 6.0 * r) / (13.0 + 8.0 * r)
    Td = L * 4.0 / (11.0 + 2.0 * r)
    return dict(Kp=Kp, Ki=Kp / Ti, Kd=Kp * Td, Ti=Ti, Td=Td)


def cetak(nama, g):
    print(f"{nama:16s} Kp={g['Kp']:7.2f}  Ki={g['Ki']:6.2f}  Kd={g['Kd']:8.2f}"
          f"   (Ti={g['Ti']:.1f}s, Td={g['Td']:.1f}s)")


def main():
    df = pd.read_csv(DATA)
    t = df["Waktu_detik"].astype(float).to_numpy()
    pv = df["PV_C"].astype(float).to_numpy()

    K, tau, L = estimasi_fopdt(t, pv, step_input=100.0, sp=100.0)

    print("=" * 64)
    print("TUNING PID DARI KURVA REAKSI TCN4S")
    print("=" * 64)
    print("Parameter FOPDT hasil estimasi (ASUMSI step input = 100):")
    print(f"   K (gain)        = {K:.3f}")
    print(f"   tau (time const)= {tau:.0f} s   (waktu PV capai SP)")
    print(f"   L (dead time)   = {L:.0f} s   (1 interval sampling)")
    print("-" * 64)
    cetak("Ziegler-Nichols", ziegler_nichols(K, tau, L))
    cetak("Cohen-Coon", cohen_coon(K, tau, L))
    print("=" * 64)
    print("CATATAN: K di sini ESTIMASI, bukan ukuran. Berlaku hanya bila")
    print("asumsi besar step input benar. Untuk hasil sahih: rekam kolom MV")
    print("(% output) dan biarkan suhu mencapai steady-state.")


if __name__ == "__main__":
    main()
