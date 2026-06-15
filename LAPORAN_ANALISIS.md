# Laporan Analisis Data Suhu TCN4S — Pencarian Parameter PID

**Sumber data:** `data/data_suhu_tcn4s_14Juni2026.csv` (TCN4S Autonics, 14 Juni 2026)
**Script:** `analisis_tcn4s.py` · **Grafik:** `output/kurva_reaksi_tcn4s.png`

---

## 1. Ringkasan masalah (kendala)

> "PID-nya mau dicari, tapi **K-nya nggak bisa dicari** dan **datanya nggak stabil**."

Setelah dianalisis, **keluhan itu benar — dan itu BUKAN salah hitung, tapi keterbatasan data.**

---

## 2. Apa isi datanya

| Kolom | Isi |
|------|-----|
| Waktu | 0 – 490 detik, interval 10 s (50 titik) |
| PV (suhu aktual) | 39 °C → naik terus → **158 °C** |
| SV (setpoint) | **100 °C konstan** |
| Fase | keterangan teks (Naik Perlahan, Stabil, dll.) |

**Yang TIDAK ada:** kolom **output kontroler / MV** (persen daya heater). Ini krusial (lihat §4).

---

## 3. Hasil analisis numerik

| Besaran | Nilai | Cara (metode numerik) |
|--------|-------|------------------------|
| Laju reaksi rata-rata dPV/dt | **0.242 °C/s** (std 0.018) | Diferensiasi numerik (beda hingga pusat) |
| Model garis | **PV = 0.248·t + 38.1**, **R² = 0.9998** | Regresi kuadrat terkecil |
| Titik infleksi / slope maks | t=30 s, R=0.25 °C/s | argmax dari turunan |
| PV memotong SP 100 °C | t ≈ **250 s** | Interpolasi linear |
| Dead time L (ekstrapolasi) | ≈ 3.5 s | Interpolasi |
| Jumlah osilasi | **0** (monoton naik) | Hitung pergantian tanda turunan |
| IAE / ISE | 14945 / 606375 | Integrasi numerik (trapesium) |

**Kesimpulan bentuk data:** R² ≈ 1 → kurvanya **garis lurus (RAMP)**, bukan kurva-S
orde-satu. PV menembus SP dan terus naik (overshoot 58%, *runaway*), **tidak pernah
mendatar di nilai akhir**. Tidak ada osilasi.

---

## 4. Kenapa "K" tidak bisa dicari

`K` (gain proses statis) didefinisikan:

```
K = ΔPV (steady-state) / ΔMV (perubahan output kontroler)
```

Untuk menghitungnya **wajib** ada dua hal yang **tidak ada** di data ini:

1. **Output / MV (% heater).** File hanya berisi PV dan SV. Penyebut `ΔMV` tidak diketahui.
2. **Nilai akhir (steady-state).** PV berbentuk ramp linear yang tak pernah mendatar,
   jadi `ΔPV steady-state` tidak terdefinisi.

Karena itu **K memang tidak bisa dihitung** dari data ini — bukan karena rumus salah.

Metode alternatif juga buntu dengan data ini:
- **Ziegler-Nichols open-loop (kurva reaksi / FOPDT):** butuh plateau + ΔMV → tidak ada.
- **Ziegler-Nichols closed-loop (Ultimate Gain Ku, Pu):** butuh **osilasi berkelanjutan**
  → data 100% monoton, tidak ada osilasi.

**Diagnosis akhir:** SV = 100 °C tapi PV cuek menembus ke 158 °C tanpa ada aksi koreksi.
Artinya ini **respons open-loop / manual (loop kontrol tidak menutup)**, bukan respons
PID tertutup. Maka "datanya nggak stabil" itu **wajar dan benar**.

---

## 5. Apa yang harus dilakukan (solusi)

Supaya K (dan PID) bisa didapat, ambil ulang datanya dengan benar — pilih salah satu:

**A. Metode kurva reaksi (open-loop) — paling cocok untuk tugas metode numerik**
1. Set TCN4S ke **mode manual**, beri **step output tetap** (mis. dari 0% → 50%),
   dan **catat nilai % output-nya** (tambahkan kolom MV).
2. Tunggu sampai suhu **mendatar (steady-state)**.
3. Dari kurva: `K = ΔPV/ΔMV`, `L` = dead time, `τ` = konstanta waktu.
4. Z-N FOPDT (PID): `Kp = 1.2·τ/(K·L)`, `Ti = 2L`, `Td = 0.5L`.

**B. Metode Ultimate Gain (closed-loop):** naikkan gain P sampai suhu **berosilasi tetap**,
catat `Ku` dan periode `Pu`, lalu pakai tabel Z-N.

**C. Praktis:** pakai fitur **Auto-Tuning (AT)** bawaan TCN4S — kontroler menghitung
P, I, D otomatis.

> Catatan istilah: TCN4S memakai **Proportional Band P (%)**, bukan gain langsung.
> Hubungannya `Kp = 100 / P(%)`. Satuan I dan D dalam detik.

Jika nanti sudah punya ΔMV, script `analisis_tcn4s.py` sudah menyiapkan rumus Z-N
untuk proses *integrating* — tinggal ganti variabel `dMV` dengan nilai output sebenarnya
(saat ini diasumsikan 100% hanya sebagai contoh).

---

## 5b. Mendapatkan Kp, Ki, Kd lewat ESTIMASI (cara reaction curve) — `tuning_pid.py`

K **tidak bisa diukur** dari data ini, tetapi **bisa diestimasi** kalau kita berani
mengambil **asumsi** besar step input/output. Inilah yang dipakai bila tetap ingin
angka PID dari kurva reaksi (model FOPDT `G(s)=K·e^(−Ls)/(τs+1)`).

**Parameter FOPDT hasil estimasi (asumsi step input = 100):**

| Parameter | Nilai | Asal dari data |
|-----------|-------|-----------------|
| `K` (gain) | ≈ **1.18–1.19** | (PV_akhir − PV_awal)/100 = (158 − 40)/100 |
| `τ` (time constant) | **250 s** | waktu PV mencapai setpoint 100 °C (interpolasi) |
| `L` (dead time) | **10 s** | 1 interval sampling |

**Rumus & hasil:**

| Metode | Rumus Kp | Kp | Ki=Kp/Ti | Kd=Kp·Td |
|--------|----------|----|----------|----------|
| **Ziegler-Nichols** | `1.2·τ/(K·L)`, Ti=2L, Td=0.5L | ≈ **25.4** | ≈ **1.27** | ≈ **127** |
| **Cohen-Coon** | `(1/K)(τ/L)(4/3+L/4τ)` | ≈ **28.5** | ≈ **1.18** | ≈ **103** |

> Catatan akurasi: pakai `PV_awal = 40` → `K = 1.18` → Kp Z-N = 25.41 (sama persis
> hasil umum di kelas). Pakai data mentah `PV_awal = 39` → `K = 1.19` → Kp = 25.21.
> Selisihnya hanya pembulatan suhu awal; metodenya identik.

**⚠️ Wajib ditulis di laporan:** nilai K ini **estimasi, bukan ukuran**. Sah hanya
bila asumsi step input benar. Untuk PID yang benar-benar sahih, tetap perlu kolom
MV dan steady-state asli (lihat §5). Jalankan: `python tuning_pid.py`.

---

## 5c. Kenapa PID itu "salah" saat diuji, dan koreksinya — `simulasi_pid.py`

Saat di-coba (Z-N **nyangkut**, Cohen-Coon **naik terus nggak ngaruh**), itu **gejala
model salah**, bukan sekadar salah angka:

1. **Model tak konsisten dengan kurva.** Pasangan (K=1.18, τ=250) memprediksi laju
   `K·step/τ = 0.472 °C/s`, padahal laju nyata data hanya **0.248 °C/s** (meleset ~2×).
2. **Prosesnya INTEGRATING, bukan FOPDT.** Kurva data = garis lurus → τ tak terdefinisi,
   sehingga rumus Z-N/Cohen-Coon versi FOPDT memang **tidak tepat** dipakai.
3. **Diuji di plant tak fisis.** Pada model integrator murni (tanpa rugi panas,
   heater-only), begitu PV lewat 100 °C output=0 tapi heater **tak bisa mendinginkan**
   → PV nyangkut di atas / divergen. Itulah "stuck" & "naik terus".

**Koreksi (diverifikasi lewat simulasi loop tertutup):** model proses dibuat realistis
(self-regulating, ada rugi panas, dikalibrasi agar laju awal = 0.248 °C/s; asumsi
Tamb=30 °C, suhu maks pada 100% ≈ 200 °C, dead time 10 s), lalu PID diuji:

| Tuning | PV akhir | Overshoot | Settling (±2%) |
|--------|----------|-----------|-----------------|
| Z-N kemarin (Kp 25.4) | 100 °C | 1.2% | ~306 s |
| Cohen-Coon kemarin | 100 °C | 0.6% | ~306 s |
| **PID KOREKSI: Kp=12, Ki=0.08, Kd=72** (Ti=150s, Td=6s) | **100 °C** | **0%** | ~339 s |

PID koreksi: **overshoot 0%** di model realistis dan **paling kecil (1.2%)** bahkan di
kasus terburuk (integrator murni). Grafiknya: `output/respons_pid.png`.

> **Kesimpulan jujur:** angka PID yang benar-benar tepat **tidak bisa dipastikan** dari
> data ini karena tak ada kolom MV & steady-state — semua nilai bergantung asumsi model.
> Solusi sahih: (1) rekam **% output (MV)** + tunggu suhu **mendatar** lalu identifikasi
> ulang, atau (2) pakai **Auto-Tuning (AT)** bawaan TCN4S. Jalankan: `python simulasi_pid.py`.

---

## 6. Kaitan dengan mata kuliah Metode Numerik

| Materi metode numerik | Penerapan di analisis ini |
|-----------------------|----------------------------|
| **Diferensiasi numerik** (beda hingga maju/pusat) | Hitung laju reaksi dPV/dt, cari titik infleksi |
| **Regresi kuadrat terkecil** | Fit model PV = m·t + c, ukur R² → deteksi bentuk ramp |
| **Interpolasi linear** | Cari waktu PV = SV dan estimasi dead time L |
| **Integrasi numerik** (trapesium) | Hitung error integral IAE = ∫\|e\|dt, ISE = ∫e²dt |
| **Pencocokan model / akar persamaan** | Estimasi parameter FOPDT (K, τ, L) saat data lengkap |

Jadi pencarian PID dari kurva suhu adalah **kasus terapan langsung** dari metode numerik:
turunan, regresi, interpolasi, dan integrasi numerik dipakai berurutan untuk
mengekstrak parameter dinamika proses.
