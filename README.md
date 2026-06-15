# Symastic-data-analysis

Analisis data suhu **TCN4S (Autonics)** untuk mencari parameter **PID**, dikaitkan
dengan materi mata kuliah **Metode Numerik**.

## Isi

| Berkas | Keterangan |
|--------|-----------|
| `analisis_tcn4s.py` | Script analisis (diferensiasi, regresi, interpolasi, integrasi numerik) |
| `tuning_pid.py` | Estimasi Kp, Ki, Kd via Ziegler-Nichols & Cohen-Coon (model FOPDT) |
| `simulasi_pid.py` | Simulasi loop tertutup: kenapa PID lama gagal + PID koreksi terverifikasi |
| `data/data_suhu_tcn4s_14Juni2026.csv` | Data suhu TCN4S 14 Juni 2026 |
| `LAPORAN_ANALISIS.md` | Laporan lengkap + diagnosis kendala "K tidak bisa dicari" |
| `output/kurva_reaksi_tcn4s.png` | Grafik kurva reaksi & laju reaksi |

## Cara menjalankan

```bash
pip install pandas numpy matplotlib
python analisis_tcn4s.py
```

## Ringkasan temuan

Data berbentuk **ramp linear** (R² = 0.9998), monoton naik tanpa osilasi, dan PV
menembus setpoint 100 °C terus naik ke 158 °C tanpa pernah mendatar.

**Akibatnya gain proses `K` tidak bisa dihitung** karena (1) tidak ada kolom output
kontroler/MV dan (2) tidak ada nilai steady-state. Ini respons *open-loop*, bukan
respons PID tertutup. Solusi & langkah pengambilan data ulang ada di
[`LAPORAN_ANALISIS.md`](LAPORAN_ANALISIS.md).
