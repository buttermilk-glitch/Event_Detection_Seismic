# Literatur dan opsi lokasi: monitoring seismik 7 hari

Dibuat 23 Nov 2026. Metadata (judul, penulis, jurnal, DOI) dicocokkan dengan halaman penerbit/indeks hasil pencarian. Kolom "Temuan" hanya dari abstrak/ringkasan yang terbaca; teks lengkap belum dibaca. Baca teks lengkap sebelum dikutip di slide.

## A. Paper (jurnal)

| # | Penulis, tahun | Judul | Jurnal | DOI | Metode | Temuan (dari abstrak) | Langkah workflow |
|---|---|---|---|---|---|---|---|
| 1 | Allen, 1978 | Automatic earthquake recognition and timing from single traces | BSSA 68(5):1521-1532 | 10.1785/BSSA0680051521 | Program deteksi dan penentuan waktu otomatis pada satu trace (dasar STA/LTA dan logika pasca-deteksi) | Program komputer untuk deteksi dan timing gempa pada single trace | Deteksi STA/LTA, aturan penolakan false trigger |
| 2 | Earle & Shearer, 1994 | Characterization of global seismograms using an automatic-picking algorithm | BSSA 84(2):366-376 | 10.1785/BSSA0840020366 | STA/LTA pada fungsi envelope; keluaran waktu tiba + kualitas pick | Dibandingkan dengan hand pick pada seismogram yang sama | Characteristic function (envelope), picking P |
| 3 | Gibbons & Ringdal, 2006 | The detection of low magnitude seismic events using array-based waveform correlation | GJI 165(1):149-166 | 10.1111/j.1365-246X.2006.02865.x | Korelasi waveform dengan master event; stacking antar kanal | Korelasi satu kanal mendeteksi sinyal sekitar 0,7 unit magnitudo lebih kecil daripada STA/LTA pada array beam; array penuh menambah sekitar 0,4 | Metode tambahan untuk event kecil |
| 4 | Zhang & Wen, 2015 | An effective method for small event detection: match and locate (M&L) | GJI 200(3):1523-1537 | 10.1093/gji/ggu466 | Template + stacking korelasi dengan koreksi waktu tempuh relatif (scan 3-D) | Pada uji foreshock, 4x lebih banyak event daripada template dan 9% lebih banyak daripada matched filter | Template matching, lokasi relatif |
| 5 | Gibbons dkk., 2007 | The detection and location of low magnitude earthquakes in northern Norway using multi-channel waveform correlation at regional distances | Phys. Earth Planet. Inter. 160(3-4):285-309 | 10.1016/j.pepi.2006.11.008 | Korelasi multi-kanal jarak regional | Event M2,0-3,5 dengan SNR tinggi dipakai sebagai master | Syarat dan batas template matching |
| 6 | Chamberlain dkk., 2018 | EQcorrscan: Repeating and near-repeating earthquake detection and analysis in Python | SRL 89(1):173-181 | 10.1785/0220170151 | Paket Python matched-filter dan subspace, berbasis ObsPy | Implementasi terbuka, teruji | Implementasi template matching |
| 7 | Zhu & Beroza, 2019 | PhaseNet: a deep-neural-network-based seismic arrival-time picking method | GJI 216(1):261-273 | 10.1093/gji/ggy423 | U-Net; input 3 komponen; keluaran probabilitas P, S, noise; dilatih dengan pick NCEDC | Akurasi dan recall picking lebih tinggi dari metode yang dibandingkan | Picking P/S otomatis |
| 8 | Mousavi dkk., 2020 | Earthquake transformer: an attentive deep-learning model for simultaneous earthquake detection and phase picking | Nature Commun. 11:3952 | 10.1038/s41467-020-17591-w | Deep learning, deteksi + picking, data satu stasiun jarak lokal (<300 km) | Pada 5 minggu data kontinu Tottori, mendeteksi sekitar 2x lebih banyak event | Alternatif deteksi/picking |
| 9 | Zhang, Ellsworth & Beroza, 2019 | Rapid earthquake association and location | SRL 90(6):2276-2284 | 10.1785/0220190052 | Asosiasi fase + lokasi (delay-and-sum), dibandingkan katalog INGV | Pendekatan asosiasi dan lokasi otomatis | Lokasi dari pick multi-stasiun |
| 10 | Krischer dkk., 2015 | ObsPy: a bridge for seismology into the scientific Python ecosystem | Comput. Sci. Discov. 8(1):014003 | 10.1088/1749-4699/8/1/014003 | Pustaka Python untuk data seismik | Dasar I/O, klien FDSN, pemrosesan | Seluruh workflow |
| 11 | Beyreuther dkk., 2010 | ObsPy: a Python toolbox for seismology | SRL 81(3):530-533 | 10.1785/gssrl.81.3.530 | Pustaka Python | Pembuka ObsPy | Seluruh workflow |
| 12 | Turner, Latto & Reading, 2021 | An ObsPy library for event detection and seismic attribute calculation: preparing waveforms for automated analysis | J. Open Res. Softw. 9(1):29 | 10.5334/jors.365 | Pustaka ObsPy: unduh, preprocessing, deteksi STA/LTA, atribut event | Dirancang untuk sumber non-gempa (gletser, longsor); fungsi unduh/deteksi dapat dipakai | Pipeline unduh-deteksi |

Catatan jurnal: poin 10 dan 12 bukan jurnal seismologi utama (tetap sah sebagai rujukan perangkat lunak).

## B. Dokumentasi teknis (bukan paper jurnal)

| Sumber | DOI/tautan | Dipakai untuk |
|---|---|---|
| Trnkoczy, 2012, IS 8.1: Understanding and parameter setting of STA/LTA trigger algorithm, NMSOP-2 | 10.2312/GFZ.NMSOP-2_IS_8.1 | Memilih STA, LTA, trigger/detrigger |
| Küperkoch, Meier & Diehl, 2012, Ch. 16 Automated event and phase identification, NMSOP-2 | 10.2312/GFZ.NMSOP-2_ch16 | Deteksi, picker otomatis, preprocessing filter |
| ObsPy docs (`classic_sta_lta`, `trigger_onset`, `remove_response`) | https://docs.obspy.org | Implementasi |
| AusPass FDSN Guide | https://auspass.edu.au/xwiki/bin/view/Data/FDSN%20Guide/ | Cara unduh |

## C. Sumber data (bukan literatur)

AusPass (FDSN), katalog Geoscience Australia (salinan di AusPass, diperbarui bulanan, tanpa phase pick), EarthScope/IRIS, GeoNet.

## D. Celah literatur (belum ditutup)

Belum ada paper terverifikasi untuk: QC/analisis noise (PSD), penolakan false trigger secara kuantitatif selain Allen 1978 dan NMSOP, dan perbandingan katalog. Perlu dicari dan diverifikasi sebelum dimasukkan ke slide. Jangan mengutip dari ingatan.

## E. Opsi lokasi

Terverifikasi dari halaman AusPass: S1 (AUSIS), 54 stasiun, 2011-sekarang, permanen, akses terbuka; M8 (MATE), 38 stasiun, 2024-sekarang, permanen, terbuka. Jaringan AU **tidak lagi** di AusPass (ambil dari EarthScope/IRIS). Jaringan VW/VX/DU kini di server Univ. Melbourne.

| | Opsi A (disarankan) | Opsi B | Opsi C |
|---|---|---|---|
| Lokasi | Woods Point, Victoria, Australia | Sama, tetapi stasiun jaringan AU | Selandia Baru |
| Periode | 21-28 Sep 2021 (UTC), 7 hari | Sama | Dipilih dari katalog GeoNet |
| Sumber | AusPass, jaringan S1 | EarthScope/IRIS, jaringan AU | GeoNet FDSN (`Client("GEONET")`) |
| Seismisitas (terverifikasi) | Mw 5,9 pada 2021-09-21 23:15:53 UTC; GA mencatat 10 aftershock M2,4-4,1 dalam 24 jam pertama | Sama | Aktif; belum dicek per minggu |
| Stasiun | Contoh S1: AUGRF, AUMTS, AUMTC (muncul di panduan AusPass, kanal BHZ pada 21 Sep 2021). Koordinat, jarak, kanal, cakupan 7 hari **belum terverifikasi** | Belum diidentifikasi | Belum diidentifikasi |
| Risiko | Stasiun sekolah: noise antropogenik siang hari; aftershock kecil mungkin tak terlihat di jarak jauh | Kode stasiun dan respons belum dicek | Kurang sesuai prioritas Australia |
| Cakupan 7 hari | **Belum terverifikasi** (jalankan langkah `availability`) | Belum | Belum |

Status verifikasi: availability service tidak bisa dicek dari lingkungan kerja saya (jaringan shell diblokir, WebFetch menolak URL query yang dibuat sendiri). Pengecekan dilakukan lewat skrip `01_cek_dan_unduh_auspass.py` di komputer Anda.
