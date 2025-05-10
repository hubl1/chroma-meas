# 🐳 Chroma-Meas Docker Image for Lattice QCD

This Docker image provides a precompiled Chroma environment specialized for computing baryon and meson two-point (2pt) correlation functions in lattice Quantum Chromodynamics (QCD). All measurements in this image are performed using CPU. It supports loading precomputed propagators to generate 2-point functions, and can also generate propagators from gauge configurations.

---

## 📆 Docker Image Setup

Pull the Docker image:

```bash
docker pull ghcr.io/hubl1/chroma-meas:latest
```

Run the container, mounting your data directory:

```bash
docker run --rm -it -v /your/local/path/lqcd_data:/lqcd_data_path ghcr.io/hubl1/chroma-meas:latest bash
```

Ensure `/your/local/path/lqcd_data` contains necessary propagators and gauge configurations.

---

## 📂 Input Data Structure

Your `/lqcd_data_path` directory should be structured as:

```
lqcd_data/
├── precomputed_Propagators        # Precomputed propagators (light, strange, charm)
├── Propagators                    # Output directory (created automatically)
├── Configuration                  # Gauge configurations for propagator generation (CPU-based)
```

The precomputed propagators enable quick and direct calculation of baryon and meson two-point functions.

---

## 🚀 Running Chroma Measurements

Execute the measurements inside the Docker container using:

```bash
bash run.sh
```

* To measure baryon/meson two-point functions:

```bash
mpirun --allow-run-as-root -n 4 chroma -i ./baryon_2pt.xml
```

* To generate propagators from gauge configurations (CPU):

```bash
mpirun --allow-run-as-root -n 24 chroma -i ./cpu_inverter.xml
```

The provided `run.sh` script sets up directories and launches these measurements.

---
## 📦 Getting the Input Data

You can download a complete archive of `lqcd_data` (\~32 GB) with the following command:

```bash
curl -f -C - -o lqcd_data.tar.gz \
  https://ksefile.hpccube.com:65241/efile/share/L3B1YmxpYy9ob21lL2h1Ymw=/8iadVJCqq

# Then extract it

tar -xvzf lqcd_data.tar.gz
```

The unpacked structure looks like:

```
lqcd_data/
├── baryon_meson_data_h5           # Raw 2-point correlators (~25 GB)
├── precomputed_Propagators        # Precomputed propagators (~3.3 GB, for chroma-meas Docker image)
├── precomputed_pydata_eff_mass    # Effective mass fit results (~3.1 GB)
├── precomputed_pydata_global_fit  # Global fit data (~18 MB)
```

Set the path in `params/env.sh`:

```bash
export lqcd_data_path=/absolute/path/to/lqcd_data
```

---

## 🧚 Supported Baryons and Mesons

**Note:** All measurements are performed using CPUs in this image.

For **vector mesons** and **spin-3/2 baryons**, add `_x`, `_y`, or `_z` to the particle name to specify the polarization component.

### Light Sector

* `PROTON`
* `LAMBDA`
* `SIGMA`
* `XI`
* `DELTA`
* `SIGMA_STAR`
* `XI_STAR`
* `OMEGA`
* `PION`
* `ETA_S`

### Charm Sector

* `LAMBDA_C`
* `SIGMA_C`
* `XI_C`
* `XI_C_PRIME`
* `XI_CC`
* `OMEGA_C`
* `OMEGA_CC`
* `SIGMA_STAR_C`
* `XI_STAR_C`
* `OMEGA_STAR_C`
* `XI_STAR_CC`
* `OMEGA_STAR_CC`
* `OMEGA_CCC`
* `D`
* `D_S`
* `ETA_C`
* `D_STAR`
* `DS_STAR`
* `JPSI`
* `CHI_C0`
* `CHI_C1`
* `H_C`

Use naming conventions like `JPSI_x`, `OMEGA_STAR_CC_z` etc. for vector or spin-3/2 analysis.

---

🎯 **Note**: The Chroma-Meas Docker image streamlines lattice QCD analysis, from propagator generation to baryon and meson 2pt measurements, ensuring consistency and reproducibility on CPU-based systems.
