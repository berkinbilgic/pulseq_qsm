# pulseq_qsm: Open-Source Self-Navigated Multi-Echo GRE for R2* and QSM Mapping

Vendor-neutral 3D multi-echo gradient echo (ME-GRE) sequence written in [Pulseq](https://pulseq.github.io), with MATLAB code to read, B0-correct and reconstruct the data for R2* and QSM mapping.

The sequence (`func/utils_pulseq/writeGradientEcho3D_ME_FreeSpacing_CAIPI.m`) provides:

- **CAIPI undersampling shifted across echoes**: each echo samples a different ky-kz pattern, so the echoes carry complementary k-space information for joint or model-based reconstruction.
- **B0 navigator**: a 1D readout during the spoiler (crusher) gradient of every TR records respiration-induced B0 fluctuations.
- **Motion navigators** (optional, `useMCnavigator`): a low-resolution (default 8 mm, `res_nav`), fully sampled k-space center spread over the acquisition, for motion estimation and coil sensitivity calibration.
- RF spoiling with a 50° phase increment, 2x readout oversampling, optional water excitation (`WaterExcite`), optional minimum readout bandwidth (`MinBwPerPixel`), and free echo spacing.

In the ISMRM 2023 abstract, the 1 mm isotropic whole-brain protocol (3T Prisma, 32-channel coil, 6 echoes, TE1/ΔTE/TE6 = 3/5/28 ms) was acquired in 4 min at R = 9 and 3 min at R = 12.

## Repository contents

| Path | Content |
|---|---|
| `script_Macro_GRE_MEsequencePulseq.m` | Generates `.seq` files for example protocols (3 mm R = 1 and R = 9, with motion navigators or variable echo spacing; 1 mm R = 9 with motion navigators) and plots the motion navigator k-space masks |
| `script_read_megre_pulseq_example_1mm.m` | Reconstruction of the 1 mm R = 9 data: reads the Siemens raw data, sorts k-space, removes oversampling, applies the navigator B0 correction, estimates ESPIRiT coil maps from the fully sampled center, and runs a slice-by-slice Tikhonov-regularized SENSE reconstruction (LSQR) followed by ARLO R2* fitting |
| `func/utils_pulseq/` | Sequence writers (`writeGradientEcho3D_ME_FreeSpacing_CAIPI.m`, earlier `writeGradientEcho3D_ME_CAIPI.m`), k-space ordering, and sequence display |
| `func/utils_jose/` | Raw data ordering and B0 correction (`ReadOrderRemoveOsB0corrSeparateNav.m`), coil maps (`CreateReceiveCoilsFromDataToData.m`), SENSE reconstruction (`ReconstructTick*.m`), ARLO R2* fitting, and display tools |
| `func/utils/` | Twix reader, ESPIRiT and SENSE helpers, FFT utilities, and `script_recon.m` (an earlier SENSE reconstruction script for the 3 mm and 1 mm data) |
| `pulseq/` | Pulseq submodule ([pulseq/pulseq](https://github.com/pulseq/pulseq) @ `8e1e3dd`) |
| `bart/` | BART submodule ([mrirecon/bart](https://github.com/mrirecon/bart) @ `5947a32`) |
| `abstracts/` | ISMRM 2023 abstract and ISMRM 2023 reconstruction workshop abstract (PDF) |

The model-based reconstruction (MOBA) of M0, R2* and B0 frequency described in the abstracts was implemented in BART. BART is included as a submodule, but the MOBA reconstruction scripts are not part of this repository.

## Usage

Clone with the submodules:

```bash
git clone --recurse-submodules https://github.com/berkinbilgic/pulseq_qsm.git
```

1. **Sequence**: in `script_Macro_GRE_MEsequencePulseq.m`, set `CurrentDir` to the repository folder and point the Pulseq path to `pulseq/matlab` (the script expects `pulseq-master\matlab`). Each protocol writes a `.seq` file and a `.mat` file with the effective sequence parameters (`SeqParamsEff`) to `compiled_seq_files/For2022Abstract/`.
2. **Reconstruction**: in `script_read_megre_pulseq_example_1mm.m`, set the data paths. The reconstruction loads the sequence parameters as variable `SeqParamsUpdated` from `SeqParamsUpdated_<name>.mat`, so save the second output of the sequence writer under that name. Set `AtDondersHPC = 0` to run locally with a parallel pool.

Raw data are not included in the repository.

## Requirements

- MATLAB with the Parallel Computing Toolbox. The 1 mm reconstruction needs about 60 GB of RAM.
- [SEPIA](https://github.com/kschan0214/sepia) for NIfTI input/output (`sepia_addpath`).
- Optional: FieldTrip `qsub` for job submission on the Donders HPC cluster.

## Third-party code

`func/utils/` includes `mapVBVD` and `twix_map_obj` by Philipp Ehses for reading Siemens raw data, and ESPIRiT calibration functions (`dat2Kernel`, `kernelEig`) from Michael Lustig's ESPIRiT code. `func/utils_jose/ARLO.m` is by Kwok-Shing Chan (Pei et al., MRM 2015).

## Reference

X Wang, B Bilgic, D Gallichan, KS Chan, JP Marques. An Open-Source Self-navigated Multi-Echo Gradient Echo Acquisition for R2* and QSM mapping using Pulseq and Model-Based Reconstruction. *ISMRM*, 2023, abstract 0402. [doi:10.58530/2023/0402](https://doi.org/10.58530/2023/0402)

Follow-up work with wave encoding (code not in this repository): X Wang, J Cho, Y Jun, B Bilgic, JP Marques. Rapid T2* and Susceptibility Mapping using Poisson Wave Encoding and Model-Based Reconstruction. *ISMRM*, 2024, abstract 3831. [doi:10.58530/2024/3831](https://doi.org/10.58530/2024/3831)

Contact: berkin AT nmr.mgh.harvard.edu
