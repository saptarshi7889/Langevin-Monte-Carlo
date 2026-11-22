# Changes

## Unreleased

- Group the initial and continued programs by incommensurate and commensurate models.
- Add a C++17 Make build, numerical and restart regression tests, and GCC/Clang CI.
- Replace hard-coded run parameters with a validated command-line interface.
- Add versioned checkpoints containing the unshifted field, RNG and Gaussian state,
  sampling schedule, counters, and accumulated order parameter and roughness.
- Restore and update roughness correctly; support rectangular lattices.
- Preserve the field phase for commensurate order parameters and checkpoints.
- Evaluate the generalized beta-function kernel using reflected positive arguments;
  preserve the original half-open image range, and support 0 < s < 2.
- Refuse existing output directories and require explicit seeds for distinct trajectories.
- Add MIT licensing and software citation metadata.

The old filename-based restart files are not compatible with the complete checkpoint
format. The new CLI `g` denotes the paper's action coupling; the local force coefficient
is calculated internally. See README.md for migration details and verification limits.
