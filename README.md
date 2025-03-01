# Description

This project is meant to do fits for the rare mode. It should be replaced eventually by the RX code used in the 
Run1 and Run2 analyses. Thus, these fits are preliminary and meant as a quick crosscheck.

# Usage

To run the fits to the muon mass distribution do:

```bash
rx_fit -q central
```

For e.g. the central $q^2$ bin.

# Outputs

Once good enough fits are obtained, the parameters should be saved in versioned
directories in `rx_fitter_data/rare_fit`. Thus, other projects can pick up these
fits.

# Resonant mode fits



## Mass resolutions and scales

These can be obtained by running:

```python
reso_scale
```

which will pick the latest versions of the data and MC fits
