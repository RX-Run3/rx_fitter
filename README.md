# Description

This project is meant to do fits for the rare mode. It should be replaced eventually by the RX code used in the 
Run1 and Run2 analyses. Thus, these fits are preliminary and meant as a quick crosscheck.

# Outputs

Once good enough fits are obtained, the parameters should be saved in versioned
directories in `rx_fitter_data/rare_fit`. Thus, other projects can pick up these
fits.

# Model choice

In order to find the best fitting model one can run:

```bash
model_tester -m mod_001 -b  0 -o B_M -v v2
```

which is going to fit the `B_M` observable with the model `mod_001` as defined in the config file

```
rx_fitter_data/model_tester/v2/reso_ee.yaml
```

thus, one can run this command for each model with cluster jobs.

# Resonant mode fits

## Mass resolutions and scales

These can be obtained by running:

```python
reso_scale
```

which will pick the latest versions of the data and MC fits

