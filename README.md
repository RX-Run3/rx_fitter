[TOC]

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

```bash
rx_fitter_data/model_tester/v2/reso_ee.yaml
```

thus, one can run this command for each model with cluster jobs.

# Resonant mode fits

## No mass constraint

For these fits run:

```bash
rx_reso_ee -m B_M_brem_track_2 -v no_dtf -b 0
```

where:

`-m` Signals which version of the mass to use    
`-v` Signals the configuration file   
`-b` Signals the brem category   

## Mass resolutions and scales

These can be obtained by running:

```bash
reso_scale
```

which will pick the latest versions of the data and MC fits

# Combinatorial

## Model validation

In order to validate the model used to fit the combinatorial, these fits are done to the SS data
with the full selection and multiple MVA working points. To do this run:

```bash
validate_cmb -q low     -c validation -s "DATA*" -t Hlt2RD_BuToKpEE_SameSign_MVA -m HypExp
validate_cmb -q central -c validation -s "DATA*" -t Hlt2RD_BuToKpEE_SameSign_MVA -m ModExp
validate_cmb -q high    -c validation -s "DATA*" -t Hlt2RD_BuToKpEE_SameSign_MVA -m SUJohnson
```

where the configuration is specified through both the arguments and the config `validation.yaml`. The latter
is specified with the `-c` flag and is part of the project itself.

# Partially reconstructed 

## PDFs
The partially reconstructed PDFs can be retrieved with:

```python
from rx_fitter import components as cmp

cmp_prc = cmp.get_kde(obs=obs, sample=sample, nbrem=nbrem, cfg=cfg)
pdf     = cmp_prc.pdf
```

`nbrem:` Brem category, in `0,1,2,None`, the last value, will put the three categories together.
`sample:` Name of the MC sample, e.g. `Bu_Kstee_Kpi0_eq_btosllball05_DPC`.
`cfg:` Dictionary with configuration.


## Scale factors

These scale factors $K^x$ are used to reparametrize the background yields as $N_{PRec}^x = K^x N_{Signal}$.
The value and error of these factors are obtained with:

```python
from rx_fitter.prec_scales    import PrecScales

obj      = PrecScales(proc=process, q2bin=q2bin)
val, err = obj.get_scale(signal=signal)
```
