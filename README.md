# PSR_paleo
summer research; paleoclimate variability in AMOC from Proxy Surrogate Reconstruction; O18 and Mg/Ca

This read me file describes the scripts and work flow for generating pseudoproxy Mg/Ca and O18 data. The files below are listed in order of work flow and should be executed as follows. (need to add something about purpose of pseudoproxy data creation)

sim_PSMbin_SP.R – annual temperature, salt, and H2O18 data forward models calcite d18O and Mg/Ca. For pseudo proxy construction, only lines 1-210 are needed, this will generate data with the form XX_YY, where XX is either piC or LM, YY is O18, MgCa, AMOC, NAM or NAO. Output with O18 or MgCa will be used to create the pseudoproxy data.

Proxybin_anomIND_SP.R – pulls out metadata for O18 and Mg/Ca needed for binning and smoothing; only lines 1-31 need to be executed

pseudoprox_PSR_SP.R – this file generates pseudoproxies from the model data created in sim_PSMbin_SP.R. Noise is created with a signal-to-noise ratio of 0.5 and added to the model data. Output data is in the form of XX_YY_PP. Resulting pseudo proxy data should be 1-to-1 with the original model data.

PPsmoothing_SP.R – marine sediment rows (sites) in each model are smoothed using a weighted moving average and guassian filter to represent bioturbation. The length of the window in each site is chosen based on the real proxy data sample resolution of the corresponding row. Non-marine sediment rows are pulled out prior to smoothing and added back after (but before binning). Output data is in the form of XX_YY_PPsmooth. Note: sample resolutions close to 1 sample/year are indicative of annual resolution and subsequently ignored as they already represent bioturbation.

PPaggregate_SP.R – this script adds non marine sediment rows back into the pseudoproxies into their respective rows. Output data is in the form of XX_YY_PPagg.

PPaveraging_SP.R – smoothed pseudoproxy data is binned to mimic the frequency of the real proxy data. Bins were calculated based on the real proxy data sample resolutions for Mg/Ca and O18 data. Note: sample resolutions close to 1 sample/year are indicative of annual resolution and subsequently ignored as they already represent bioturbation.
