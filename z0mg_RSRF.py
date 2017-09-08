"""
Author:
    I-Da Chiang (idchiang@ucsd.edu)

Name:
    z0mg_RSRF.py

Purpose:
    Giving an IR spectrum, convolve it with selected instruments' RSRF to
    predict the observational outcome

Current available instruments:
    PACS_70, PACS_100, PACS_160, SPIRE_250, SPIRE_350, SPIRE_500
    MIPS_24, MIPS_100, MIPS_160

Current problems:
    No extrapolation of response function
        Response function range: MIPS_24  (18.01, 32.21)
                                 MIPS_70  (49.96, 111.02)
                                 MIPS_160 (100.09, 199.92)
                                 PACS_70  (2.50, 500.0)
                                 PACS_100 (13.49, 500.0)
                                 PACS_160 (13.49, 500.0)
                                 SPIRE    (166.67, 1000)
    Unable to reproduce color correction factors in PACS_160 for T <= 7K

Note:
    Need the folder 'RSRF' in the same directory with the code
"""
import os
import numpy as np
import astropy.units as u
from scipy.integrate import cumtrapz
from astropy.constants import h, k_B, c
# Reference wavelength for each instrument
ref_wl = {'SPIRE_500': 500, 'SPIRE_350': 350, 'SPIRE_250': 250,
          'PACS_160': 160, 'PACS_100': 100, 'PACS_70': 70,
          'MIPS_160': 155.899, 'MIPS_70': 71.440, 'MIPS_24': 23.675}


# Main function
def z0mg_RSRF(wl, s, bands=['SPIRE_350', 'SPIRE_500']):
    """
    <Inputs>
        wl (list or 1-d numpy.ndarray; shape (n)):
            Wavelengths of the input spectrum in micrometer.
        s (list or 1-d numpy.ndarray; shape (n)):
            Input spectrum in Jy/sr or equivalent (The point is, /frequncy
            instead of /wavelength).
        bands (list or 1-d numpy.ndarray; shape (m)):
            Selected output bands, e.g., ['PACS_100', 'PACS_160', 'SPIRE_250']
    <Outputs>
        results (numpy.ndarray; shape (m)):
            Resulting spectrum with sequence corresponding to bands.
    """
    #
    # Making sure WL and S are having the same 1-d, non-zero length
    assert len(wl) == len(s)
    assert len(wl) > 0
    if type(wl) == np.ndarray:
        assert wl.ndim == 1
    elif type(wl) == list:
        wl = np.array(wl)
    if type(s) == np.ndarray:
        assert s.ndim == 1
    elif type(s) == list:
        s = np.array(s)
    #
    # Due to turning the RSRF integral from d\nu to d\lambda,
    # there is a \frac{-c}{\lambda^2} factor correction.
    # (-c will be canceled in the process)
    s /= wl**2
    #
    # Making sure bands in expected formats
    assert len(bands) > 0
    for i in range(len(bands)):
        for j in range(len(bands[i]))[::-1]:
            try:
                int(bands[i][j])
            except ValueError:
                num = bands[i][j + 1:]
                break
        if bands[i][:4].upper() in ['PACS', 'MIPS']:
            bands[i] = bands[i][:4].upper() + '_' + num
        elif bands[i][:5].upper() in ['SPIRE']:
            bands[i] = bands[i][:5].upper() + '_' + num
        if bands[i] not in ref_wl.keys():
            raise ValueError('The instrument \"' + bands[i] +
                             "\" is currently not supported")
    #
    # Generating dir
    rsrf_dir, _ = os.path.split(__file__)
    rsrf_dir += '/RSRF/'
    #
    # Results array
    results = []

    # Simple 10000K blackbody for MIPS reference
    def B_10000K(wl):
        temp = (h * c / (wl * u.um) / k_B / (10000 * u.K)).decompose().value
        return 1 / (np.exp(temp) - 1) / wl**3
    #
    for band in bands:
        #
        # Reading RSRF data
        fn = rsrf_dir + band + '.csv'
        data = np.genfromtxt(fn, delimiter=',', skip_header=1)
        #
        # Some of the response functions are not flat at the boundary, thus
        # extrapolation might give negative value.
        mask = (np.min(data[:, 0]) <= wl) * (wl <= np.max(data[:, 0]))
        if np.sum(mask) == 0:
            raise ValueError("There is no available wavelength points inside",
                             band + ', response function range, which is (' +
                             str(round(np.min(data[:, 0]), 3)) + ', ' +
                             str(round(np.max(data[:, 0]), 3)) + ') micron.')
        wlb = wl[mask]
        #
        # Interpolate the RSRF to match input wavelength
        rsrf = np.interp(wlb, data[:, 0], data[:, 1])
        if band[:4] == 'MIPS':
            results.append(cumtrapz(s[mask] * rsrf, wlb)[-1] /
                           cumtrapz(rsrf * B_10000K(wlb) /
                                    B_10000K(ref_wl[band]) / wlb**2, wlb)[-1])
        else:
            results.append(cumtrapz(s[mask] * rsrf, wlb)[-1] /
                           cumtrapz(rsrf / wlb / ref_wl[band], wlb)[-1])
    return np.array(results)


def idc_test(SigmaD, T, beta, bands):
    def B(T, freq):
        """Return blackbody SED of temperature T(with unit) in MJy"""
        with np.errstate(over='ignore'):
            return (2 * h * freq**3 / c**2 / (np.exp(h * freq / k_B / T) - 1)
                    ).to(u.Jy).value * 1E-6

    def model(wl, SigmaD, T, beta):
        """Return fitted SED in MJy"""
        kappa160 = 9.6 * np.pi
        const = 2.0891E-4 * kappa160
        freq = (c / wl / u.um).to(u.Hz)
        return const * (160.0 / wl)**beta * SigmaD * B(T * u.K, freq)

    wl = np.linspace(2.5, 1000, 10**5)
    s = model(wl, SigmaD, T, beta)
    m = model(np.array([ref_wl[b] for b in bands]), SigmaD, T, beta)
    return z0mg_RSRF(wl, s, bands) / m


if __name__ == "__main__":
    print('\nExcuted as \"__main__\". Entering test mode.')
    print('\nBegin Color correction factor (CCF) test.')
    SigmaD = 1E-2
    # MIPS: Table 4.16 in handbook v3 (p.98)
    # MIPS RSRF is originally in 1/wavelength
    beta = 0
    bands = ['MIPS_24', 'MIPS_70', 'MIPS_160']
    MIPS_T = [1E4, 5000, 1000, 500, 300, 200, 150, 100, 70, 50, 30, 20]
    MIPS_CCF = [[1.0, 1.0, 1.0], [0.999, 0.999, 1.0], [0.992, 0.995, 0.999],
                [0.983, 0.989, 0.997], [0.970, 0.980, 0.996],
                [0.957, 0.970, 0.993], [0.948, 0.959, 0.991],
                [0.947, 0.938, 0.986], [0.986, 0.914, 0.979],
                [1.119, 0.893, 0.971], [2.031, 0.901, 0.954],
                [7.005, 1.052, 0.944]]
    assert len(MIPS_T) == len(MIPS_CCF)
    print('\nMIPS_24, MIPS_70, MIPS_160: Diff in CCF / handbook CCF')
    for i in range(len(MIPS_T))[::-1]:
        Cal = idc_test(SigmaD, MIPS_T[i], beta, bands)
        p = (np.array(MIPS_CCF[i]) - Cal) / np.array(MIPS_CCF[i])
        print('T =', str(MIPS_T[i]) + ':', p, np.sum(np.abs(p) <= 0.05) == 3)
    # SPIRE: Table 5.8 for ext sources in handbook v3.1 (p.105)
    beta = 2  # Didn't check the beta=1.5 part. One should be enough.
    bands = ['SPIRE_250', 'SPIRE_350', 'SPIRE_500']
    SPIRE_T = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40]
    SPIRE_CCF = [[0.5479, 0.7986, 0.9023], [0.7679, 0.9361, 1.0119],
                 [0.8915, 0.9939, 1.0470], [0.9586, 1.0169, 1.0470],
                 [0.9944, 1.0246, 1.0523], [1.0130, 1.0253, 1.0464],
                 [1.0220, 1.0227, 1.0395], [1.0254, 1.0188, 1.0329],
                 [1.0153, 0.9984, 1.0079], [0.9998, 0.9846, 0.9938],
                 [0.9881, 0.9757, 0.9853], [0.9796, 0.9697, 0.9796],
                 [0.9734, 0.9655, 0.9756], [0.9688, 0.9623, 0.9727]]
    assert len(SPIRE_T) == len(SPIRE_CCF)
    print('\nSPIRE_250, SPIRE_350, SPIRE_500: Diff in CCF / handbook CCF')
    for i in range(len(SPIRE_T)):
        Cal = 1 / idc_test(SigmaD, SPIRE_T[i], beta, bands)
        p = (np.array(SPIRE_CCF[i]) - Cal) / np.array(SPIRE_CCF[i])
        print('T =', str(SPIRE_T[i]) + ':', p, np.sum(np.abs(p) <= 0.05) == 3)
    # PACS: Table 1 in Muller, T. et al 2011 (PICC-ME-TN-038)
    beta = 0
    bands = ['PACS_70', 'PACS_100', 'PACS_160']
    PACS_T = [10000, 5000, 1000, 500, 250, 100, 50, 40, 30, 20, 19, 18, 17,
              16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5]
    PACS_CCF = [[1.016, 1.034, 1.074], [1.016, 1.033, 1.074],
                [1.013, 1.031, 1.072], [1.011, 1.029, 1.068],
                [1.005, 1.023, 1.062], [0.989, 1.007, 1.042],
                [0.982, 0.985, 1.010], [0.992, 0.980, 0.995],
                [1.034, 0.982, 0.976], [1.224, 1.036, 0.963],
                [1.269, 1.051, 0.964], [1.325, 1.069, 0.967],
                [1.396, 1.093, 0.972], [1.488, 1.123, 0.979],
                [1.607, 1.162, 0.990], [1.768, 1.213, 1.005],
                [1.992, 1.282, 1.028], [2.317, 1.377, 1.061],
                [2.816, 1.512, 1.110], [3.645, 1.711, 1.184],
                [5.175, 2.024, 1.300], [8.497, 2.554, 1.491],
                [17.815, 3.552, 1.833], [58.391, 5.774, 2.528],
                [456.837, 12.259, 4.278]]
    assert len(PACS_T) == len(PACS_CCF)
    print('\nPACS_70, PACS_100, PACS_160: Diff in CCF / handbook CCF')
    for i in range(len(PACS_T))[::-1]:
        # Inversed definition
        Cal = idc_test(SigmaD, PACS_T[i], beta, bands)
        p = (np.array(PACS_CCF[i]) - Cal) / np.array(PACS_CCF[i])
        print('T =', str(PACS_T[i]) + ':', p, np.sum(np.abs(p) <= 0.05) == 3)
