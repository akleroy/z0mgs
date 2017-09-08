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

To do's:
    Try if I can reproduce the color correction factors

Note:
    Need the folder 'RSRF' in the same directory with the code
"""
import os
import numpy as np
from scipy.integrate import cumtrapz
# Some definitions
repr_wl = {'SPIRE_500': 500, 'SPIRE_350': 350, 'SPIRE_250': 250,
           'PACS_160': 160, 'PACS_100': 100, 'PACS_70': 70, 'MIPS_160': 160,
           'MIPS_70': 70, 'MIPS_24': 24}


def z0mg_RSRF(wl, s, bands=['SPIRE_350', 'SPIRE_500']):
    """
    <Inputs>
        wl (list or 1-d numpy.ndarray; shape (n)):
            Wavelengths of the input spectrum in micrometer.
        s (list or 1-d numpy.ndarray; shape (n)):
            Input spectrum.
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
        if bands[i] not in repr_wl.keys():
            raise ValueError('The instrument \"' + bands[i] +
                             "\" is currently not supported")
    #
    # Generating dir
    rsrf_dir, _ = os.path.split(__file__)
    rsrf_dir += '/RSRF/'
    #
    # Results array
    results = []
    #
    for band in bands:
        #
        # Reading RSRF data
        fn = rsrf_dir + band + '.csv'
        data = np.genfromtxt(fn, delimiter=',', skip_header=1)
        #
        # Interpolate the RSRF to match input wavelength
        rsrf = np.interp(wl, data[:, 0], data[:, 1])
        results.append(cumtrapz(s * rsrf, wl)[-1] /
                       cumtrapz(rsrf / wl / repr_wl[band], wl)[-1])
    return np.array(results)


def idc_test(SigmaD, T, beta, bands):
    from astropy.constants import c, h, k_B
    import astropy.units as u

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

    wl = np.linspace(1.0, 600.0, 1000)
    s = model(wl, SigmaD, T, beta).tolist()
    results = z0mg_RSRF(wl.tolist(), s, bands)
    print('Result:', results)
    print('Model :', model(np.array([repr_wl[b] for b in bands]), SigmaD, T,
                           beta))


if __name__ == "__main__":
    print('\nExcuted as \"__main__\". Entering test mode.')
    bands = ['MIPS_24', 'MIPS_70', 'PACS_70', 'PACS_100', 'MIPS_160',
             'PACS_160', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500']
    SigmaD, T, beta = 1E-2, 20, 2
    print('Generating test spectrum in', bands)
    print('SigmaD, T, beta =', str(SigmaD) + ',', str(T) + ',', str(beta))
    idc_test(SigmaD, T, beta, bands)
