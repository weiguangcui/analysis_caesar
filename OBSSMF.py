import numpy as np
#every value here is log10'd

#    Input: z:redshift
#           take:all: SMF for all galaxies
#                 sf: SMF for star forming galaxies
#                  q: quenched galaxies
#    use:
#    import OBSSMF as obs
#    z = 0.5
#    smf = obs.SMF(z)
#    if smf.cond: errorbar(smf.x, smf.y, yerr=smf.yerr, label=smf.name)


#Baldry-12
BalSMF = {
      'name':'Baldry+12',
      'x':np.array([6.25, 6.75, 7.10,  7.30,  7.50,  7.70 , 7.90,  8.10,  8.30,  8.50,  8.70,  8.90,  9.10,  9.30,
                    9.50, 9.70, 9.90, 10.10, 10.30, 10.50, 10.70, 10.90, 11.10, 11.30, 11.50, 11.70, 11.90]
                  ),
      'all':np.array([
                      [ [-1.51, -1.74, -1.75, -1.37, -1.50, -1.46, -1.56, -1.55, -1.63, -1.72, -1.74, -1.84, -1.99, -2.02,
                         -2.13, -2.21, -2.24, -2.26, -2.26, -2.29, -2.45, -2.62, -2.90, -3.47, -4.38, -4.68, -4.38],
                        [ 0.23,  0.14,  0.12,  0.08,  0.11,  0.09,  0.06,  0.04,  0.05,  0.03,  0.06,  0.05,  0.02,  0.02,
                          0.02,  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.04,  0.05,  0.10,  0.23,  0.30,  0.23],
                        [ 0.23,  0.14,  0.12,  0.08,  0.11,  0.09,  0.06,  0.04,  0.05,  0.03,  0.06,  0.05,  0.02,  0.02,
                          0.02,  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.04,  0.05,  0.10,  0.23,  0.30,  0.23] ]
                     ]
                    )
       }

#Bernardi+17
BernSMF = {
      'name':'Bernardi+17',
      'x':np.arange(9.05,12.16,0.1),
      'all': np.array([
                      [ [-2.051, -2.075, -2.092, -2.089, -2.118, -2.159, -2.183, -2.208, -2.213, -2.231, -2.252, -2.240, -2.252, -2.250, -2.274, -2.314, -2.357, -2.421, -2.504, -2.611, -2.732, -2.885, -3.055, -3.252, -3.472, -3.715, -3.950, -4.254, -4.554, -4.881, -5.321, -5.742],
                        [0.022, 0.020, 0.018, 0.016, 0.015, 0.014, 0.012, 0.011, 0.010, 0.010, 0.009, 0.008, 0.007, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.007, 0.008, 0.009, 0.011, 0.014, 0.018, 0.024, 0.037, 0.055],
                        [0.022, 0.020, 0.018, 0.016, 0.015, 0.014, 0.012, 0.011, 0.010, 0.010, 0.009, 0.008, 0.007, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.007, 0.008, 0.009, 0.011, 0.014, 0.018, 0.024, 0.037, 0.055] ]
                      ]
                     )
       }

#Wright+17
WrightSMF = {
      'name':'Wright+17',
      'x':np.arange(12.0,5.24,-0.15),
      'all': np.array( [
                     [ [1.929e-05, 2.448e-05, 2.527e-04, 1.098e-04, 3.412e-04, 7.370e-04, 1.179e-03, 1.868e-03, 2.599e-03, 3.740e-03, 4.361e-03, 5.054e-03, 4.959e-03, 5.287e-03, 5.547e-03, 6.146e-03, 6.807e-03, 7.930e-03, 8.685e-03, 9.734e-03, 1.116e-02, 1.288e-02, 1.352e-02, 1.559e-02, 1.973e-02, 2.272e-02, 3.129e-02, 3.522e-02, 3.801e-02, 4.161e-02, 7.313e-02, 1.453e-01, 1.654e-01, 1.788e-01, 1.375e-01, 1.530e-01, 2.524e-01, 3.203e-01, 9.584e-02, 2.228e-01, 8.100e-01, 3.815e-01, 3.239e-02, 3.742e-02, 1.581e-02, 6.482e-03],
                     [2.623e-05, 1.783e-05, 1.778e-04, 3.004e-05, 6.006e-05, 1.122e-04, 1.220e-04, 1.594e-04, 1.688e-04, 2.112e-04, 2.050e-04, 2.468e-04, 2.334e-04, 2.191e-04, 2.183e-04, 2.460e-04, 2.752e-04, 3.172e-04, 3.244e-04, 3.789e-04, 4.573e-04, 5.924e-04, 7.481e-04, 8.825e-04, 1.278e-03, 1.896e-03, 3.190e-03, 4.228e-03, 7.141e-03, 9.480e-03, 2.022e-02, 4.561e-02, 4.349e-02, 4.503e-02, 3.401e-02, 4.164e-02, 1.498e-01, 2.135e-01, 4.514e-02, 1.855e-01, 7.009e-01, 3.857e-01, 3.304e-02, 3.572e-02, 2.140e-02, 6.445e-03],
                     [2.912e-05, 2.178e-05, 1.675e-04, 3.749e-05, 1.026e-04, 1.504e-04, 1.358e-04, 1.733e-04, 1.819e-04, 2.216e-04, 2.229e-04, 3.705e-04, 2.899e-04, 2.327e-04, 2.336e-04, 2.698e-04, 3.883e-04, 6.464e-04, 4.938e-04, 4.391e-04, 5.017e-04, 5.940e-04, 1.213e-03, 1.036e-03, 1.592e-03, 2.137e-03, 8.558e-03, 9.164e-03, 8.313e-03, 1.331e-02, 2.708e-02, 6.180e-02, 5.548e-02, 6.684e-02, 6.365e-02, 1.601e-01, 1.732e+00, 2.681e+00, 1.986e+00, 1.549e+00, 1.022e+00, 3.700e-01, 3.401e-02, 4.519e-02, 5.194e-02, 6.068e-03] ]
                      ]
                     )
      }

#Tomczak 2014
z_lim = np.asarray([0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
TomSMF = {
      'name':'Tomczak+14',
      'x':np.arange(8,11.75,0.25),
      'all':   np.array([
                         [ [-1.37,   -1.53,   -1.71,   -1.86,   -2.03,   -2.01,   -2.10, -2.17, -2.24, -2.31, -2.41, -2.53, -2.91, -3.46, np.nan],
                           [ 0.06,    0.06,    0.07,    0.07,    0.08,    0.07,    0.07,  0.08,  0.08,  0.08,  0.08,  0.09,  0.11,  0.14, np.nan],
                           [ 0.07,    0.07,    0.08,    0.08,    0.09,    0.08,    0.09,  0.10,  0.10,  0.09,  0.10,  0.11,  0.15,  0.18, np.nan] ],
                         [ [np.nan,  -1.53,   -1.60,   -1.76,   -1.86,   -2.00,   -2.12, -2.21, -2.25, -2.35, -2.45, -2.55, -2.82, -3.32, np.nan],
                           [np.nan,   0.06,    0.05,    0.06,    0.06,    0.06,    0.07,  0.06,  0.06,  0.07,  0.07,  0.08,  0.09,  0.10, np.nan],
                           [np.nan,   0.07,    0.06,    0.06,    0.07,    0.07,    0.08,  0.07,  0.08,  0.08,  0.09,  0.09,  0.11,  0.13, np.nan] ],
                         [ [np.nan, np.nan,   -1.70,   -1.86,   -2.01,   -2.10,   -2.23, -2.39, -2.45, -2.45, -2.52, -2.59, -2.93, -3.47, np.nan],
                           [np.nan, np.nan,    0.05,    0.05,    0.06,    0.06,    0.06,  0.07,  0.07,  0.07,  0.08,  0.08,  0.10,  0.11, np.nan],
                           [np.nan, np.nan,    0.06,    0.06,    0.06,    0.07,    0.07,  0.08,  0.09,  0.09,  0.09,  0.10,  0.13,  0.15, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,   -1.99,   -2.14,   -2.24,   -2.29, -2.48, -2.59, -2.73, -2.64, -2.72, -3.01, -3.62, np.nan],
                           [np.nan, np.nan,  np.nan,    0.06,    0.06,    0.06,    0.06,  0.07,  0.08,  0.08,  0.07,  0.08,  0.10,  0.11, np.nan],
                           [np.nan, np.nan,  np.nan,    0.06,    0.07,    0.07,    0.07,  0.08,  0.09,  0.10,  0.09,  0.10,  0.12,  0.15, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,   -2.02,   -2.14,   -2.28,   -2.46, -2.53, -2.61, -2.68, -2.71, -2.84, -3.12, -3.65,  -4.99],
                           [np.nan, np.nan,  np.nan,    0.06,    0.06,    0.06,    0.07,  0.07,  0.08,  0.08,  0.08,  0.08,  0.10,  0.12,   0.30],
                           [np.nan, np.nan,  np.nan,    0.07,    0.07,    0.07,    0.08,  0.08,  0.09,  0.09,  0.09,  0.10,  0.13,  0.16,   0.41] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,   -2.20,   -2.31,   -2.41, -2.54, -2.67, -2.76, -2.87, -3.03, -3.13, -3.56,  -4.27],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.05,    0.05,    0.05,  0.06,  0.06,  0.06,  0.07,  0.08,  0.08,  0.10,   0.12],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.06,    0.06,    0.06,  0.06,  0.07,  0.07,  0.08,  0.09,  0.10,  0.13,   0.15] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,   -2.53,   -2.50, -2.63, -2.74, -2.91, -3.07, -3.35, -3.54, -3.89,  -4.41],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.06,    0.06,  0.06,  0.07,  0.08,  0.09,  0.10,  0.12,  0.12,   0.14],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.07,    0.07,  0.07,  0.08,  0.09,  0.10,  0.13,  0.16,  0.17,   0.19] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,   -2.65, -2.78, -3.02, -3.21, -3.35, -3.74, -4.00, -4.14,  -4.73],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.06,  0.07,  0.08,  0.09,  0.10,  0.13,  0.18,  0.17,   0.31],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.07,  0.08,  0.09,  0.10,  0.13,  0.17,  0.25,  0.28,   2.00] ]
                         ]
                         ),

      'sf':  np.array([
                         [ [-1.42,   -1.59,   -1.76,   -1.91,   -2.08,   -2.06,   -2.17, -2.25, -2.36, -2.50, -2.63, -2.91, -3.43, -4.39, np.nan],
                           [ 0.06,    0.06,    0.07,    0.07,    0.08,    0.07,    0.07,  0.08,  0.08,  0.08,  0.09,  0.10,  0.13,  0.30, np.nan],
                           [ 0.07,    0.07,    0.08,    0.08,    0.09,    0.08,    0.08,  0.10,  0.10,  0.09,  0.11,  0.12,  0.18,  0.41, np.nan] ],
                         [ [np.nan,  -1.60,   -1.67,   -1.83,   -1.92,   -2.09,   -2.19, -2.28, -2.39, -2.55, -2.76, -3.00, -3.46, -4.30, np.nan],
                           [np.nan,   0.06,    0.05,    0.06,    0.06,    0.06,    0.07,  0.06,  0.07,  0.07,  0.08,  0.08,  0.10,  0.20, np.nan],
                           [np.nan,   0.07,    0.06,    0.06,    0.07,    0.07,    0.08,  0.07,  0.08,  0.08,  0.09,  0.10,  0.14,  0.25, np.nan] ],
                         [ [np.nan, np.nan,   -1.72,   -1.88,   -2.04,   -2.14,   -2.27, -2.47, -2.55, -2.60, -2.77, -2.91, -3.37, -4.17, np.nan],
                           [np.nan, np.nan,    0.05,    0.05,    0.06,    0.06,    0.06,  0.07,  0.08,  0.07,  0.08,  0.09,  0.10,  0.16, np.nan],
                           [np.nan, np.nan,    0.06,    0.06,    0.06,    0.07,    0.07,  0.08,  0.09,  0.09,  0.09,  0.11,  0.13,  0.20, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,   -2.00,   -2.16,   -2.26,   -2.32, -2.52, -2.68, -2.88, -2.81, -2.99, -3.29, -4.21, np.nan],
                           [np.nan, np.nan,  np.nan,    0.06,    0.06,    0.06,    0.06,  0.07,  0.08,  0.09,  0.07,  0.08,  0.10,  0.15, np.nan],
                           [np.nan, np.nan,  np.nan,    0.06,    0.07,    0.07,    0.07,  0.08,  0.07,  0.10,  0.09,  0.10,  0.13,  0.20, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,   -2.03,   -2.15,   -2.29,   -2.48, -2.55, -2.68, -2.75, -2.87, -3.07, -3.39, -3.95,  -5.17],
                           [np.nan, np.nan,  np.nan,    0.06,    0.06,    0.06,    0.07,  0.07,  0.08,  0.08,  0.08,  0.08,  0.10,  0.13,   0.37],
                           [np.nan, np.nan,  np.nan,    0.07,    0.07,    0.07,    0.08,  0.08,  0.09,  0.10,  0.09,  0.10,  0.13,  0.17,   0.52] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,   -2.20,   -2.32,   -2.42, -2.56, -2.73, -2.89, -3.07, -3.26, -3.35, -3.85,  -4.78],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.05,    0.05,    0.05,  0.06,  0.06,  0.07,  0.07,  0.09,  0.09,  0.10,   0.17],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.06,    0.06,    0.06,  0.06,  0.07,  0.07,  0.09,  0.10,  0.11,  0.13,   0.21] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,   -2.53,   -2.51, -2.67, -2.78, -3.00, -3.26, -3.54, -3.69, -4.00,  -4.59],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.06,    0.06,  0.06,  0.07,  0.08,  0.09,  0.11,  0.13,  0.13,   0.15],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.07,    0.07,  0.07,  0.08,  0.09,  0.11,  0.14,  0.17,  0.17,   0.21] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,   -2.66, -2.79, -3.06, -3.32, -3.59, -3.97, -4.16, -4.32,  -4.94],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.06,  0.07,  0.08,  0.09,  0.11,  0.16,  0.20,  0.18,   0.32],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.07,  0.08,  0.09,  0.11,  0.14,  0.20,  0.28,  0.29,   2.00] ]
                         ]
                         ),
      'q':   np.array([
                         [ [np.nan,  -2.41,   -2.62,   -2.82,   -2.96,   -2.96,   -2.98, -2.91, -2.86, -2.78, -2.80, -2.76, -3.07, -3.52, np.nan],
                           [np.nan,   0.08,    0.10,    0.12,    0.14,    0.08,    0.09,  0.09,  0.09,  0.08,  0.09,  0.09,  0.12,  0.14, np.nan],
                           [np.nan,   0.10,    0.11,    0.14,    0.16,    0.10,    0.10,  0.11,  0.11,  0.10,  0.11,  0.12,  0.16,  0.19, np.nan] ],
                         [ [np.nan, np.nan,   -2.42,   -2.58,   -2.77,   -2.75,   -2.94, -2.99, -2.83, -2.78, -2.75, -2.75, -2.93, -3.37, np.nan],
                           [np.nan, np.nan,    0.07,    0.07,    0.09,    0.09,    0.10,  0.07,  0.07,  0.07,  0.08,  0.08,  0.09,  0.11, np.nan],
                           [np.nan, np.nan,    0.08,    0.08,    0.10,    0.10,    0.11,  0.08,  0.08,  0.09,  0.09,  0.10,  0.11,  0.14, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,   -3.19,   -3.17,   -3.33, -3.16, -3.16, -2.97, -2.89, -2.87, -3.12, -3.57, np.nan],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.11,    0.10,    0.12,  0.11,  0.11,  0.08,  0.08,  0.09,  0.10,  0.12, np.nan],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.12,    0.12,    0.14,  0.12,  0.12,  0.09,  0.10,  0.11,  0.13,  0.15, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,   -3.46,   -3.65,   -3.46, -3.57, -3.37, -3.26, -3.11, -3.05, -3.33, -3.75, np.nan],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.12,    0.15,    0.13,  0.14,  0.12,  0.11,  0.08,  0.08,  0.10,  0.12, np.nan],
                           [np.nan, np.nan,  np.nan,  np.nan,    0.14,    0.17,    0.14,  0.16,  0.14,  0.13,  0.09,  0.09,  0.13,  0.16, np.nan] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,   -3.97,   -3.79, -3.75, -3.45, -3.52, -3.24, -3.23, -3.46, -3.95,  -5.47],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.21,    0.17,  0.16,  0.12,  0.13,  0.08,  0.09,  0.10,  0.13,   0.52],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,    0.24,    0.19,  0.18,  0.14,  0.15,  0.10,  0.11,  0.13,  0.17,   0.90] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,   -4.41, -3.95, -3.55, -3.35, -3.30, -3.40, -3.54, -3.87,  -4.44],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.17,  0.14,  0.09,  0.08,  0.08,  0.09,  0.09,  0.10,   0.13],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,    0.19,  0.15,  0.11,  0.09,  0.09,  0.11,  0.11,  0.13,   0.16] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan, -3.72, -3.76, -3.64, -3.53, -3.82, -4.08, -4.54,  -4.89],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  0.11,  0.11,  0.11,  0.10,  0.13,  0.17,  0.15,   0.19],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  0.12,  0.13,  0.12,  0.12,  0.16,  0.22,  0.21,   0.26] ],
                         [ [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan, -4.16, -4.08, -3.89, -3.74, -4.12, -4.51, -4.61,  -5.14],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  0.17,  0.16,  0.13,  0.12,  0.18,  0.37,  0.19,   0.34],
                           [np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  0.20,  0.18,  0.15,  0.15,  0.22,  0.38,  0.32,   2.00] ]
                         ]
                         )
      }

# Song+15
z_song = [4.0, 5.0, 6.0]
SongSMF = {
      'name':'Song+15',
      'x':np.arange(7.75,11.26,0.5),
      'all':   np.array([
                         [ [-1.77,-2.00,-2.22,-2.52,-2.91,-3.41,-4.11,-5.00],
                           [0.15,0.13,0.09,0.09,0.12,0.13,0.22,0.24],
                           [0.24,0.10,0.09,0.09,0.05,0.08,0.21,0.97] ],
                         [ [-1.72,-2.01,-2.33,-2.68,-3.12,-3.63,-4.40,-5.96],
                           [0.20,0.16,0.15,0.07,0.09,0.13,0.15,0.0],
                           [0.20,0.16,0.10,0.14,0.11,0.11,0.35,0.0] ],
                         [ [-1.81,-2.26,-2.65,-3.14,-3.69,-4.55,-5.96,-7.50],
                           [0.23,0.21,0.15,0.12,0.12,0.19,0.52,0.0],
                           [0.28,0.16,0.15,0.11,0.13,0.24,0.32,0.0] ]
                        ])
        }

def sel_obsmf(z,name=None):
    print('OBSSMF:',z,z_lim[0])
    if z<z_lim[0]: 
        if name == 'Baldry+12': return BernSMF, 0, 0
        elif name == 'Bernardi+17': return BalSMF, 0, 0
        elif name == 'Wright+17': return WrightSMF, 0, 0
        else: return WrightSMF, 0, 0
    if z<z_lim[-1]+0.5*(z_song[0]-z_lim[-1]): 
        i = np.argmin(np.abs(z_lim-z))
        if z < z_lim[0] + 0.5*(z_lim[1]+z_lim[0]):
            return TomSMF, i-1, i-1
        else:
            return TomSMF, i-1, np.amin([i,len(z_lim)-2])
# must be high-z... use Song+15
    i = np.argmin(np.abs(z-z_song))
    return SongSMF, i, i


class SMF(object):

    '''
    Input: z:redshift
           take:all: SMF for all galaxies
                 sf: SMF for star forming galaxies
                  q: quenched galaxies
    use:
    import OBSSMF as obs
    z = 0.5
    smf = obs.SMF(z)
    if smf.cond: errorbar(smf.x, smf.y, yerr=smf.yerr, label=smf.name)
           
    '''
    def __init__(self, z, take='all', name=None):


        cond = 1
        if z<z_lim[0] or z>z_lim[-1]:
            if take!='all': cond = 0

        self.cond = cond

        if not cond: print('No Observational data for %s galaxies at z=%.2f'%(take,z))

        if cond:
            self.z     = z
            self.xerr  = 0
            obsmf      = sel_obsmf(z,name)
            SMF_   = obsmf[0]
            index_ = obsmf[1:]
            self.name = SMF_['name']
            self.x    = SMF_['x']
            y0        = SMF_[take][index_[0]]
            y1        = SMF_[take][index_[1]]
            if z <= 0.12 and self.name == 'Wright+17':  # Wright is already in non-log
                self.y    = np.nanmean(np.array([y0[0], y1[0]]), axis=0)
            else:
                self.y    = np.nanmean(np.array([10**y0[0], 10**y1[0]]), axis=0)
            self.y    = np.log10(self.y)
            if z <= 0.12 and self.name == 'Wright+17':  # Wright is already in non-log
                err_up    = np.log10(y0[0]+y0[2])-np.log10(y0[0])
                err_down  = np.log10(y0[0])-np.log10(y0[0]-y0[1])
            else: 
                err_up    = np.nanmean(np.array([y0[1], y1[1]]), axis=0)
                err_down  = np.nanmean(np.array([y0[2], y1[2]]), axis=0)
            self.yerr  = [err_down,err_up]


