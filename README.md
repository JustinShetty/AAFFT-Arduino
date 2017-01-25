# AAFFT-Arduino
Implementation of the AAFFT data sampling algorithm on the Arduino platform for applications ranging from simple experimental data collection to efficient collection and transmission of IoT (Internet of Things) data. 

[Survey article detailing the algorithm](http://users.cms.caltech.edu/~jtropp/papers/GST08-Tutorial-Fast.pdf)

[Code upon which this is based](https://github.com/annacgilbert/Simple-sublinear-Fourier-sampling)

# Libraries Utilized
- [StandardCplusplus](https://github.com/maniacbug/StandardCplusplus)
- [Complex](https://github.com/RobTillaart/Arduino/tree/master/libraries/Complex)
- [FFT](http://wiki.openmusiclabs.com/wiki/ArduinoFFT)

# Component Checklist
- [ ] test_AAFFT_with_presampling (driver)
- [ ] sample_shattering
- [x] sample_residual
- [ ] identify_frequencies
- [x] generate_tspairs
- [x] generate_signal
- [ ] generate_sample_set
- [x] eval_sig
- [ ] estimate_coeffs