# AAFFT-Arduino
Implementation of the AAFFT data sampling algorithm on the Arduino platform for applications ranging from simple experimental data collection to efficient collection and transmission of IoT (Internet of Things) data. 

[Survey article detailing the algorithm](http://users.cms.caltech.edu/~jtropp/papers/GST08-Tutorial-Fast.pdf)

[Code upon which this is based](https://github.com/annacgilbert/Simple-sublinear-Fourier-sampling)

# Libraries
- [StandardCplusplus](https://github.com/maniacbug/StandardCplusplus)
- [Complex](https://github.com/RobTillaart/Arduino/tree/master/libraries/Complex)
- [fix_fft](https://github.com/TJC/arduino/tree/master/sketchbook/libraries/fix_fft) (prog_int8_t in fix_fft.cpp was changed to int8_t)

# Components checklist
- [ ] test_AAFFT_with_presampling (driver)
- [x] sample_shattering
- [x] sample_residual
- [x] eval_sig
- [ ] eval_sig_vect
- [ ] identify_frequencies
- [x] generate_tspairs
- [x] generate_signal
- [x] generate_sample_set
- [x] eval_sig
- [x] estimate_coeffs
