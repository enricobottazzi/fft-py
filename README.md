# fft-py

Experiements on FFT

```
export DPS=41
python3 -m unittest test.py
```

`DPS` is the number of decimal places for `mpmath` library. It is used to set the precision of the floating point numbers during FFT operations. When working with 60-bit coefficients, a precision of 41 is enough to avoid any overflow.

