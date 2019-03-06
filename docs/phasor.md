
---

[Home](https://engineerjoe440.github.io/ELECTRICALPYTHON/index)
 | 
[Example Usage](https://engineerjoe440.github.io/ELECTRICALPYTHON/example)
 | 
[Function List](https://engineerjoe440.github.io/ELECTRICALPYTHON/functionlist)

---

# ELECTRICALPYTHON
*The Python Package for all Electrical Engineers!*

### Example Usage:
The following presents a simple usage of the *eepower.phasor* function.

#### Function Usage:
The *phasor* function will accept the magnitude and angle (in degrees) of a complex variable (i.e. a Voltage or Current) and return the complex value corresponding to the inputs
represented in Python's native (R + jI) scheme.

#### Required Function Arguments:
- mag: (FLOAT) The Magnitude of the Complex Variable
- ang: (FLOAT) The Angle (in degrees) of the Complex Variable

#### Optional Function Arguments:
- NONE

#### Returns:
- phasor: (COMPLEX) Standard Pythonic Complex Representation of the specified voltage or current.

#### Function Demonstration:


```python
>>> import eepower as eep # Import EEPOWER
>>> mag = 120 # |V|
>>> ang = 30 # ∠V
>>> x = eep.phasor(mag, ang)
>>> x
(103.92304845413264+59.99999999999999j)

>>> eep.cprint(x)
 120.0 ∠ 30.0°
```
