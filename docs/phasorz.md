
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
The following presents a simple usage of the *eepower.phasorz* function.

#### Function Usage:
The *phasorz* function is designed to generate the phasor representation of a given inductance or capacitance at a specified frequency.

#### Required Function Arguments:
- f: (FLOAT) The specified system frequency at which the phasor impedance will be calculated

#### Optional Function Arguments:
*(Exclusive Requirement: Either C or F must be specified, not both.)*
- C: (FLOAT) The capacitance for which the impedance should be calculated (specified in Farads). default=None
- F: (FLOAT) The inductance for which the impedance should be calculated (specified in Henrys). default=None
- complex: Control argument to specify whether the returned value should be returned as a complex value. default=True

#### Returns:
- Z: (COMPLEX/FLOAT) The calculated impedance of the passive element (C or L) at the specified system frequency.

#### Function Demonstration:


```python
>>> import eepower as eep # Import EEPOWER
>>> frequency = 60 # Hz
>>> capacitance = 1e-6 # Farads
>>> eep.phasorz(frequency, C=capacitance)
(-0-2652.5823848649225j)

>>> inductance = 4e-3 # Henrys
>>> eep.phasorz(frequency, L=inductance,complex=False)
1.5079644737231006
```
