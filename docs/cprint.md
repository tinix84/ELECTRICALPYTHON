
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
The following presents a simple usage of the *eepower.cprint* function.

#### Function Usage:
The *cprint* function will print a complex variable in the standard format of MAG ∠ ANG° with formatting prepended text and appended unit text. The function
can also return a set of values corresponding to the magnitude and angle of the complex value.

#### Required Function Arguments:
- val: (COMPLEX) The Complex Value to be Printed

#### Optional Function Arguments:
- unit: (STRING) The string to be printed corresponding to the unit mark. default=""
- label: (STRING) The prepended string used as a descriptive labeling string. default=""
- printval: (BOOL) Control argument enabling/disabling printing of the string. default=True
- ret: (BOOL) Control argument allowing the evaluated value to be returned. default=False
- decimals: (INTEGER) Control argument specifying how many decimals of the complex value to be printed. May be negative to round to spaces
to the left of the decimal place (follows standard round() functionality). default=3

#### Returns:
- numarr: (NUMPY.NDARRAY) The array of values corresponding to the magnitude and angle, values are returned in the form:
            
    [[ mag, ang ],
     [ mag, ang ],
          ...    ,
     [ mag, ang ]]

#### Function Demonstration:


```python
>>> import eepower as eep # Import EEPOWER
>>> V = eep.phasor(60,-30)
>>> eep.cprint(V)
 60.0 ∠ -30.0° 

>>> eep.cprint(V,"Volts","My Output:")
My Output: 60.0 ∠ -30.0° Volts

>>> import numpy as np
>>> Vset = np.array([[eep.phasor(50,0)],[eep.phasor(50,-120)]])
>>> eep.cprint(Vset)
[[' 50.0 ∠ 0.0° ']
 [' 50.0 ∠ -120.0° ']]

>>> eep.cprint(Vset,"V")
[[' 50.0 ∠ 0.0° V']
 [' 50.0 ∠ -120.0° V']]

>>> VIset = np.array([[eep.phasor(100, 0)],
		      [eep.phasor(10, -30)]])
>>> x = eep.cprint(VIset,["V","A"],["Voltage:","Current"],ret=True)
[['Voltage: 100.0 ∠ 0.0° V']
 ['Current 10.0 ∠ -30.0° A']]
>>> x
array([[100.,   0.],
       [ 10., -30.]])
```
