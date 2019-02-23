
---

[Home](https://engineerjoe440.github.io/ELECTRICALPYTHON/index)
 | 
[Example Usage](https://engineerjoe440.github.io/ELECTRICALPYTHON/example)
 | 
[Function List](https://engineerjoe440.github.io/ELECTRICALPYTHON/functionlist)

---

# ELECTRICALPYTHON
*The Python Package for all Electrical Engineers!*

Electrical Engineers use a wide variety of functions, formulas, and constants. This package contains many of these devices, organized in a logical flow in the goal of expiditing calculations for electrical engineers.

### Contained Packages:
This package (*electricalpython*) includes both *eepower* and *eesignal* which each contain several sub-packages that are discretely broken up in a logical fashion. The package-breakdown is described below.

1. EEPOWER
    - CAPACITOR
    - PERUNIT
    - SYSTEMSOLUTION
2. EESIGNAL
    - BODE
    - FILTER
    - FILTERSIM

### Dependencies:
The following packages are required dependencies for this package. They are required for most, if not all, of the functions contained herein.

- NUMPY [https://pypi.org/project/numpy/](https://pypi.org/project/numpy/)
- SCIPY [https://pypi.org/project/scipy/](https://pypi.org/project/scipy/)
- SYMPY [https://pypi.org/project/sympy/](https://pypi.org/project/sympy/)
- MATPLOTLIB [https://pypi.org/project/matplotlib/](https://pypi.org/project/matplotlib/)

### Installation:
#### Install Dependencies:
The required dependencies can be installed using the following pip-install commands.


```python
pip install numpy
pip install scipy
pip install sympy
pip install matplotlib
```

#### Install electricalpython:
Following successful installation of the afore-mentioned package-dependencies, installation of *electricalpython* can be accomplished by using the following command line


```python
pip install electricalpython
```

#### Validate Installation:


```python
import eepower as eep # Import the EEPOWER package
eep.ver # Display the version information for EEPOWER
```




    '1.1.2'




```python
import eesignal as ees # Import the EESIGNAL package
ees.ver # Display the version information for EESIGNAL
```




    '2.0.1'



### Credits and Thanks:
Contributions regarding function usage, function defenition and construction were made by several generous individuals. Special thanks to:

- Stephen Weeks (B.S.E.E. University of Idaho)
- Dr. Dennis Sullivan (Ph.D Electrical Engineering \| Professor University of Idaho)
- Paul Ortmann (B.S.E.E. \| Professor University of Idaho)
