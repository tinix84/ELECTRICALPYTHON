---

title: README

---

# ELECTRICALPYTHON
*Electrical-Engineering-for-Python*

Python Libraries with functions and constants related to electrical engineering.

The functions and constants that make up these modules represent a library of material compiled with the intent of being used primarily
for research, development, education, and exploration in the realm of electrical engineering.

### Special thanks to:
- Steven Weeks | Student - University of Idaho
- Dr. Dennis Sullivan | Proffessor - University of Idaho
- StackOverflow user gg349
- Paul Ortman | Power Quality Engineer - Idaho Power | Instructor - University of Idaho

In this repository, there are:
- EESIGNAL (Used primarily for signals and systems or signal/power quality)
- EEPOWER (Used primarily for electrical power calculations)
  - CAPACITOR (Special functions related to capacitors in electrical engineering)
  - PERUNIT (Special functions related to Per-Unit calculations)

### Dependencies:
- NUMPY
- MATPLOTLIB
- SCIPY


## INSTALLATION:
 1. Install required dependencies (NUMPY, SCIPY, and MATPLOTLIB)
    - `pip install numpy`
    - `pip install scipy`
    - `pip install matplotlib`
  
 2. Install *electricalpython*
    - `pip install electricalpython`
  
 3. Check installation success in Python environment

   ```
   import eepower
   eepower.ver
   ```
   
   ```
   import eesignal
   eesignal.ver
   ```
