
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
The following presents a simple usage of the *eepower.powertriangle* function.

#### Function Usage:
The *powertriangle* function will accept two arguments to generate a Power-Triangle. In addition to it's primary goal of generating and plotting the triangle, the function can also return the remaining unknowns and print all values on the plot.

#### Required Function Arguments (any two of the following):
- P: (FLOAT) Real Power (Watts, absolute value - real, not complex)
- Q: (FLOAT) Reactive Power (VARs, absolute value - real, not complex)
- S: (FLOAT) Apparent Power (VAs, absolute value - real, not complex)
- PF: (FLOAT) Power Factor (expressed as a decimal value - i.e. 0.8)

#### Optional Function Arguments:
- color: (STRING) Line color of the plotted triangle, default is red.
- text: (STRING) Text to appear as the plot title, default is "Power Triangle".
- printval: (BOOL) Argument specifying whether the values should be printed on the plot, default is False.
- ret: (BOOL) Argument specifying whether the values should be returned after the function call, default is False.
- plot: (BOOL) Argument specifying whether the plot should be generated, default is True.

#### Function Demonstration:


```python
>>> import eepower as eep # Import EEPOWER
>>> real = 100 # W
>>> reactive = 85 # VAR
>>> eep.powertriangle(P=real,Q=reactive,text="Power Triangle Demo",printval=True)
```

![Power Triangle Output](https://raw.githubusercontent.com/engineerjoe440/ELECTRICALPYTHON/master/docs/example_1.png)



#### Help Printout for eepower.powertriangle:


```python
>>> help(eep.powertriangle)

Help on function powertriangle in module eepower:

powertriangle(P=None, Q=None, S=None, PF=None, color='red', text='Power Triangle', printval=False)
    POWERTRIANGLE Function
    
    Purpose:
    --------
    This function is designed to draw a power triangle given
    values for the complex power system.
    
    Required Arguments:
    -------------------
    NONE; a minimum of two of the optional arguments
          must be entered for proper execution.
    
    Optional Arguments:
    -------------------
    P:          Real Power, unitless; default=None
    Q:          Reactive Power, unitless; default=None
    S:          Apparent Power, unitless; default=None
    PF:         Power Factor, unitless, provided as a
                decimal value, lagging is positive,
                leading is negative; default=None
    color:      The color of the power triangle lines;
                default="red"
    text:       The title of the power triangle plot,
                default="Power Triangle"
    printval:   Control argument to allow the numeric
                values to be printed on the plot,
                default="False"
    
    Returns:
    --------
    NONE;   plots generated.

``` 


