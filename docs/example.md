---

[MAIN](https://engineerjoe440.github.io/electrical-engineering-python/index)
 | 
[Example Usage](https://engineerjoe440.github.io/electrical-engineering-python/example)

---

# ELECTRICALPYTHON
*The Python Package for all Electrical Engineers!*

### Example Usage:
The following presents a simple usage of the *eepower.powertriangle* function.

#### Function Usage:
The *powertriangle* function will accept two arguments to generate a Power-Triangle. In addition to it's primary goal of generating and plotting the triangle, the function can also return the remaining unknowns and print all values on the plot.

#### Required Function Arguments (any two of the following):
- P: Real Power (Watts, absolute value - real, not complex)
- Q: Reactive Power (VARs, absolute value - real, not complex)
- S: Apparent Power (VAs, absolute value - real, not complex)
- PF: Power Factor (expressed as a decimal value - i.e. 0.8)

#### Optional Function Arguments:
- color: Line color of the plotted triangle, default is red.
- text: Text to appear as the plot title, default is "".
- figure: Figure number, soon to be depreciated, default is 1.
- printval: Argument specifying whether the values should be printed on the plot, default is False.
- ret: Argument specifying whether the values should be returned after the function call, default is False.
- plot: Argument specifying whether the plot should be generated, default is True.

#### Function Demonstration:


```python
import eepower as eep # Import EEPOWER
real = 100 # W
reactive = 85 # VAR
eep.powertriangle(P=real,Q=reactive,text=" Demo",printval=True)
```


![Power Triangle Output](https://raw.githubusercontent.com/engineerjoe440/electrical-engineering-python/master/docs/example_1.png)


#### Help Printout for eepower.powertriangle:


```python
print(help(eep.powertriangle))
```

    Help on function powertriangle in module eepower.eepower:
    
    powertriangle(P=False, Q=False, S=False, PF=False, color='red', text='', figure=1, printval=False, ret=False, plot=True)
        ###################################################################
        #   Define Power Triangle Function
        #
        #   Draws a Power Triangle given a set of two values from set of:
        #   { P, Q, S, PF }. Can also return all values from said set.
        #   Plot title can be added to, and plot may be turned off.
        #   Color of plot may be changed, and multiple figures may be shown.
        #
        #   Requires two values from set: { P, Q, S, PF }
        #   All values given must be given as absolute value, not complex.
        ###################################################################
    
    None
    


```python

```
