---
title: 'KAT: A Python-based GUI for the analysis of enzyme kinetic data using both classical and complex models. '
tags:
  - Python
  - biochemisty
  - enzyme kinetics
  - kinetics
  - fluorescent data
  - Michaelis-Menten Equation
  - Hill Equation
authors:
  - name: Evan R. Jones
    orcid: 0009-0008-2914-8558
    affiliation: 1
affiliations:
 - name: Department of Chemistry and Biochemistry, University of Oklahoma, United States
   index: 1
   ror: 02aqsxs83
date: 20 July 2025
bibliography: paper.bib

---

# Summary

When studying previously unknown or uncharacterized enzymes, kinetics serves as the basis
for understanding the key properties of the enzyme. These properties can give 
rise to vital information regarding binding affinity of substrate, energetics or 
speed of the reaction, and even potential allosteric effects within the enzyme 
itself. These parameters often help guide novel therapeutic discovery.
Often, kinetic assays are the simplest way to begin to ascertain key
details about the function of a given enzyme, allowing properties to be compared 
between wild-type and mutant enzymes. Currently, kinetic data is processed in a 
concerted, step-wise fashion that allows it to be swiftly automated and fitted 
to a number of kinetics equations, a necessary step for understanding novel enzymes. 
A fully-functional GUI that accepts raw data and outputs both kinetic parameters 
and a graph of the fit not only significantly accelerates data analysis, but also 
guarantees the same treatment of data across replicas and experiments.


# Statement of need

`KAT` is a toolkit dedicated to parsing both fluorescent and absorbance kinetic 
data and fitting the data to several classical (Michaelis-Menten and Hill)
[@johnson:2011; @goutelle:2008] and complex (Monod-Wyman-Changeux and Koshland-Nemethy-Filmer)
[@monod:1965; @koshland:1966] models. Each of these models can also be used to analyze either
replicate data or data arising from several mutations. These replica fitting
functions output either averages and standard deviations for each kinetic 
parameter or side-by-side analysis for mutant data sets. The Python-based GUI
allows for simple inputting of necessary data: the CSV file with fluorescent or
absorbance data, the substrate information (number of substrate concentrations, 
dilution factor, and maximum concentration), and the window within which to 
calculate the velocity data. There is a built-in function to auto calculate the
linear range of the assay data, which is a necessary parameter for many kinetic
models. This is calculated through analysis of subsequent slope values and only
beginning data analysis once the difference in slopes is less than 5% on average.
Once the user has provided the necessary data, s/he can easily fit the data to 
multiple models, allowing for the comparison of fits and kinetic parameters.

While many other tools effectively fit Michaelis-Menten kinetics to a data set, 
fewer fit data to the non-linear form of the Hill Equation, which adds the extra
Hill coefficient parameter. Further, software like EnzFitter do not have complex
models built-in, requiring the user to input complex algebraic equations by hand.
[@leatherbarrow:1988] `KAT` already has these models integrated and directly outputs 
easily-modifiable SVG graphs, as well as typical PNG-formatted graphs. For the 
classical models, the kinetic parameters are conviently displayed upon completing 
the calculation, whereas the complex models display each parameter in a pop-up window.

`KAT` utilizes both numerical solving (classical models) and optimization (complex
models) using the lowest residual square sum for the determination of kinetic 
parameters. For fitting to classical models, a standard 12 substrate concetrations 
is enough to provide substantial confidence in the model fitting; however, the use 
of complex models that solve for up to 6 parameters are difficult to fit when using 
under ~30 substrate concentrations. Therefore, several statistical techiniques have 
been implemented within `KAT` to test the confidence interval of the model fit and
each parameters. Should one parameter fall without these intervals, a warning message 
is displayed in the GUI citing potential poor confidence. 

With an efficient method of analyzing raw enzyme kinetic data using a simple GUI,
`KAT` will undoubtedly be useful for a wide-range of enzyme types and significantly
standardize the fitting of data to both classical and complex kinetic models.

# Mathematical Basis
## Kinetic Models
### *Classical*
Based on the simple model of an enzyme existing in three discreet states, the Michaelis-Menten 
equation identifies a *V~max~* and a *K~M~* that define the maximum velocity and 
the Michaelis constant, which correlates to binding affinity between substrate and 
enzyme, respectively.[@johnson:2011] The Michaelis-Menten equation assumes that the three 
states are discreet and independent. 

\begin{equation}
E + S \longleftrightarrow ES \longleftrightarrow E + P
\end{equation}

After rearranging for the rate-constants that govern each transition, one solves for the 
complete Michaelis-Menten equation, where [S] is the substrate concentration.

\begin{equation}
v = \frac{V_{max}[S]}{K_{M}+[S]}
\end{equation}

The Hill equation is an extension of the Michaelis-Menten equation that accounts for
positive coopertivity, which is often the result of allostery within the enzyme and
is the effect where binding of one substrate molecule accelerates the binding of another.
The Hill coefficient, or n, is a measure of the coopertivity, where values greater than 1
exhibit positive coorpertivity (*i.e.* increased binding of subsequent ligands) and values
less than 1 exhibit negative coopertivity (*i.e* decreased binding of subsequent ligands).
[@goutelle:2008]

\begin{equation}
v = \frac{V_{max}[S]^n}{K_{M}^n+[S]^n}
\end{equation}

### *Complex*

While the Michaelis-Menten and Hill  models often explain the majority of enzymes, the Hill 
equation in particular is limited in its explanation for the mechanism of allostery within 
the enzyme. Monod, Wyman, and Changeux further expanded upon the equation and introduced a 
two-stage existance of an allosterically-regulated enzyme: a tensed or "T" state, 
where binding of substrate and catalysis is limited, and a relaxed or "R" state, where both
binding and catalysis are accelerated. These additional states each have a *V* and a *K* 
parameter that govern velocity of catalyis and affinity, respectively. Further, there are
two additional parameters (for a total of 6): *L~0~*, or the cooperativity coefficient as a ratio
of the enzyme in the T state vs. the R state, and *N*, or typically the number of allosteric
sites in the enzyme.[@monod:1965] These parameters combine to give the following 
Monod-Wyman-Changeux equation:

\begin{equation}
v = \frac{V_{T}L_{0}(1 + \frac{[S]}{K_{T}}) + V_{R}(1 + \frac{[S]}{K_{R}})}{L_{0}(1 + \frac{[S]}{K_{T}})(1 + \frac{[S]}{K_{R}})}
\end{equation}

Although the Monod-Wyman-Changeux model accurately predicts the transitions between the T 
and R states, the model requires that the T and R states exist in discreet environments, 
where part of the enzyme cannot exist in both states at once. Since not all enzymes exhibit
a complete transition from T to R states at once, Koshland, Nemethy, and Filmer developed a 
model that accounts for enzymes that exist in states T and R simultaneously by introducing a 
$\gamma$ term. This term provides the basis for coopertivity in the enzyme. Additionally,
the Koshland-Nemethy-Filmer model utilizes a term for each acitve site of the enzyme (denoted
as *i* here).[@koshland:1966] Upon summing all possible enzyme states (from unbound to bound) 
in each active site, *i*:

\begin{equation}
v = \frac{E_{total}\sum_{i=0}^{j} k_{i} \binom{j}{i} (\frac{[S]}{K_{d}})^i}{\sum_{i=0}^{j} \binom{j}{i} (\frac{[S]}{K_{d}})^i}
\end{equation}

where *j* is the total number of active sites, *E~total~* can be built into a *V~max~* term,
and $k_{i} = k_{basal} + (V_{max} - k_{basal})(\frac{i}{j})^\gamma$. Due to the large number 
of parameters, these complex equations are difficult to assess confidence of fit given the 
estimated starting guess ("Best-Fit") when the number of substrate concentrations is less than 30.
Therefore, Cross-Validation of the solved parameters using the KFold technique with 10 splits 
is used. Further, if the number of substrate concentrations is below 30, Bayesian bootstrapping 
is implemented to assess the 99% confidence intervals of each parameter, and the "best-fit" 
data is tested to be within these confidence intervals. If a best-fit parameter falls outside 
of the 99% confidence interval, the cross-validation parameters are provided instead of "best-fit." 
Otherwise, the "best-fit" values are provided.

## Model Fitting
### *Using Classical Models*

For fitting data to a classical model, a simple minimization of residual sum of squares is 
performed using a data-driven starting guess. This starting guess uses an average of the
three largest velocity values for *V~max~* and the substrate concentration at half *V~max~*
as a guess for *K~M~*. If the Hill model is used, a Hill coefficient of 2 is used as the
starting guess for *n*. `sympy.nsolve` is used to numerically solve the three partial 
derivatives, $\frac{\partial Q}{\partial V_{max}}$, $\frac{\partial Q}{\partial K_{M}}$,
$\frac{\partial Q}{\partial n}$, of the following equation:[@meurer:2017]

\begin{equation}
Q = \sum_{s=1}^{S} (V_{data} - \frac{V_{max}s}{K_{M}+s})^2
\end{equation}

for the Micaelis-Menten model, or 

\begin{equation}
Q = \sum_{s=1}^{S} (V_{data} - \frac{V_{max}s^n}{K_{M}^n+s^n})^2
\end{equation}

for the Hill model. The relatively few parameters when compared to the number of typical
substrate concentrations and good, data-specific guess allows for high confidence in the 
initial fit

### *Using Complex Models*

Unlike fitting the classical kinetics models, the complex models have too many parameters
to effiently and accurately numerically solve the partial derivatives of the extended 
residual sums. Therefore, `scipy.optimize.minimize` has been implemented to minimize the
loss function (defined similarly as above in (6) and (7)). The minimization is bound by 
typical biological constraints, where *V~T~*, *V~R~*, *K~T~*, and *K~R~* and bound between
0 and 1,000, *L~0~* is bound between 0.001 and 500, and *N* is bound between 0.5 and 14 for
the Monod-Wyman-Changeux equation. Similarly for the Koshland-Nemethy-Filmer model, *V~max~*,
*K~d~*, and *k~basal~* are bound between 0 and 10,000 and $\gamma$ is bound between -50 and 50.
Due to the unstable structure of fitting relatively few data points to complex models, 
`scipy.optimize.differential_evolution` is first used to help identify global minima.
[@virtanen:2020]

# Acknowledgements

I acknowledge contributions from Dr. Christina Bourne and the University of Oklahoma Department
of Chemistry and Biochemistry for their resources.

# References
