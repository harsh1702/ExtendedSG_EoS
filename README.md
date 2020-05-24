# ExtendedSG_EoS
Finding the extended stiffened gas EoS parameters for compressible diphasic flow under mechanical equilibrium. 
This code will find the Stiffened gas EoS parameters according to the paper
Alexandre Chiapolino, Richard Saurel. Extended Noble–Abel Stiffened-Gas Equation of State for
Sub-and-Supercritical Liquid-Gas Systems Far from the Critical Point. Fluids, MDPI, 2018, 3 (3),
pp.48. 10.3390/fluids3030048.hal-01836245<br/>
Works for when the vapour is not considered an ideal gas the values for pInf for it are non-zero the given initial test case will not work for water vapour as it is approximately an ideal gas water was used as the test case to verify the code as the paper by Chiapolino et al aldready had given results for water.

Compile by using<br/>
$make

Input files<br/>
There are five input files in the input folder namely 'expData.txt', 'refStateAtm.txt', 'refStateCrit.txt', 'refStateLiq.txt' and 'refStateVap.txt'<br/>
'expData.txt' contains the experimental data for the liquid and vapour as given in the format of the aldeready existing file<br/>
'refStateAtm.txt' contains the atmospheric reference state data used to compute the value of CL and CG<br/>
'refStateCrit.txt' conatains the required values of the substance at its critical point<br/>
'refStateLiq.txt' conatins the values of the liquid parameters at the choosen liquid refernce state<br/>
'refStateVap.txt' conatins the values of the vapour parameters at the choosen vapour refernce state<br/>

Execute the program by using <br/>
$./exe

Plotting results<br/>
In the 'res' folder execute the command given below for plots of theoretical vs experimental inputs<br/>
$./runPlot.sh

Things to add:-
1)Calculation of liquid and gaseous entropies<br/>

