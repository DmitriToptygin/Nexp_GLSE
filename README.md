# Nexp_GLSE
Nexp_GLSE is a Global Least Square Estimator for fitting lifetime or dwell time distributions by a linear combination of N exponential terms. The experimental setup is assumed to have an infinitely narrow Impulse Response Function (IRF), therefore no convolution is carried out. Other programs by the same author are designed specifically for the experimental data with wide IRF, for example, for Time Correlated Single Photon Counting (TCSPC) data; there convolution is included in the model function. The weighting of squared residuals is based on the assumption of Poissonian distribution of the counts in the bins. The program can handle bins of equal or unequal widths. The program can analyze multiple data sets simultaneously, with some of the model parameters being shared between two or more data sets (global parameters), while the other model parameters being relevant to only one data set (local parameters). The bin widths can be different for different data sets. Minimization of weighted least squares is performed using a version of Levenberg-Marquardt algorithm.

TO INSTALL AND TEST THE PROGRAM FOLLOW THESE STEPS:

1. Compile the source code using the "make" build automation tool in the current directory, which contains Makefile. Use the ls command to make sure that the Makefile is present. Then use the make command with no command-line options.

2. Move the executable files Nexp_GLSE and Nexp_dwelltime_simulator to a directory that is contained in the PATH. To print a list of all directories contained in the PATH use the command echo $PATH. If you are not the owner of the destination directory, then you will need to use the sudo command, for example:

          sudo mv Nexp_GLSE Nexp_dwelltime_simulator /usr/local/bin

3. Before using the program(s) please read the user manual, it is in the file ./doc/Nexp_GLSE_manual.pdf . Answers to most questions can be found in that pdf document.

4. To test the operation of the program Nexp_GLSE on a single data set, change to the directory test_single and run the program:

          cd test_single
          Nexp_GLSE < Nexp_GLSE_test.inp

5. To test the operation of the program Nexp_GLSE for global analysis, change to the directory test_global and run the program:

          cd test_global
          Nexp_GLSE < Nexp_GLSE_test.inp
