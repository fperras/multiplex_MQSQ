# multiplex_MQSQ
Pulse sequence and processing program for acquiring multiplex MQSQ spectra

Multiplex Processing Program

This program is used for the processing of multiplex MQ/SQ Bruker NMR data 
acquired using the R12v5n2_multiplex.scs pulse program. The data are
reorganized into a 3D spectrum with the MQ/SQ (M=0,1,2...) 2D spectra along the
F3-F1 planes, starting with the 0Q/SQ spectrum.
-written by S. Chandra Shekar, 01/19/2023

The program can be compiled as: g++ multiplex_processing.cpp -lm

To convert the R12v5n2_multiplex.scs 3D spectrum the following steps are
to be followed:

(1) Copy the dataset into a new directory

(2) Paste the executable for the Multiplex Processing Program into the newly
    created directory and run it. 

(3) The program will then request some information on the original acquisition
    (TD1, TD2, TD3, file names) which must be entered as directed. The program
    will read the original ser file (usually named 'ser') and generate a new ser
    file (with a name you provide) with the processed MQ/SQ 3D data. The program
    will write two new files, namely the new processed ser file and a log file
    that contains additional information for reading the data in Bruker Topspin

(4) Open the created log.txt file. Using Bruker Topspin update the values of 
    '1 td', '2 td', 's 1 td', and 's 2 td' to the TD1 and TD2 values in the log
    file.

(5) In the Topspin 'edp' window, set the following processing values:

               F3    F2    F1
    phmod      pk    no    pk            (repeat for "s phmod")
    ftmod     fqc    no    no            (repeat for "s ftmod")
    mc2              QF    echo-antiecho (repeat for "s mc2")

(6) Fourier transform the data with the commands 'tf3' and 'tf1' or with 'xfb'
    to process a particular F1-F3 slice.

(7) The F1 referencing will be different in each of the MQ/SQ datasets. One way
    to ensure a correct referencing is to reference the center of the spectrum
    to M-times the value of the chemical shift in the center of F3.
