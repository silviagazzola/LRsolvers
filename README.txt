  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHOR:

       Silvia Gazzola
       University of Bath, United Kingdom
       Email: s.gazzola@bath.ac.uk


   REFERENCES:

       Silvia Gazzola, Chang Meng, and James Nagy
       Krylov Methods for Low-Rank Regularization
       arXiv:1910.10664 [math.NA]
       2019

       Silvia Gazzola, Per Christian Hansen, and James Nagy	
       IR Tools: A MATLAB Package of Iterative Regularization
       Methods and Large-Scale Test Problems
       NUMERICAL ALGORITHMS, 81(3), pp. 773-811, 2019


   SOFTWARE LANGUAGE:

       MATLAB 9.3 (R2017b)


=====================================================================
SOFTWARE
=====================================================================

The present MATLAB software codes some of the new available low-rank solvers. Should be used in connection with IR Tools (whose basic functionalities are provided within the present repository). 

The software has been developed and tested using MATLAB version 9.3.
No other MathWorks products or toolboxes are required.



=====================================================================
HOW TO INSTALL AND CHECK THE INSTALLATION
=====================================================================

Please follow these steps:

- Extract the archive in a directory of your choice. This creates a
  directory called IRTools with all files of the toolbox.

- Start MATLAB.

- Add the directory IRTools containing the toolbox as well as all
  subdirectories to the MATLAB search path.
  - The easiest way to do this is by using the function
    IRtools_setup which updates the path permanently.
  - The alternative is to use the "Set Path" tool in the "Home" tab
    of the main MATLAB window, by clicking "Add with Subfolders..."
    and choosing the IRTools directory.

=====================================================================
GENERAL ORGANIZATION OF THE REPOSITORY
=====================================================================

Some of the new low-rank (LR) Krylov solvers based on the Arnoldi algorithm are available in the `LRsolvers' folder. The other folders contain essential IR Tools functionalities, that are needed to run the new LR solvers and to perform comparisons. The MATLAB script `Example1' reproduces the result of the first numerical test described in the arXiv document arXiv:1910.10664 [math.NA].
