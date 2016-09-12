README

Before running the "modelC.py" script you must do the following:

1.  Download the model file from febio.org.  The file is
called "hip_n10rb.feb" and is available in the example problems.
You can access the example problems as follows:
- Log in to febio.org
- From the menu bar select "FEBio" --> "FEBio Downloads".
- Click on the "Example Problems" link under "Test Suites"
*Note: we are unable to redistribute the FEB model file here
due to license restrictions.


2.  Add output requests to the FEB file as follows: 
- Open "hip_n10rb.feb" in any text or XML editor.
- Find the <Output> tag near the bottom of the file.
- Adjust the <logfile> contents to output all element stesses
and strains as follows:

<logfile>
    <element_data data="Ex;Ey;Ez;Exy;Eyz;Exz;sx;sy;sz;sxy;syz;sxz">1:149300</element_data>
</logfile>


3.  In the ### USER VARIABLES ### section of "modelC.py" set the
"path2febio" variable to the path of the "FEBio2" executable on your
system.
!!! NOTE !!!
The example file "hip_n10rb.feb" exits with errors if using FEBio 
version 2.5.0. An earlier version of FEBio is needed.
Suggested:  version 2.4.2.


4. In the ### USER VARIABLES ### section of "modelC.py" set the
following variables:
- "fnameFEB"
- "fnameTEMP" 
- "K" material parameter value

