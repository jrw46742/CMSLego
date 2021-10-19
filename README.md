# CMSLego

![Lego](https://i.gyazo.com/d4d037e2fd4f4bfc06e9b2fb53c3315e.png)

Use [this document](https://docs.google.com/document/d/10PcXfrXLM_u_z3NdOkStSNnz8gG4Fx1YDnN6xqQEulE/edit?usp=sharing) as a guide through [this spreadsheet](https://docs.google.com/spreadsheets/d/1fynQMzaGjha22tBcTz--OduF70qkugXDlZUlLmLmJWQ/edit?usp=sharing).

- [HitstoLego.py](https://github.com/jrw46742/CMSLego/blob/main/HitstoLego.py):
  - Also Located in ~jwesterl/nobackup/LegoPflow/CMSSW_10_2_22/src/HitstoLego.py
  - For ECAL hits:
    - From Untime.cc
    - Printout ECAL hits in order iEta,iPhi,hitE and save into text file ehits.txt in the same directory as HitstoLego.py
    - Change boolean flag isECAL to True
  - For HCAL hits:
    - From Untime.cc
    - Printout HCAL hits in same order and save into hhits.txt
    - Change boolean flag isECAL to False
  - $ python HitstoLego.py
    - Outputs 2d array of your input file (for sanity checks)
    - Outputs # of hits
    - Outputs # of blue/red lego's you will need (accounts for HCAL hits being 5x5)
    - Outputs ECAL.csv or HCAL.csv
  - Open ECAL.csv or HCAL.csv in a spreadsheet
  - Apply conditional formatting based on # of legos
  - Which you can set for yourself in the colorgradiant() function
