# mpSchwarz
Code for experiments with multiprecision versions of the classical Schwarz methods.

-------------------------
The code can be run by one of the three main scripts:
- run_mpSchwarz_NonSym.m
- run_mpSchwarz_Sym.m
- run_mpSchwarz_SymDiagDom.m    
    
each of which will take some time to run (tens of seconds to overnight run), depending on 
- the used software
- the used precisions to run
- the used mesh-size
- the number of problems we run the code for
- the number of methods we run the code for

The outputs are saved as a cell array together with brief explanations of what each output is.

The outputs can then be loaded and plotted using the scripts
- TexPlot_mpSM_Direct_AllPrblms.m
- TexPlot_mpSM_Precs_AllPrblms.m

where (some of) the specifications used for generating the data need to eb specified to load this data for plotting.
