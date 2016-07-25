Matlab function package to "Algorithms of the Ising Model" 
------------------------------------------------------------
		by Mario Ullrich	(28.05.2010)


This package contains all programs that are used for the simulations 
	in the article.
	I hope the names are self-explanatory to see, what the programs do.


In Matlab you could read the instructions, if you change your path 
	to the directory where you saved the files and type

		help "file name"

	or you look in the file and read the documentation.


The directory ".../SpontanMagnetization/" and its subdirectories must not 
	be renamed if you want to use the program "SpontanMagnetization.mat"
	
	
In the directory ".../PerfectSampling/" are some generated random cluster 
	states that are distributed exactly according to the stationary 
	distribution.
	They are of the form:  PW_seed-'SIZE'_'BETA'_'TIME'.mat,  where 
	'SIZE' is the size of the lattice, 'BETA' is the inverse temperature and
	'TIME' is the number of steps that we had to go in the past to finish.
	
	With this states you can (after loading in Matlab) generate other exactly 
	distributed initial Ising states for SpontanMagnetization.mat by 
	using the RcToIsing.m routine.
	
	
	
I would be deeply grateful for any kind of suggestions for improvement.

For questions or suggestions, send a email to

	mario.ullrich@uni-jena.de