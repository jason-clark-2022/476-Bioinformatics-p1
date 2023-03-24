COMPILE AND RUN:
g++ 'Sequence Alignment.cpp'
./a.out


COMPILE AND RUN ALT:
g++ 'Sequence Alignment.cpp' -fconcepts
./a.out


COMPILATION NOTES
With regards to file compilation, there is only one file that needs to be compiled
which is the 'Sequence Alignment.cpp' file, however due to the use of autos as
return types for some functions, it will throw warnings, though the program should 
still function (at least from my experience). The addition of the -fconcepts 
command seems to handle those warnings.

RUNTIME NOTES
With regards to the program itself, the outputs for some files / submatrix variations
seem to vary when compared to the the outputs of the emboss software. I am not sure
exactly where the discrepency is comming from, but the assumption that I made is that
there is variance between the substitution matricies being used between the program
that I made, and the emboss software. I would like to note that the discrepencies
also only seem to arise when using polypeptides sequences, for the nucleotide
sequences appear to run without issue. 


