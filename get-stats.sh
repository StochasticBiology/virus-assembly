# hacky bash script to extract summary statistics from a set of simulation outputs in the same directory

# use tail to get the last line of each summary file along with the file label, then processes this pair to produce rows containing simulation parameters and summary output (latest time, latest energy, latest count of 12-mers)
tail -n1 virusout*sizes* | awk '{if($1 == "==>") printf("%s ", $2); else if($1 > 0) printf("%i %f %i\n", $1, $2, $15);}' | sed 's/virusout-//g' | sed 's/[.]sizes.txt//g' | sed 's/\([0-9]\)-/\1 /g' > sim-summary.txt

# construct a Gnuplot script plotting the energy trajectory with time for each simulation summary found in this directory
ls virusout*sizes* | awk 'BEGIN{printf("reset\nset xlabel \"Cycle\"\nset ylabel \"Energy\"\nset key outside right\nset term svg size 1600,480\nset output \"plot-energies.svg\"\nplot ");}{printf("\"%s\" u 1:2 w l, ", $0);}' | sed 's/, $/\n/g' > plot-energies.sh

# run that script
gnuplot "plot-energies.sh"


