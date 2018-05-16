#!/bin/bash
number=10000
i=0
while [ "$i" -le "$number" ]; do
#
# cleaning
#
echo ""
echo "Cleaning out old output files..."
rm -f iso_massf*DAT
rm -f x-time.dat
rm -f *png
rm -f screen.out
#
# preparing new run
#
echo "Preparing new run..."
network=networksetup_$i.txt
# copy from a directory where you have created imc=1,nmc networksetup_imc.txt files
# to your ppn run directory
cp /Users/dpa/NuGrid_github/NuPPN/examples/prep_onezone_MC_sims/$network ./networksetup.txt
echo $network
echo ""
echo "==========> RAWD_ppn run $i"
echo ""
./ppn.exe > screen.out
wait $!
printf -v i0 "%05d" $i
# the abundances from the cycle 891 of your ppn run are stored in the directory
# onezone_MC_sims_results
mv iso_massf00891.DAT onezone_MC_sims_results/onezone_MC_sims_iso_massf$i0.DAT
echo onezone_MC_sims_iso_massf$i0.DAT
i=$(($i + 1))
done
