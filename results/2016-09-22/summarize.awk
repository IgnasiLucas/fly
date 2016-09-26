function InbreedingCoefficient(S,  F0, F, LOCI){
   # S is a string of 0s and 1s, indicating whether the two alleles at a locus
   # are identical by descent (0) or not (1).
   F0 = length(S)
   F = F0
   split(S, LOCI, "")
   for (locus in LOCI) {
      F = F - LOCI[locus]
   }
   return F / F0
}
function CountTracts(S,R,  N,L,T){
   # S is the string of 0s and 1s. R is an array, that will be modified to write
   # in order, the number of tracts of 1s, their mean length, and their individual lengths.
   split(S, A, /0+/)
   for (a in A) {
      if (length(A[a]) > 0) {
         N++
         R[N + 2] = length(A[a])
         T += length(A[a])
      }
   }
   R[1] = N
   R[2] = 0
   if (N > 0) R[2] = T/N
}
BEGIN{
   print "Gen\tF1\tF2\tFMean\tN.Tracts.1\tN.Tracts.2\tMean.Length.1\tMean.Length.2\tLengths.1\tLengths.2"
}(/^Gen/){
   split($1,A,/:/)
   GEN = A[2]
   F1  = InbreedingCoefficient($2)
   F2  = InbreedingCoefficient($3)
   delete(TRACTS1); TRACTS1[1] = 0
   delete(TRACTS2); TRACTS2[1] = 0
   CountTracts($2, TRACTS1)
   CountTracts($3, TRACTS2)
   LIST1 = "-"; LIST2 = "-"
   if (TRACTS1[1] > 0) {
      LIST1 = TRACTS1[3]
      for (i = 4; i <= length(TRACTS1); i++) LIST1 = LIST1 "," TRACTS1[i]
   }
   if (TRACTS2[1] > 0) {
      LIST2 = TRACTS2[3]
      for (i = 4; i <= length(TRACTS2); i++) LIST2 = LIST2 "," TRACTS2[i]
   }
   printf "%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%.2f\t%.2f\t%s\t%s\n", GEN, F1, F2, (F1 + F2)/2, TRACTS1[1], TRACTS2[1], TRACTS1[2], TRACTS2[2], LIST1, LIST2
}

