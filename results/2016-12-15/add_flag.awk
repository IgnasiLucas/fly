BEGIN{
   OFS = "\t"
   WRITTEN = 0
}
(/^#/){
   if ((/^##FORMAT/) && (WRITTEN == 0)){
      print "##INFO=<ID=BPF,Number=A,Type=Integer,Description=\"Binary presence flag, indicates what samples each allele is present in.\">"
      WRITTEN = 1
   }
   print
}(/^[^#]/){
   delete FLAG
   FLAG[0] = 0
   split($5, ALLELES, /,/)
   ALLELES[0] = $4
   GT = 0
   split($9,ORDER,/:/)
   for (o in ORDER) {
      if (ORDER[o] == "GT") GT = o
   }
   for (i = 0; i <= (NF - 10); i++) {
      split($(i+10), FORMAT, /:/)
      for (a in ALLELES) {
         if (FORMAT[GT] ~ a) FLAG[a] += 2^i
      }
   }
   BPF = ";BPF=" FLAG[0]
   for (a = 1; a < length(ALLELES); a++) {
      BPF = BPF "," FLAG[a]
   }
   $8 = $8 BPF
   print $0
}
